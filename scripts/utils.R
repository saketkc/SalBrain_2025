suppressPackageStartupMessages({
  library(Seurat)
  library(cluster)
  library(tidyverse)
  library(diffusionMap)
  library(princurve)
  library(ggpubr)
  library(FNN)
  library(grDevices)
  library(RColorBrewer)
  library(tradeSeq)
  library(patchwork)
  library(pheatmap)
  library(sparseMatrixStats)
})
theme_set(theme_pubr())

stages <- paste0("ST", c("30", "33", "36", "41", "46", "50", "53", "55"))
Spectral <- colorRampPalette(rev(brewer.pal(length(stages), "Spectral")))
L1colors <- list(MP = "#e41a1c", DP = "#377eb8", LP = "#4daf4a", VP = "#984ea3", AMY = "#ff7f00", Unknown = "#969696")

KnnSmooth <- function(y, coords, k) {
  knn.out <- get.knn(data = coords, k = k)
  w <- 1 / (knn.out$nn.dist + .Machine$double.eps)
  w <- w / apply(X = w, MARGIN = 1, FUN = sum)
  v <- apply(X = knn.out$nn.index, MARGIN = 2, FUN = function(i) y[i])
  return(apply(v * w, 1, sum))
}


GetPCAVar <- function(object) {
  eigs <- (object@reductions$pca@stdev)^2
  df.pervar <- data.frame(
    variance = eigs,
    percent.var = eigs * 100 / sum(eigs),
    PCs = 1:length(eigs)
  )
  return(df.pervar)
}


GetVaryingGenes <- function(object, meta.data, tech.covariate = NULL, genes.to.test = NULL) {
  pseudotime <- data.frame(Curve1 = meta.data$maturation.score.smooth, row.names = rownames(meta.data))
  cellWeights <- pseudotime
  cellWeights$Curve1 <- 1.0

  counts <- object@assays$RNA@counts
  if (!is.null(genes.to.test)) {
    counts <- counts[intersect(rownames(counts), genes.to.test), ]
  }
  if (!is.null(tech.covariate)) {
    tech.model.matrix <- model.matrix(as.formula(paste0("~", paste0(tech.covariate, collapse = "+"))), data = meta.data)
  } else {
    tech.model.matrix <- NULL
  }



  z <- BiocParallel::bpparam()
  z$workers <- 42
  gam.results <- fitGAM(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    U = tech.model.matrix,
    parallel = TRUE,
    BPPARAM = z,
    nknots = 3,
    verbose = TRUE
  )

  # which genes change along this trajectory
  startRes <- startVsEndTest(gam.results)
  startRes$gene <- rownames(startRes)
  startRes.shortlist <- startRes %>%
    filter(pvalue < 0.01) %>%
    arrange(desc(waldStat))
  return(startRes.shortlist)
}


SmoothSpline <- function(gene.data, time.data, df = 3) {
  spl <- smooth.spline(x = time.data, y = gene.data, df = df, tol = 1e-6)
  smoothed_y <- predict(spl, time.data)$y
  return(smoothed_y)
}

SmoothSpline.Matrix <- function(time.data, gene.matrix, df = 3) {
  z <- apply(X = gene.matrix, MARGIN = 1, FUN = function(row) SmoothSpline(row, time.data))
  z <- t(z)
  rownames(z) <- rownames(gene.matrix)
  colnames(z) <- colnames(gene.matrix)
  return(z)
}



GetFASCT <- function(object) {
  fa <- object@assays$SCT@SCTModel.list[[1]]@feature.attributes
  fa$gene <- rownames(fa)
  var.genes.postregress <- rowVars(x = object@assays$SCT@scale.data)
  var.genes.postregress.df <- data.frame(
    residual_variance_postregress = var.genes.postregress,
    gene = rownames(object@assays$SCT@scale.data)
  )
  fa2 <- left_join(fa, var.genes.postregress.df)
  fa2 <- fa2 %>% arrange(desc(residual_variance), desc(residual_variance_postregress))
  return(fa2)
}

cell.cycle.genes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

GetPercentRibo <- function(object) {
  ribo.genes <- grep(
    pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", x = rownames(object@assays$RNA@counts),
    ignore.case = T, value = T
  )
  object[["percent.ribo"]] <- PercentageFeatureSet(
    object = object,
    features = ribo.genes
    # pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^rp[sl][[:digit:]]|^rplp[[:digit:]]|^rpsa"
  )
  return(object)
}

DoDiffusion <- function(object, npcs) {
  object <- object %>%
    FindNeighbors(dims = 1:npcs, verbose = FALSE) %>%
    RunUMAP(dims = 1:npcs, verbose = FALSE) %>%
    CreateDiffusionMap(npcs = npcs)
  return(object)
}

DoFeatureSelection <- function(object, var.threshold = 1.5, cell.cycle.genes = NULL, ribo.genes = NULL) {
  fa <- GetFASCT(object)
  fa$is_cell_cycle <- FALSE
  fa$is_cell_cycle[fa$gene %in% cell.cycle.genes] <- TRUE
  fa$is_ribo <- FALSE
  fa$is_ribo[fa$gene %in% ribo.genes] <- TRUE
  fa.shortlist <- fa %>%
    filter(residual_variance_postregress >= var.threshold) %>%
    filter(!is_cell_cycle) %>%
    filter(!is_ribo)
  VariableFeatures(object = object) <- fa.shortlist$gene
  return(object)
}
DetermineNPCs <- function(object, var.thresh = 2.5) {
  npcs <- GetPCAVar(object) %>%
    filter(percent.var >= var.thresh) %>%
    nrow()
  return(npcs)
}


ProcessObject <- function(object, pca.var.thresh, align_PCs_to_gene = NULL, do.diffusion = FALSE) {
  object.subset <- DoFeatureSelection(object,
    cell.cycle.genes = cell.cycle.genes,
    var.threshold = 1
  )
  object.subset <- GetResidual(
    object = object.subset,
    features = VariableFeatures(object.subset)
  )
  object.subset <- RunPCA(object = object.subset)
  npcs <- DetermineNPCs(object = object.subset, var.thresh = pca.var.thresh)
  object.subset <- RunPCA(object = object.subset, npcs = npcs)

  if (!is.null(align_PCs_to_gene)) {
    feature.loadings <- object.subset@reductions$pca@feature.loadings
    cell.embeddings <- object.subset@reductions$pca@cell.embeddings
    gene.embedding <- object.subset@reductions$pca@feature.loadings[align_PCs_to_gene, ]
    PCs.negative <- gene.embedding[gene.embedding < 0]

    # for these PCs flip the sign of feature.loadings and gene.embeddings
    for (PC in names(PCs.negative)) {
      print(PC)
      feature.loadings[, PC] <- -1 * feature.loadings[, PC]
      cell.embeddings[, PC] <- -1 * cell.embeddings[, PC]
    }
    colnames(feature.loadings) <- paste0("alignedPCA_", seq(1, ncol(feature.loadings)))
    colnames(cell.embeddings) <- paste0("alignedPCA_", seq(1, ncol(cell.embeddings)))

    allreduc <- object.subset@reductions$pca
    allreduc@key <- "alignedPCA_"
    allreduc@feature.loadings <- feature.loadings
    allreduc@cell.embeddings <- cell.embeddings
    object.subset[["alignedPCA"]] <- allreduc
  }

  if (do.diffusion) {
    object.subset <- DoDiffusion(object = object.subset, npcs = npcs)
  }
  object.subset <- FindNeighbors(object = object.subset, reduction = "pca", dims = 1:npcs)
  object.subset <- FindClusters(object = object.subset, resolution = 0.3)
  object.subset <- RunUMAP(object = object.subset, dims = 1:npcs)
  return(object.subset)
}

LearnPrinCurve <- function(object,
                           reduction.name = "DC",
                           dimnames.to.use = NULL,
                           dim1 = 1,
                           dim2 = 2,
                           plot.density = FALSE,
                           tech.covariate = NULL,
                           pricu.f = 1 / 3,
                           stretch = 2,
                           only.prin = T, genes.to.test = NULL) {
  if (reduction.name == "pca") {
    reduction.name.1 <- paste0("PC_", dim1)
    reduction.name.2 <- paste0("PC_", dim2)
  } else {
    reduction.name.1 <- paste0(reduction.name, "_", dim1)
    reduction.name.2 <- paste0(reduction.name, "_", dim2)
  }

  cell.embeddings <- object@reductions[[reduction.name]]@cell.embeddings
  if (!is.null(dimnames.to.use)) {
    cell.embeddings.forprincurv <- cell.embeddings[, dimnames.to.use]
  } else {
    cell.embeddings.forprincurv <- cell.embeddings
  }

  pricu <- principal_curve(cell.embeddings.forprincurv,
    smoother = "lowess",
    trace = TRUE, f = pricu.f, stretch = stretch
  )
  pc.line <- as.data.frame(pricu$s[order(pricu$lambda), ])

  meta.data <- cbind(object@meta.data, cell.embeddings.forprincurv)
  meta.data$maturation.score <- pricu$lambda / max(pricu$lambda)
  meta.data$maturation.score.smooth <- KnnSmooth(meta.data$maturation.score, cell.embeddings.forprincurv[, 1:2], 20)


  p <- ggplot(meta.data, aes_string(reduction.name.1, reduction.name.2)) +
    geom_point(aes(color = maturation.score.smooth), size = 1, shape = 16) +
    scale_color_gradientn(colours = Spectral(100), name = "Maturation score") +
    geom_line(data = pc.line, color = "red", size = 0.77) #+
  if (plot.density) {
    p <- p + stat_density2d(n = 111, na.rm = TRUE, color = "black", size = 0.33, alpha = 0.5)
  }

  if (only.prin) {
    return(list(meta.data = meta.data, plot = p, pc.line = pc.line))
  }
  pseudotime <- data.frame(Curve1 = meta.data$maturation.score.smooth, row.names = rownames(meta.data))
  cellWeights <- pseudotime
  cellWeights$Curve1 <- 1.0

  counts <- object@assays$RNA@counts
  if (!is.null(genes.to.test)) {
    counts <- counts[intersect(rownames(counts), genes.to.test), ]
  }
  if (!is.null(tech.covariate)) {
    tech.model.matrix <- model.matrix(as.formula(paste0("~", paste0(tech.covariate, collapse = "+"))), data = meta.data)
  } else {
    tech.model.matrix <- NULL
  }



  z <- BiocParallel::bpparam()
  z$workers <- 42
  gam.results <- fitGAM(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    U = tech.model.matrix,
    parallel = TRUE,
    BPPARAM = z,
    nknots = 3,
    verbose = TRUE
  )

  # which genes change along this trajectory
  startRes <- startVsEndTest(gam.results)
  startRes$gene <- rownames(startRes)
  startRes.shortlist <- startRes %>%
    filter(pvalue < 0.01) %>%
    arrange(desc(waldStat))
  return(list(
    meta.data = meta.data, plot = p, startRes = startRes,
    startRes.shortlist = startRes.shortlist, pc.line = pc.line
  ))
}

BreakDownPCA <- function(object, npcs, metadata.axes,
                         title = "Pearson's correlation between PCs and characterized axis of variation",
                         cor.threshold = 0.25,
                         filename = NULL, ...) {
  for (column in colnames(metadata.axes)) {
    if (class(metadata.axes[[column]]) != "numeric") {
      # convert to numeric
      metadata.axes[[column]] <- as.numeric(as.factor(metadata.axes[[column]]))
    }
  }
  correlations.PC <- abs(x = cor(object[["pca"]]@cell.embeddings[, 1:npcs], metadata.axes))
  correlations.PC.round <- round(x = correlations.PC, 2)
  breaksList <- seq(0, 1, by = 0.01)
  color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
  plot <- pheatmap::pheatmap(t(correlations.PC),
    display_numbers = t(correlations.PC.round),
    color = color,
    fontsize = 14, fontsize_number = 13,
    cluster_rows = F,
    cluster_cols = F,
    cellheight = 42,
    cellwidth = 42,
    main = title,
    breaks = breaksList
  )
  if (!is.null(filename)) {
    pheatmap::pheatmap(t(correlations.PC),
      display_numbers = t(correlations.PC.round),
      color = color,
      fontsize = 14, fontsize_number = 13,
      cluster_rows = F,
      cluster_cols = F,
      cellheight = 42,
      cellwidth = 42,
      main = title,
      breaks = breaksList,
      filename = filename, ...
    )
  }
  PCs.selected <- seq(1:nrow(correlations.PC))[colSums(t(correlations.PC) >= cor.threshold) == 0]
  plot2 <- pheatmap::pheatmap(t(correlations.PC[PCs.selected, ]),
    display_numbers = t(correlations.PC.round[PCs.selected, ]),
    color = color,
    fontsize = 14, fontsize_number = 13,
    cluster_rows = F,
    cluster_cols = F,
    cellheight = 42,
    cellwidth = 42,
    breaks = breaksList,
    main = title
  )
  return(list(correlations = correlations.PC, plot = plot, plot2 = plot2, PCs.selected = PCs.selected))
}


GetVaryingGenes2 <- function(object, col.name,
                            tech.covariate = NULL, genes.to.test = NULL, verbose = FALSE) {
  library(tradeSeq)
  meta.data <- object@meta.data
  pseudotime <- data.frame(
    Curve1 = meta.data[[col.name]],
    row.names = rownames(meta.data)
  )
  cellWeights <- pseudotime
  cellWeights$Curve1 <- 1.0

  counts <- object@assays$RNA@counts
  if (!is.null(genes.to.test)) {
    counts <- counts[intersect(rownames(counts), genes.to.test), ]
  }
  if (!is.null(tech.covariate)) {
    tech.model.matrix <- model.matrix(as.formula(paste0("~", paste0(tech.covariate, collapse = "+"))), data = meta.data)
  } else {
    tech.model.matrix <- NULL
  }



  z <- BiocParallel::bpparam()
  z$workers <- 42
  gam.results <- fitGAM(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    U = tech.model.matrix,
    parallel = TRUE,
    BPPARAM = z,
    nknots = 3,
    verbose = verbose
  )

  # which genes change along this trajectory
  startRes <- startVsEndTest(gam.results)
  startRes$gene <- rownames(startRes)
  startRes.shortlist <- startRes %>%
    filter(pvalue < 0.01) %>%
    arrange(desc(waldStat))
  return(list(
    startRes = startRes,
    startRes.shortlist = startRes.shortlist
  ))
}
