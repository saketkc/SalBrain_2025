source("scripts/utils.R")
object <- readRDS("objects/SLCHigh.rds")
pct.expressed <- rowSums(object[["RNA"]]@counts >= 1) %>% as.data.frame()
colnames(pct.expressed) <- "ncounts"
pct.expressed <- pct.expressed %>% arrange(desc(ncounts))
DefaultAssay(object) <- "RNA"
features.to.keep <- pct.expressed %>%
  filter(ncounts > 0) %>%
  rownames()

object <- DietSeurat(object = object, assays = "RNA")
object <- subset(x = object, features = features.to.keep)
object <- SCTransform(
  object = object,
  vst.flavor = "v2",
  vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt", "G2M.Score", "S.Score", "Stage")
) %>% RunPCA()
object <- ProcessObject(object = object, pca.var.thresh = 1, align_PCs_to_gene = "WNT8B")
p.learncurve <- LearnPrinCurve(
  object = object, reduction.name = "alignedPCA",
  tech.covariate = NULL,
  dimnames.to.use = colnames(object@reductions$alignedPCA@cell.embeddings), only.prin = F,
  genes.to.test = VariableFeatures(object),
  pricu.f = 2 / 3,
  stretch = 0
)
p.learncurve

meta.data <- p.learncurve$meta.data

# correlate with WNT8B
if (cor(object@assays$SCT@data["WNT8B", ], meta.data$maturation.score.smooth) > 0) {
  meta.data$maturation.score <- 1 - meta.data$maturation.score
  meta.data$maturation.score.smooth <- 1 - meta.data$maturation.score.smooth

  meta.data$maturation.score <- (meta.data$maturation.score - min(meta.data$maturation.score)) / (max(meta.data$maturation.score) - min(meta.data$maturation.score))
  meta.data$maturation.score.smooth <- (meta.data$maturation.score.smooth - min(meta.data$maturation.score.smooth)) / (max(meta.data$maturation.score.smooth) - min(meta.data$maturation.score.smooth))
}

pc.line <- p.learncurve$pc.line
reduction.name <- "alignedPCA"
reduction.name.1 <- paste0(reduction.name, "_1")
reduction.name.2 <- paste0(reduction.name, "_2")
p <- ggplot(meta.data, aes_string(reduction.name.1, reduction.name.2)) +
  geom_point(aes(color = maturation.score.smooth), size = 1, shape = 16) +
  scale_color_gradientn(colours = Spectral(100), name = "Mediolateral score") +
  geom_line(data = pc.line, color = "red", size = 0.77) #+
p

startRes.shortlist <- p.learncurve$startRes.shortlist
object <- AddMetaData(object, metadata = meta.data)

var.genes <- GetVaryingGenes(object = object, meta.data = object@meta.data, genes.to.test = VariableFeatures(object))
genes.of.interest <- var.genes %>% pull(gene)

slot <- "data"
object <- GetResidual(object = object, features = genes.of.interest)
data.matrix <- GetAssayData(object = object, assay = "SCT", slot = slot)[genes.of.interest, ]
data.matrix.smooth <- SmoothSpline.Matrix(time.data = object$maturation.score, gene.matrix = data.matrix)

cells.cluster <- pam(as.dist((1 - cor(data.matrix.smooth, method = "spearman"))), k = 4)
cluster.info <- data.frame(cluster = paste0("cluster", cells.cluster$clustering), maturation.score = object$maturation.score)
boundaries <- aggregate(maturation.score ~ cluster, cluster.info, mean) %>%
  pull(maturation.score) %>%
  sort()

regions <- rev(c("MP", "DP", "LP", "VP", "AMY"))
z <- Hmisc::cut2(x = cluster.info$maturation.score, cuts = sort(boundaries))
levels(z) <- regions
object$DV.Partition <- z
Idents(object) <- "DV.Partition"


features.to.plot <- c("WNT8B", "SLIT2", "EMX1", "FOXP1", "SFRP1", "NEUROD4")
feature.plots <- list()

for (feature in features.to.plot) {
  p <- FeatureScatter(
    object = object,
    feature1 = "maturation.score",
    paste0("sct_", feature)
  ) + geom_smooth(method = "loess") + scale_color_manual(values = L1colors) + xlab("Mediolateral Score")
  feature.plots[[feature]] <- p
}
wrap_plots(feature.plots) + plot_layout(guides = "collect") & theme(legend.position = "right")


saveRDS(object = object, "objects/SLCHigh_stage_regressed_MedioLateralScored.rds")
