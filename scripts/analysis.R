# NOTE: Update file paths to your local environment before running

## =======================
## Libraries
## =======================

###### =========================
###### Single-cell / Seurat
###### =========================

library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(scDblFinder)

###### =========================
###### Manipulación de datos / IO
###### =========================

library(dplyr)
library(tidyr)
library(tibble)
library(readr)

###### =========================
###### Estadística / modelos
###### =========================

library(Matrix)
library(apeglm)
library(DESeq2)

###### =========================
###### Enriquecimiento funcional
###### =========================

library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(GO.db)

###### =========================
###### Visualización general
###### =========================

library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(scales)

###### =========================
###### Visualización específica
###### =========================

library(pheatmap)
library(ggVennDiagram)
library(circlize)
library(RColorBrewer)


## =======================
## Data loading and preprocessing
## =======================

# Cargamos los RDS de Seven Bridges
sample1 <- read_rds("/home/anabel/Proyectos/Samanta_gyn/Samanta/data/seurats/sample1-tags_Seurat.rds")
sample2 <- read_rds("/home/anabel/Proyectos/Samanta_gyn/Samanta/data/seurats/sample2-tags_Seurat.rds")

# Añadimos un prefijo a los nombres de células para evitar solapamientos
sample1 <- RenameCells(sample1, add.cell.id = "sample1")
sample2 <- RenameCells(sample2, add.cell.id = "sample2")

# Normalización y selección de genes variables a nivel de muestra
sample1 <- NormalizeData(sample1)
sample1 <- FindVariableFeatures(sample1)
sample2 <- NormalizeData(sample2)
sample2 <- FindVariableFeatures(sample2)

# Merge simple (no integración todavía)
integrated_data <- merge(
  sample1,
  y = sample2,
  add.cell.ids = c("sample1", "sample2"),
  project = "AbSeq"
)

## =======================
## Metadata and condition assignment
## =======================

# Inicializamos como "otro"
integrated_data@meta.data$condition <- "otro"

# WT
integrated_data@meta.data$condition[
  integrated_data@meta.data$Sample_Name %in%
    c("sample1-1", "sample1-2", "sample1-3", "sample1-4")
] <- "WT"

# Casp1-/-
integrated_data@meta.data$condition[
  integrated_data@meta.data$Sample_Name %in%
    c("sample2-3", "sample2-4", "sample2-5", "sample2-6")
] <- "Casp1-/-"

# Nos quedamos sólo con WT y Casp1-/-
integrated_data <- subset(integrated_data, subset = condition != "otro")

# Opcional: comprobar distribución
# table(integrated_data@meta.data$Sample_Name, integrated_data@meta.data$condition)

## =======================
## Quality control
## =======================

integrated_data$percent.mt   <- PercentageFeatureSet(integrated_data, pattern = "^mt-")
integrated_data$percent.ribo <- PercentageFeatureSet(integrated_data, pattern = "^Rp[sl]")

VlnPlot(
  integrated_data,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol = 4
)

plot1 <- FeatureScatter(integrated_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(integrated_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(integrated_data, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot1 + plot2 + plot3

# Filtro de calidad
integrated_data <- subset(
  integrated_data,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 10000 &
    percent.mt < 25 &
    percent.ribo < 25
)

## =======================
## CD45 filtering
## =======================

integrated_data <- PercentageFeatureSet(
  integrated_data,
  features = "Ptprc",
  col.name = "percent_Ptprc"
)

noPtprc <- integrated_data[, integrated_data$percent_Ptprc == 0]

integrated_data_filter <- subset(
  integrated_data,
  subset = percent_Ptprc > 0
)

## =======================
## PCA and initial clustering (sin renormalizar)
## =======================

# integrated_data_filter <- NormalizeData(integrated_data_filter)

integrated_data_filter <- FindVariableFeatures(integrated_data_filter)
integrated_data_filter <- ScaleData(integrated_data_filter, verbose = FALSE)
integrated_data_filter <- RunPCA(integrated_data_filter, npcs = 30, verbose = FALSE)

# Necesitamos clústeres preliminares antes de modelHomotypic
integrated_data_filter <- FindNeighbors(integrated_data_filter, dims = 1:30)
integrated_data_filter <- FindClusters(integrated_data_filter, resolution = 0.5)


## =======================
## Doublet removal per sample (scDblFinder)
## =======================

DefaultAssay(integrated_data_filter) <- "RNA"
integrated_data_filter[["RNA"]] <- as(integrated_data_filter[["RNA"]], Class = "Assay")

samples_list <- SplitObject(integrated_data_filter, split.by = "Sample_Name")

# sample1-1
sce1 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample1-1"]]))
samples_list[["sample1-1"]]$doublet_class <- sce1$scDblFinder.class
samples_list[["sample1-1"]] <- subset(samples_list[["sample1-1"]], subset = doublet_class == "singlet")

# sample1-2
sce2 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample1-2"]]))
samples_list[["sample1-2"]]$doublet_class <- sce2$scDblFinder.class
samples_list[["sample1-2"]] <- subset(samples_list[["sample1-2"]], subset = doublet_class == "singlet")

# sample1-3
sce3 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample1-3"]]))
samples_list[["sample1-3"]]$doublet_class <- sce3$scDblFinder.class
samples_list[["sample1-3"]] <- subset(samples_list[["sample1-3"]], subset = doublet_class == "singlet")

# sample1-4
sce4 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample1-4"]]))
samples_list[["sample1-4"]]$doublet_class <- sce4$scDblFinder.class
samples_list[["sample1-4"]] <- subset(samples_list[["sample1-4"]], subset = doublet_class == "singlet")

# sample2-3
sce5 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample2-3"]]))
samples_list[["sample2-3"]]$doublet_class <- sce5$scDblFinder.class
samples_list[["sample2-3"]] <- subset(samples_list[["sample2-3"]], subset = doublet_class == "singlet")

# sample2-4
sce6 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample2-4"]]))
samples_list[["sample2-4"]]$doublet_class <- sce6$scDblFinder.class
samples_list[["sample2-4"]] <- subset(samples_list[["sample2-4"]], subset = doublet_class == "singlet")

# sample2-5
sce7 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample2-5"]]))
samples_list[["sample2-5"]]$doublet_class <- sce7$scDblFinder.class
samples_list[["sample2-5"]] <- subset(samples_list[["sample2-5"]], subset = doublet_class == "singlet")

# sample2-6
sce8 <- scDblFinder(as.SingleCellExperiment(samples_list[["sample2-6"]]))
samples_list[["sample2-6"]]$doublet_class <- sce8$scDblFinder.class
samples_list[["sample2-6"]] <- subset(samples_list[["sample2-6"]], subset = doublet_class == "singlet")

# Merge otra vez
integrated_data_filter <- merge(
  samples_list[[1]],
  y = samples_list[2:length(samples_list)]
)

## =======================
## Data integration
## =======================

# Nos aseguramos de que condition sea factor
integrated_data_filter$condition <- factor(integrated_data_filter$condition)

# Dividimos por condición (WT y Casp1-/-)
samples_list <- SplitObject(integrated_data_filter, split.by = "Sample_Name")

samples_list[[1]] <- NormalizeData(samples_list[[1]])
samples_list[[1]] <- FindVariableFeatures(samples_list[[1]])

samples_list[[2]] <- NormalizeData(samples_list[[2]])
samples_list[[2]] <- FindVariableFeatures(samples_list[[2]])

samples_list[[3]] <- NormalizeData(samples_list[[3]])
samples_list[[3]] <- FindVariableFeatures(samples_list[[3]])

samples_list[[4]] <- NormalizeData(samples_list[[4]])
samples_list[[4]] <- FindVariableFeatures(samples_list[[4]])

samples_list[[5]] <- NormalizeData(samples_list[[5]])
samples_list[[5]] <- FindVariableFeatures(samples_list[[5]])

samples_list[[6]] <- NormalizeData(samples_list[[6]])
samples_list[[6]] <- FindVariableFeatures(samples_list[[6]])

samples_list[[7]] <- NormalizeData(samples_list[[7]])
samples_list[[7]] <- FindVariableFeatures(samples_list[[7]])

samples_list[[8]] <- NormalizeData(samples_list[[8]])
samples_list[[8]] <- FindVariableFeatures(samples_list[[8]])

# Encontrar anclajes
anchors <- FindIntegrationAnchors(object.list = samples_list, dims = 1:30)

# Integrar
integrated_data_filter <- IntegrateData(anchorset = anchors, dims = 1:30)

# Escalado y PCA en datos integrados
DefaultAssay(integrated_data_filter) <- "integrated"

integrated_data_filter <- ScaleData(integrated_data_filter, verbose = FALSE)
integrated_data_filter <- RunPCA(integrated_data_filter, npcs = 30, verbose = FALSE, reduction.name = "pca_integrated")

# Vecinos y clústeres sobre PCA integrado
integrated_data_filter <- FindNeighbors(integrated_data_filter, dims = 1:30, reduction = "pca_integrated")
integrated_data_filter <- FindClusters(integrated_data_filter, resolution = 0.5)

## =======================
## UMAP and clustering
## =======================

integrated_data_filter <- RunUMAP(
  integrated_data_filter,
  reduction = "pca_integrated",
  dims = 1:30,
  reduction.name = "umap_rna5"
)

## =======================
## Cluster filtering
## =======================

clusters.validos <- names(which(table(Idents(integrated_data_filter)) >= 100))

integrated_data_filter <- subset(integrated_data_filter, idents = clusters.validos)

table(Idents(integrated_data_filter))

# MUY IMPORTANTE: recalcular tras eliminar células
integrated_data_filter <- FindNeighbors(
  integrated_data_filter,
  dims = 1:30,
  reduction = "pca_integrated"
)

integrated_data_filter <- FindClusters(
  integrated_data_filter,
  resolution = 0.5
)

integrated_data_filter <- RunUMAP(
  integrated_data_filter,
  reduction      = "pca_integrated",
  dims           = 1:30,
  reduction.name = "umap_rna5"
)

table(Idents(integrated_data_filter))

# Si quieres guardar el objeto en este punto:
saveRDS(integrated_data_filter, "integrated_data_despuesUMAP_sample12.rds")


## =======================
## Visualization
## =======================

DimPlot(integrated_data_filter, reduction = "umap_rna5", label = TRUE)
ggsave("./results/DimPlot_filtrado.png", width = 1920, height = 980, units = "px", scale = 2)

DimPlot(integrated_data_filter, reduction = "umap_rna5", group.by = "condition")
ggsave("./results/DimPlot_filtrado_condition.png", width = 1920, height = 980, units = "px", scale = 2)


tabla <- table(Idents(integrated_data_filter), integrated_data_filter$condition)
prop.table(tabla, margin = 1)


DimPlot(integrated_data_filter, label = TRUE, split.by  = "condition") + labs(x = "UMAP_1", y = "UMAP_2")
ggsave("./results/DimPlot_condition_label_split_sample12.svg", width = 1920, height = 980, units = "px", scale = 2)

## =======================
## Cell type annotation
## =======================

DefaultAssay(integrated_data_filter) <- "RNA"

### T cells 
FeaturePlot(integrated_data_filter,
            features = c("Cd3d","Cd3e","Cd3g","Cd8a","Cd4","Cd28",
                         "Gzma","Gzmb","Il2rb","Il7r","Tcf7","Lef1",
                         "Tcra","Tcrb"),
            ncol = 4, label = TRUE)
ggsave("./results/FeaturePlot_Tcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Cd3d","Cd3e","Cd3g","Cd8a","Cd4","Cd28",
                       "Gzma","Gzmb","Il2rb","Il7r","Tcf7","Lef1",
                       "Tcra","Tcrb"),
          ncol = 4)
ggsave("./results/RidgePlot_Tcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Cd3d","Cd3e","Cd3g","Cd8a","Cd4","Cd28",
                     "Gzma","Gzmb","Il2rb","Il7r","Tcf7","Lef1",
                     "Tcra","Tcrb"),
        ncol = 4)
ggsave("./results/VlnPlot_Tcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### Basófilos 
FeaturePlot(integrated_data_filter,
            features = c("Mcpt8","Gata2","Cpa3","Ms4a2","Fcer1a"),
            ncol = 3, label = TRUE)
ggsave("./results/FeaturePlot_Basofilo5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Mcpt8","Gata2","Cpa3","Ms4a2","Fcer1a"),
          ncol = 3)
ggsave("./results/RidgePlot_Basofilo5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Mcpt8","Gata2","Cpa3","Ms4a2","Fcer1a"),
        ncol = 3)
ggsave("./results/VlnPlot_Basofilo5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### B cells 
FeaturePlot(integrated_data_filter,
            features = c("Cd79a","Cd79b","Cd19","Ly6d","Ebf1",
                         "Ms4a1","Ighm","Ighd"),
            ncol = 4, label = TRUE)
ggsave("./results/FeaturePlot_Bcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Cd79a","Cd79b","Cd19","Ly6d","Ebf1",
                       "Ms4a1","Ighm","Ighd"),
          ncol = 4)
ggsave("./results/RidgePlot_Bcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Cd79a","Cd79b","Cd19","Ly6d","Ebf1",
                     "Ms4a1","Ighm","Ighd"),
        ncol = 4)
ggsave("./results/VlnPlot_Bcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### Macrófagos 
FeaturePlot(integrated_data_filter,
            features = c("Itgam","Adgre1","Cd68","Csf1r","Lyz2",
                         "Clec4f","Apoe","C1qa","C1qb","C1qc","Cd5l","Ctss"),
            ncol = 4, label = TRUE)
ggsave("./results/FeaturePlot_Macrophage5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Itgam","Adgre1","Cd68","Csf1r","Lyz2",
                       "Clec4f","Apoe","C1qa","C1qb","C1qc","Cd5l","Ctss"),
          ncol = 4)
ggsave("./results/RidgePlot_Macrophage5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Itgam","Adgre1","Cd68","Csf1r","Lyz2",
                     "Clec4f","Apoe","C1qa","C1qb","C1qc","Cd5l","Ctss"),
        ncol = 4)
ggsave("./results/VlnPlot_Macrophage5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### ILC1 
FeaturePlot(integrated_data_filter,
            features = c("Itgae","Tnfsf10","Cd27"),
            ncol = 3, label = TRUE)
ggsave("./results/FeaturePlot_ILC15_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Itgae","Tnfsf10","Cd27"),
          ncol = 3)
ggsave("./results/RidgePlot_ILC15_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Itgae","Tnfsf10","Cd27"),
        ncol = 3)
ggsave("./results/VlnPlot_ILC15_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### NK cells 
FeaturePlot(integrated_data_filter,
            features = c("Ncr1","Klrb1c","Klra8","Klrb1b","Gzma","Nkg7"),
            ncol = 4, label = TRUE)
ggsave("./results/FeaturePlot_NKcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Ncr1","Klrb1c","Klra8","Klrb1b","Gzma","Nkg7"),
          ncol = 4)
ggsave("./results/RidgePlot_NKcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Ncr1","Klrb1c","Klra8","Klrb1b","Gzma","Nkg7"),
        ncol = 4)
ggsave("./results/VlnPlot_NKcell5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### DC 
FeaturePlot(integrated_data_filter,
            features = c("Itgax","Siglech","Cd80","Xcr1","Irf8",
                         "Flt3l","Batf3","Itgae","Sirpa"),
            ncol = 4, label = TRUE)
ggsave("./results/FeaturePlot_DC5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Itgax","Siglech","Cd80","Xcr1","Irf8",
                       "Flt3l","Batf3","Itgae","Sirpa"),
          ncol = 4)
ggsave("./results/RidgePlot_DC5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Itgax","Siglech","Cd80","Xcr1","Irf8",
                     "Flt3l","Batf3","Itgae","Sirpa"),
        ncol = 4)
ggsave("./results/VlnPlot_DC5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

### Neutrófilos 
FeaturePlot(integrated_data_filter,
            features = c("Ly6g","Csf3r","Itgam","S100a9","S100a8",
                         "Lcn2","Ccrl2","Cxcr2","Retnlg"),
            ncol = 4, label = TRUE)
ggsave("./results/FeaturePlot_Neutrofilos5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

RidgePlot(integrated_data_filter,
          features = c("Ly6g","Csf3r","Itgam","S100a9","S100a8",
                       "Lcn2","Ccrl2","Cxcr2","Retnlg"),
          ncol = 4)
ggsave("./results/RidgePlot_Neutrofilos5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

VlnPlot(integrated_data_filter,
        features = c("Ly6g","Csf3r","Itgam","S100a9","S100a8",
                     "Lcn2","Ccrl2","Cxcr2","Retnlg"),
        ncol = 4)
ggsave("./results/VlnPlot_Neutrofilos5_integrated_sample12.png", width = 1920, height = 980, units = "px", scale = 2)

## =======================
## Subsetting cell populations
## =======================

integrated_data_filter<-read_rds("integrated_data_despuesUMAP_sample12.rds")

integrated_data_filter@meta.data$condition <- as.factor(
  integrated_data_filter@meta.data$condition
)
integrated_data_filter@meta.data$condition <- relevel(
  integrated_data_filter@meta.data$condition, ref = "WT"
)

DefaultAssay(integrated_data_filter)<-"RNA"

integrated_data_filter1<- subset(integrated_data_filter, idents=1)

integrated_data_filter <- RenameIdents(integrated_data_filter, 
                                       "0" = "B",
                                       "1" = "CD8+", 
                                       "2" = "Mo/Mac",
                                       "3" = "CD8+",
                                       "4" = "Mo/Mac",
                                       "5" = "NK",
                                       "6" = "CD4+",
                                       "7" = "B",
                                       "8" = 'Neu',
                                       "9" = "DC",
                                       "10" = "Mo/Mac",
                                       "11" = "Mo/Mac",
                                       "12" = "CD4+",
                                       "13" = "Other",
                                       "14" = "NK",
                                       "15" = "DC",
                                       "16" = "Mo/Mac", 
                                       "17" = "CD4+")

### Cluster selection

integrated_data_filter_mo <- subset(integrated_data_filter, idents='Mo/Mac')

integrated_data_filter_cd8 <- subset(integrated_data_filter, idents='CD8+')

integrated_data_filternoNA <- subset(integrated_data_filter, 
                                     cells = WhichCells(integrated_data_filter, 
                                                        idents = setdiff(levels(Idents(integrated_data_filter)), "Other")))

png("./results/DimPlot_condition_label_split_sample12_label.png",width = 2000,
    height = 1600,
    res = 300)

DimPlot(integrated_data_filternoNA, reduction = "umap_rna5",label=FALSE,split.by = "condition") +
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  )

dev.off()

## =======================
## Plots
## =======================

##Violin plot Tnfa en Cluster Mo/Mac
df <- FetchData(integrated_data_filter_mo, vars = c("Tnf", "condition"))

png("./results/Violin_plot_tnf_mo_mac.png",width = 2000,
    height = 1600,
    res = 300)
ggplot(df, aes(x = condition, y = Tnf, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 3,
    color = "black"
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_fill_manual(values = c("WT" = "indianred3", "Casp1-/-" = "skyblue")) +
  theme_classic() +
  labs(
    title = "Tnf expression",
    x = "",
    y = "Expression level"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

dev.off()

## Heatmap Mo/Mac de: "Cxcl2","Ccl5","Ccl4","Ccl3","Ccl2","Cxcr4","Cxcl1","Cd83","Cd14"

obj <- integrated_data_filter_mo
Idents(obj) <- "Sample_Name"

# Average por ratón
avg.mat <- AverageExpression(
  obj,
  assays = "RNA",
  slot = "data"
)$RNA

genes <- c("Cxcl2","Ccl5","Ccl4","Ccl3","Ccl2","Cxcr4","Cxcl1","Cd83","Cd14")
mat <- avg.mat[genes, ]

# metadata
meta <- obj@meta.data
sample_condition <- unique(meta[, c("Sample_Name", "condition")])
sample_condition <- sample_condition[match(colnames(mat), sample_condition$Sample_Name), ]

# separar WT y KO
wt_samples <- sample_condition$Sample_Name[sample_condition$condition == "WT"]
ko_samples <- sample_condition$Sample_Name[sample_condition$condition == "Casp1-/-"]

mat_wt <- mat[, wt_samples]
mat_ko <- mat[, ko_samples]

# clustering dentro de cada grupo
mat_wt <- mat_wt[, hclust(dist(t(mat_wt)))$order]
mat_ko <- mat_ko[, hclust(dist(t(mat_ko)))$order]

# recombinar
mat_final <- cbind(mat_wt, mat_ko)

# annotation
annotation_col <- data.frame(
  condition = c(rep("WT", ncol(mat_wt)), rep("Casp1-/-", ncol(mat_ko)))
)
rownames(annotation_col) <- colnames(mat_final)

# colores
annotation_colors <- list(
  condition = c("WT" = "indianred3", "Casp1-/-" = "skyblue")
)

# heatmap
png("heatmap_momac.png",
    width = 2000,
    height = 1600,
    res = 300)
pheatmap(
  mat_final,
  scale = "row",
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = FALSE, 
  main = "Mo/Mac gene expression per mouse",
  show_colnames = FALSE,
  fontsize_row = 12,   
  fontsize_col = 14,
  cellwidth = 20,      
  cellheight = 20,     
  border_color = NA
)

dev.off()

##Heatmap CD8+ de: "Myc", "Cxcr4", "Ide","Ifng", "Gzmk","Tcf7", "Il7r","Pdcd1", "Havcr2", "Tox", "Ctla4", "Tigit","Batf", "Mt1", "Ikzf2", "Tnfrsf1b"

## =======================
## Preparar objeto
## =======================

obj_cd8 <- integrated_data_filter_cd8
Idents(obj_cd8) <- "Sample_Name"

## =======================
## Average por ratón
## =======================

avg.mat <- AverageExpression(
  obj_cd8,
  assays = "RNA",
  slot = "data"
)$RNA

## =======================
## Orden EXACTO de genes 
## =======================

genes_ordered <- c(
  "Myc", "Cxcr4", "Ide",
  "Ifng", "Gzmk",
  "Tcf7", "Il7r",
  "Pdcd1", "Havcr2", "Tox", "Ctla4", "Tigit",
  "Batf", "Mt1", "Ikzf2", "Tnfrsf1b"
)

## =======================
## Filtrar genes presentes
## =======================

genes_use <- genes_ordered[genes_ordered %in% rownames(avg.mat)]
mat <- avg.mat[genes_use, , drop = FALSE]

## =======================
## Metadata
## =======================

meta <- obj_cd8@meta.data
sample_condition <- unique(meta[, c("Sample_Name", "condition")])

rownames(sample_condition) <- sample_condition$Sample_Name
sample_condition <- sample_condition[colnames(mat), ]

wt_samples <- sample_condition$Sample_Name[sample_condition$condition == "WT"]
ko_samples <- sample_condition$Sample_Name[sample_condition$condition == "Casp1-/-"]

mat_wt <- mat[, wt_samples, drop = FALSE]
mat_ko <- mat[, ko_samples, drop = FALSE]

## =======================
## Orden columnas (opcional)
## =======================

if(ncol(mat_wt) > 1){
  mat_wt <- mat_wt[, hclust(dist(t(mat_wt)))$order]
}

if(ncol(mat_ko) > 1){
  mat_ko <- mat_ko[, hclust(dist(t(mat_ko)))$order]
}

mat_samples <- cbind(mat_wt, mat_ko)

## =======================
## Anotación columnas
## =======================

annotation_col <- data.frame(
  condition = c(rep("WT", ncol(mat_wt)), rep("Casp1-/-", ncol(mat_ko)))
)
rownames(annotation_col) <- colnames(mat_samples)

annotation_colors <- list(
  condition = c("WT" = "indianred3", "Casp1-/-" = "skyblue")
)

## =======================
##  HEATMAP FINAL
## =======================

png("heatmap_CD8_final.png", width = 2000, height = 1600, res = 300)

pheatmap(
  mat_samples,
  scale = "row",  # cambia a "none" si quieres valores reales
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize_row = 12,   
  fontsize_col = 14,
  cellwidth = 20,     
  cellheight = 20,    
  border_color = NA,
  main = "CD8+ gene expression per mouse"
)

dev.off()

## =======================
## Pseudobulk DESeq2 analysis
## =======================

# Agregación por suma (pseudobulk) y extracción

cts_list <- AggregateExpression(
  integrated_data_filter,
  group.by = c("ident","condition","Sample_Name"),
  assays = "RNA",
  slot = "counts",
  functions = "sum",
  return.seurat = FALSE
)
cts <- cts_list$RNA

# Transponer y convertir
cts.t <- t(cts) %>% as.data.frame()

# Split por ident (prefijo antes de "_")
splitRows <- sub("_.*", "", rownames(cts.t))
cts.split <- split.data.frame(cts.t, f = factor(splitRows))

# Arreglar rownames: de "CD8+_WT_sample1" -> "WT|sample1"
cts.split.modified.all <- lapply(cts.split, function(x){
  rn <- rownames(x)
  # quitar el ident y el primer "_"
  rn <- sub("^[^_]*_", "", rn)           # "WT_sample1"
  # convertir "_" por "|" para luego separar: "WT|sample1"
  rn <- sub("_", "|", rn, fixed = TRUE)
  rownames(x) <- rn
  t(x)                                   # genes x muestras
})

write_rds(cts.split.modified.all, file = "./results/cts.split.modified.all.rds")
cts.split.modified.all <- readRDS(file = "./results/cts.split.modified.all.rds")

### Loop por cluster 
for (i in c("Mo/Mac","CD8+")) {
  message(i)
  counts <- cts.split.modified.all[[i]]
  
  # Construir colData a partir de colnames tipo "WT|sample1"
  colData <- data.frame(condition = factor(colnames(counts)), row.names = colnames(counts))
  colData <- tidyr::separate(colData, col = condition, into = c("condition","sample"),
                             sep = "\\|", remove = FALSE)
  
  # Asegúrate de que las columnas de counts sean enteros
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                                        colData = colData,
                                        design = ~ condition)
  
  # Nivel de referencia
  dds$condition <- relevel(dds$condition, ref = "WT")
  
  # Filtrado mínimo
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  dds <- DESeq2::DESeq(dds)
  levels(dds$condition)
  
  # Contraste Casp1-/- vs WT
  res <- DESeq2::results(dds, contrast = c("condition","Casp1-/-","WT"))
  
  ## Casp1-/- vs WT
  res_CASP1_vs_WT <- DESeq2::results(dds, contrast = c("condition", "Casp1-/-", "WT"))
  
  if(i == "Mo/Mac"){
    write.csv(res_CASP1_vs_WT, file = paste0("./results/Mo_Mac_res_CASP1_vs_WT",".csv"),row.names = T,sep=";")
  } else {
    write.csv(res_CASP1_vs_WT, file = paste0("./results/",i,"_res_CASP1_vs_WT",".csv"),row.names = T,sep=";")
  }
  res_CASP1_vs_WT_df <- res_CASP1_vs_WT %>%
    as.data.frame() %>%
    dplyr::mutate(significant = !is.na(padj) &  padj < 0.05)
  
  res_CASP1_vs_WT_df <- res_CASP1_vs_WT_df %>%
    tibble::rownames_to_column(var="genes")
  
  if(i == "Mo/Mac"){
    write.csv(res_CASP1_vs_WT_df, file = paste0("./results/Mo_Mac_res_CASP1_vs_WT_df_0.05",".csv"),row.names = T,sep=";")
  } else {
    write.csv(res_CASP1_vs_WT_df, file = paste0("./results/",i,"_res_CASP1_vs_WT_df_0.05",".csv"),row.names = T,sep=";")
  }
 
  ## plotMA
  
 p <- ggmaplot(
    res_CASP1_vs_WT_df,
    main = paste(i, "plotMA_labelled", sep = "_"),
    padj = 0.05, size = 0.4,fc.cutoff = 0,
    palette = c("palevioletred4", "aquamarine4", "darkgray"),
    genenames = res_CASP1_vs_WT_df$genes,   # <- mismo largo que nrow(data)
    legend = "top", top = 0,               # <- sin auto-etiquetado
    font.label = c("bold", 11),
    label.rectangle = TRUE,
    font.legend = "bold",
    font.main = "bold",
    ggtheme = ggplot2::theme_minimal()
  ) + geom_label_repel(
    data = subset(res_CASP1_vs_WT_df, genes %in% c("Tnfaip3")),
    aes(x = log2(baseMean + 1), y = log2FoldChange, label = genes),
    size = 3,
    fontface = "bold",
    fill = "white",
    color = "black",
    # --- que la línea se vea sí o sí ---
    min.segment.length = 0,      # fuerza a dibujar el segmento
    segment.color = "black",
    segment.size = 0.7,
    segment.linetype = "solid",
    # --- separaciones para que no quede debajo de la caja ---
    box.padding   = unit(0.5, "lines"),
    point.padding = unit(0.4, "lines"),
    # mueve un poco la etiqueta para despegarla del punto
    nudge_x = 0.15,              # ajusta según tu escala
    nudge_y = 0.15,
    # línea recta y clara
    segment.curvature = 0,
    segment.ncp = 1,
    max.overlaps = Inf
  )
 
  if(i == "Mo/Mac"){
    ggsave(
      "./results/Mo_Mac_plotMA_Tnfaip3_CASP1_vs_WT.png",
      plot = p,
      width = 12,
      height = 6,
      units = "in",
      dpi = 300   # 🔥 clave
    )
  } else {
    ggsave(
      paste0("./results/", i, "_plotMA_Tnfaip3_CASP1_vs_WT.png"),
      plot = p,
      width = 12,
      height = 6,
      units = "in",
      dpi = 300
    )
  }
  
}

## =======================
## Venn diagram
## =======================

mo_data<-read_delim("./results/Mo_Mac_res_CASP1_vs_WT_df_0.05.csv", 
                    delim = ",", escape_double = FALSE, trim_ws = TRUE) %>% filter(significant==TRUE)
cd8_data<-read_delim("./results/CD8+_res_CASP1_vs_WT_df_0.05.csv", 
                     delim = ",", escape_double = FALSE, trim_ws = TRUE) %>% filter(significant==TRUE)

venn_data <- list(
  CD8 = cd8_data$genes,
  Mo  = mo_data$genes
)

p<-ggVennDiagram(venn_data, label_alpha = 0) +
  scale_fill_gradient(low = "lightskyblue", high = "lightskyblue3") +
  labs(title = "Significantly expressed genes (padj < 0.05)") +
  theme_void()


ggsave(
  "./results/Venn_diagram_Mo_vs_CD8.png",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300   # 🔥 clave para buena calidad
)

final_venn<-inner_join(mo_data,cd8_data,by='genes') %>% dplyr::select(genes)
write.table(final_venn, file = "common_genes_venn.csv",row.names = F,sep=";")

## =======================
## GSEA analysis
## =======================

## Crear los gene sets GO BP

go_bp <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys     = keys(org.Mm.eg.db, keytype = "ENTREZID"),
  columns  = c("GO", "ONTOLOGY"),
  keytype  = "ENTREZID"
) %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  dplyr::select(GO, ENTREZID) %>%
  dplyr::distinct()

go_names <- AnnotationDbi::select(
  GO.db,
  keys     = unique(go_bp$GO),
  columns  = "TERM",
  keytype  = "GOID"
) %>%
  dplyr::rename(
    GO = GOID,
    Description = TERM
  ) %>%
  dplyr::distinct()

### Global GSEA (Mo/Mac)

## Leer resultados DESeq2

res_mo <- read_delim(
  "./results/Mo_Mac_res_CASP1_vs_WT.csv",
  delim = ",",
  locale = locale(decimal_mark = "."),
  trim_ws = TRUE
)

## Preparar ranking 
res_mo <- res_mo %>%
  dplyr::filter(!is.na(stat))

geneList <- res_mo$stat
names(geneList) <- res_mo$...1
geneList <- sort(geneList, decreasing = TRUE)

## SYMBOL -> ENTREZ 
gene_df <- bitr(
  names(geneList),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

geneList_entrez <- geneList[gene_df$SYMBOL]
names(geneList_entrez) <- gene_df$ENTREZID

geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

## GSEA
gsea_ct <- GSEA(
  geneList      = geneList_entrez,
  TERM2GENE     = go_bp,
  TERM2NAME     = go_names,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

## Hacer legible
gsea_ct <- setReadable(
  gsea_ct,
  OrgDb   = org.Mm.eg.db,
  keyType = "ENTREZID"
)

write.csv2(
  as.data.frame(gsea_ct),
  paste0("results/Mo_mac_GSEA_GO_BP__0.65.csv"),
  row.names = FALSE
)

### Targeted GSEA (Mo/Mac)

go_interest <- c(
  "GO:0043065", # apoptosis
  "GO:0006955", # immune response
  "GO:0010628", # gene expression
  "GO:0000165", # MAPK
  "GO:0035556", # intracellular signaling
  "GO:0071456", # hypoxia
  "GO:0032731", # IL-1 beta
  "GO:0002474"  # antigen presentation MHC I
)

go_bp_subset <- go_bp %>%
  dplyr::filter(GO %in% go_interest)

go_names_subset <- go_names %>%
  dplyr::filter(GO %in% go_interest)

## Leer DESeq2
res_mo <- read_delim(
  "./results/Mo_Mac_res_CASP1_vs_WT.csv",
  delim = ",",
  locale = locale(decimal_mark = "."),
  trim_ws = TRUE
)

## Ranking
res_mo <- res_mo %>%
  dplyr::filter(!is.na(stat))

geneList <- res_mo$stat
names(geneList) <- res_mo$...1

geneList <- sort(geneList, decreasing = TRUE)


## SYMBOL → ENTREZ

gene_df <- bitr(
  names(geneList),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

geneList_entrez <- geneList[gene_df$SYMBOL]
names(geneList_entrez) <- gene_df$ENTREZID

geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

## GSEA dirigido
gsea_ct_subset <- GSEA(
  geneList      = geneList_entrez,
  TERM2GENE     = go_bp_subset,
  TERM2NAME     = go_names_subset,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,   
  pAdjustMethod = "BH"
)

## Legible
gsea_ct_subset <- setReadable(
  gsea_ct_subset,
  OrgDb   = org.Mm.eg.db,
  keyType = "ENTREZID"
)

write.csv2(
  as.data.frame(gsea_ct_subset),
  "results/Mo_mac_GSEA_targeted_GO_BP.csv",
  row.names = FALSE
)

## Visualización
p <- dotplot(gsea_ct_subset, showCategory = 20) +
  ggtitle("Mo/Mac – Targeted GSEA GO BP") +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text  = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title  = element_text(size = 14, face = "bold")
  )

ggsave(
  "./results/GSEA_targeted_BP_mo_mac.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

#### =========================
#### CD8+
#### =========================

### Global GSEA (CD8)

## Leer resultados DESeq2

res_cd8 <- CD8_res_CASP1_vs_WT <- read_delim("results/CD8+_res_CASP1_vs_WT.csv", locale = locale(decimal_mark = "."),
                                             delim = ",", escape_double = FALSE, trim_ws = TRUE)

## Preparar ranking
res_cd8 <- res_cd8 %>%
  dplyr::filter(!is.na(stat))

geneList <- res_cd8$stat
names(geneList) <- res_cd8$...1
geneList <- sort(geneList, decreasing = TRUE)

## SYMBOL -> ENTREZ
gene_df <- bitr(
  names(geneList),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

geneList_entrez <- geneList[gene_df$SYMBOL]
names(geneList_entrez) <- gene_df$ENTREZID

geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

## GSEA clásico
gsea_ct <- GSEA(
  geneList      = geneList_entrez,
  TERM2GENE     = go_bp,
  TERM2NAME     = go_names,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

## Hacer legible
gsea_ct <- setReadable(
  gsea_ct,
  OrgDb   = org.Mm.eg.db,
  keyType = "ENTREZID"
)

write.csv2(
  as.data.frame(gsea_ct),
  paste0("results/cd8_GSEA_GO_BP__0.65.csv"),
  row.names = FALSE
)

### Targeted GSEA (CD8)

go_interest <- c(
  "GO:0042127","GO:0043410","GO:0007166","GO:0006468",
  "GO:0000165","GO:0030154","GO:0010628","GO:0002250",
  "GO:0009749","GO:0050852","GO:0006955","GO:0006954",
  "GO:1990869","GO:0034341","GO:0071356","GO:0042149",
  "GO:0034612","GO:0001819"
)

go_bp_subset <- go_bp %>%
  dplyr::filter(GO %in% go_interest)

go_names_subset <- go_names %>%
  dplyr::filter(GO %in% go_interest)

## Leer DESeq2
res_cd8 <- read_delim(
  "./results/CD8+_res_CASP1_vs_WT.csv",
  delim = ",",
  locale = locale(decimal_mark = "."),
  trim_ws = TRUE
)

## Ranking
res_cd8 <- res_cd8 %>%
  dplyr::filter(!is.na(stat))

geneList <- res_cd8$stat
names(geneList) <- res_cd8$...1

geneList <- sort(geneList, decreasing = TRUE)

## SYMBOL -> ENTREZ
gene_df <- bitr(
  names(geneList),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

geneList_entrez <- geneList[gene_df$SYMBOL]
names(geneList_entrez) <- gene_df$ENTREZID

geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

## GSEA dirigido
gsea_cd8_subset <- GSEA(
  geneList      = geneList_entrez,
  TERM2GENE     = go_bp_subset,
  TERM2NAME     = go_names_subset,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,   
  pAdjustMethod = "BH"
)

## Legible
gsea_cd8_subset <- setReadable(
  gsea_cd8_subset,
  OrgDb   = org.Mm.eg.db,
  keyType = "ENTREZID"
)

write.csv2(
  as.data.frame(gsea_cd8_subset),
  "results/CD8_GSEA_targeted_GO_BP.csv",
  row.names = FALSE
)

## Visualización
GSEA_cd8 <- dotplot(gsea_cd8_subset, showCategory = 20) +
  ggtitle("CD8+ – Targeted GSEA GO BP") +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text  = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title  = element_text(size = 14, face = "bold")
  )

ggsave(
  "./results/GSEA_targeted_BP_cd8.png",
  plot = GSEA_cd8,
  width = 8,
  height = 6,
  dpi = 300
)

