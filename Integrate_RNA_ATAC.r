library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

rds <- readRDS("/mnt/6/anjuan/01.mouse_muscle/02.RNA_publish/Integrate.rds")
DefaultAssay(rds) <- 'RNA'
rds <- FindVariableFeatures(object = rds, nfeatures = 3000, verbose = FALSE)
#rds <- ScaleData(object = rds, verbose = FALSE)
genesUse <- VariableFeatures(object = rds)

proj <- loadArchRProject(path = "./",showLogo = F)
GeneScoreMatrix <- getMatrixFromProject(ArchRProj = proj,useMatrix = "GeneScoreMatrix")
gene_score <- assays(GeneScoreMatrix)$GeneScoreMatrix
rownames(gene_score) <- rowData(GeneScoreMatrix)$name
gene_score <- gene_score[rownames(gene_score) %in% genesUse,]
mat = log(gene_score+1)
atac <- CreateSeuratObject(
                counts = mat,
                assay = 'RNA',
                project = 'ATAC',
                min.cells = 1,
                meta.data = as.data.frame(proj@cellColData)
        )

atac$tech <- "ATAC"
rds$tech <- "RNA"
obj = merge(atac,rds)
ifnb.list <- SplitObject(obj, split.by = "tech")
features = rownames(atac)

rm(obj)
rm(atac)
rm(rds)
gc()

mMuscle.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
mMuscle.combined <- IntegrateData(anchorset = mMuscle.anchors)
DefaultAssay(mMuscle.combined) <- "integrated"
mMuscle.combined <- ScaleData(mMuscle.combined, verbose = FALSE)
mMuscle.combined <- RunPCA(mMuscle.combined, npcs = 30, verbose = FALSE)
mMuscle.combined <- RunUMAP(mMuscle.combined, reduction = "pca", dims = 1:30)
mMuscle.combined <- FindNeighbors(mMuscle.combined, reduction = "pca", dims = 1:30)
mMuscle.combined <- FindClusters(mMuscle.combined, resolution = 0.5)

saveRDS(mMuscle.combined,"Integrate_RNA_ATAC.rds")

