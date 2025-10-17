library(RColorBrewer)
library(ArchR)
library(dplyr)
library(ggplot2)
library(Seurat)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(ComplexHeatmap)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRThreads(threads = 10) 
addArchRGenome("mm10")

fiber <- loadArchRProject(path = "./",showLogo = F)

type2b <- fiber[fiber$anno_1103 %in% c("Type IIb Myonuclei","Sik1+ Myofiber")]
type2b

type2b <- addIterativeLSI(
    ArchRProj = type2b,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI1", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force =T
)
type2b <- addHarmony(
    ArchRProj = type2b,
    reducedDims = "IterativeLSI1",
    name = "Harmony",
    groupBy = "Sample",
    force =T
)
type2b <- addUMAP(
    ArchRProj = type2b, 
    reducedDims = "Harmony", 
    name = "harmony_UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
)
type2b <- addClusters(
    input = type2b,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "harmony_Clusters1",
    resolution = 0.8,
    maxClusters = NULL,
    force=T
)

pathToMacs2 <- "/home/anjuan/anaconda3/envs/stereopy/bin/macs2"
type2b <- addGroupCoverages(ArchRProj = type2b, groupBy = "harmony_Clusters1",force = T)
type2b <- addReproduciblePeakSet(
    ArchRProj = type2b,
    groupBy ="harmony_Clusters1",
    pathToMacs2 = pathToMacs2,
    force=T
)
type2b <- addPeakMatrix(type2b, force=T)
type2b <- addMotifAnnotations(ArchRProj = type2b, motifSet = "cisbp", name = "Motif")

trajectory <- c("C4", "C3", "C2","C5","C1")
type2b <- addTrajectory(
    ArchRProj = type2b, 
    name = "AgingdU", 
    groupBy = "harmony_Clusters1",
    trajectory = trajectory, 
    embedding = "harmony_UMAP", 
    force = TRUE
)

p <- plotTrajectory(type2b, trajectory = "AgingdU", colorBy = "cellColData", name = "AgingdU",embedding = "harmony_UMAP")
p[[1]]

trajMM  <- getTrajectory(ArchRProj = type2b, name = "AgingdU", useMatrix = "MotifMatrix", log2Norm = FALSE, trajectoryLabel = "harmony_Clusters1")
trajMM
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)))

trajGSM <- getTrajectory(ArchRProj = type2b, name = "AgingdU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE, trajectoryLabel = "harmony_Clusters1")
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajGSM)$label)))

trajPM  <- getTrajectory(ArchRProj = type2b, name = "AgingdU", useMatrix = "PeakMatrix", log2Norm = TRUE, trajectoryLabel = "harmony_Clusters1")
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "blueYellow"), colorColumns = TRUE, 
                            columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)))
                                                       
