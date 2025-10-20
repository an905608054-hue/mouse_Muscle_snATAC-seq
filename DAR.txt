#####R
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ArchR)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(reshape2)
library(Seurat)
library(ComplexHeatmap)
addArchRThreads(threads = 10) 
addArchRGenome("mm10")

mone <- loadArchRProject(path = "./",showLogo = F)

pathToMacs2 <- "/home/wu_xinyu/miniconda3/bin/macs2"
mone <- addGroupCoverages(ArchRProj = mone, groupBy = "anno_group",force = T)
mone <- addReproduciblePeakSet(
    ArchRProj = mone, 
    groupBy = "anno_group", 
    pathToMacs2 = pathToMacs2,
    force=T
)
mone <- addPeakMatrix(mone,force = T)
mone <- addMotifAnnotations(mone, motifSet = "cisbp", name = "Motif",force = T)
pal = c("qMuSC_Adult"="#A15DBA","qMuSC_Old"="#2fa1dd","dMuSC_Adult"="#A15DBA","dMuSC_Old"="#2fa1dd",
        "ArtEC_Adult"="#A15DBA","ArtEC_Old"="#2fa1dd","CapEC_Adult"="#A15DBA","CapEC_Old"="#2fa1dd","VenEC_Adult"="#A15DBA","VenEC_Old"="#2fa1dd",
        "Gpc3+ FAP_Adult"="#A15DBA","Gpc3+ FAP_Old"="#2fa1dd","Mme+ FAP_Adult"="#A15DBA","Mme+ FAP_Old"="#2fa1dd","Cd55+ FAP_Adult"="#A15DBA","Cd55+ FAP_Old"="#2fa1dd",
        "B cell_Adult"="#A15DBA","B cell_Old"="#2fa1dd","Macrophage_Adult"="#A15DBA","Macrophage_Old"="#2fa1dd","T cell_Adult"="#A15DBA","T cell_Old"="#2fa1dd")

marker_CapEC <- getMarkerFeatures(
  ArchRProj = mone, 
  useMatrix = "PeakMatrix",
  groupBy = "anno_group",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CapEC_Old",
  bgdGroups = "CapEC_Adult"
)
pma_CapEC <- markerPlot(seMarker = marker_CapEC, name = "CapEC_Old", cutOff = "Pval <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "MA")
pma_CapEC
peak_CapEC <- cbind(marker_CapEC@elementMetadata,pma_CapEC$data)
peak_CapEC_up <- peak_CapEC[peak_CapEC$color == "Up-Regulated",]
peak_CapEC_up <- data.frame(peak_CapEC_up)
peak_CapEC_up <- peak_CapEC_up[,c(1,3,4)]
write.table(peak_CapEC_up,"./98.Figure4/00.bed/Peak_CapEC_up.bed",sep = "\t",row.names = F,col.names = F)
peak_CapEC_do <- peak_CapEC[peak_CapEC$color == "Down-Regulated",]
peak_CapEC_do <- data.frame(peak_CapEC_do)
peak_CapEC_do <- peak_CapEC_do[,c(1,3,4)]
write.table(peak_CapEC_do,"./98.Figure4/00.bed/Peak_CapEC_down.bed",sep = "\t",row.names = F,col.names = F)
write.csv(mone@peakSet,"98.Figure4/00.bed/Trend_All_celltype_peakset.csv")

#####Python
anno = pd.read_csv("98.Figure4/00.bed/Trend_All_celltype_peakset.csv")
anno = anno[['seqnames','start','end','nearestGene']]
anno['peak'] = anno.apply(lambda row:f"{row['seqnames']}-{row['start']}-{row['end']}",axis=1)
dMuSC_down = pd.read_csv("98.Figure4/00.bed/CapEC_p2g_down.bed",sep="\t",names=["seqnames", "start", "end", "strand"])
dMuSC_down['peak'] = dMuSC_down.apply(lambda row: f"{row['seqnames']}-{row['start']}-{row['end']}", axis=1)
merged_dMuSC_down = pd.merge(anno, dMuSC_down, on='peak', how='inner')
merged_dMuSC_down = merged_dMuSC_down.drop_duplicates(subset='nearestGene', keep='first')
merged_dMuSC_down.to_csv("98.Figure4/00.bed/CapEC_p2g_down_annotation.csv")
dMuSC_up = pd.read_csv("98.Figure4/00.bed/CapEC_p2g_up.bed",sep="\t",names=["seqnames", "start", "end", "strand"])
dMuSC_up['peak'] = dMuSC_up.apply(lambda row: f"{row['seqnames']}-{row['start']}-{row['end']}", axis=1)
merged_dMuSC_up = pd.merge(anno, dMuSC_up, on='peak', how='inner')
merged_dMuSC_up = merged_dMuSC_up.drop_duplicates(subset='nearestGene', keep='first')
merged_dMuSC_up.to_csv("98.Figure4/00.bed/CapEC_p2g_up_annotation.csv")


