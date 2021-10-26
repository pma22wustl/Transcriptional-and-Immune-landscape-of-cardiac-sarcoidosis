#package install
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
## Required
BiocManager::install(c("AUCell", "RcisTarget"),ask = F,update = F) 
BiocManager::install(c("GENIE3"),ask = F,update = F)  # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"),ask = F,update = F) 
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),ask = F,update = F)
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"),ask = F,update = F)
install.packages("doMC", repos="http://R-Forge.R-project.org")
library(doMC)
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")

BiocManager::install("BiocParallel")
library(BiocParallel)


#download database
#for human
mm_url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather"
mm_url2="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather"
#for mouse
hg_url="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather"
hg_url2="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather"

wget -c $mm_url
wget -c $mm_url2
wget -c $hg_url
wget -c $hg_url2

library(RcisTarget)
# Explore tutorials in the web browser:
browseVignettes(package="RcisTarget")  
vignette("RcisTarget") # open


#use Seurat object to do the analysis
rm(list = ls())
install.packages("Seurat")
library(Seurat) 
seuratobj <- readRDS("Renamed_Myeloid_UMAP.new.ident.rds")
exprMat  <-  as.matrix(seuratobj@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4] 

# prepare cellInfo (dataframe with "CellType",'nGene' ,'nUMI')
cellInfo <-  seuratobj@meta.data[,c(19,2,3)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("0"="forestgreen",    #"Mac_res",  
                           "1"="forestgreen",
                           "2"="darkorange",     #"Mac_VCAN",
                           "3"="magenta4",     #"Mac_HLA-DR", 
                           "4"="hotpink",    #"Giant cell_GPNMB", 
                           "5"="red3",    #"DCs_CD1C",
                           "6"="darkblue"    #"Mono_CD16"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), unlist(cellInfo$CellType))]
saveRDS(colVars, file="int/colVars.Rds")

# Initialize settings
library(SCENIC)
org="hgnc" # or mgi 
dbDir="cisTarget_databases" # RcisTarget databases location
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, nCores=1)# dbs=dbs, datasetTitle=myDatasetTitle,
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# change here if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

#Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions,minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)#基因过滤/选择 scenicOptions=scenicOptions, #minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
interestingGenes <- c("GPNMB", "TPRG1", "HLA-DRA")#Check if any genes are filtered out
interestingGenes[which(!interestingGenes %in% genesKept)]
#check expression matrix
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered_log <- log2(exprMat_filtered+1)
# Run GENIE3
runGenie3(exprMat_filtered_log, scenicOptions)  #very time consuming

#settings
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
### Build and score the GRN
# Toy run settings: scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Toy run settings:coexMethod=c("top5perTarget")
#other parameters can be changed: coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))

saveRDS(scenicOptions, file="int/scenicOptions after SCENIC_2.Rds") # To save status
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log) 
saveRDS(scenicOptions, file="int/scenicOptions after SCENIC_3.Rds") # To save status

#Optional steps:transform network activity to ON/OFF binary format
nPcs <- c(5,15,50) # For toy dataset: nPcs <- c(5)
#use different parameter to run t-SNE
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
#visualization (individual files in int/ ):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, varName="CellType", cex=.5)

# Using only "high-confidence" regulons (normally similar)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, varName="CellType", cex=.5)

#project AUC and TF expression to t-SNE
logMat <- exprMat_filtered_log # Better if it is logged/normalized
scenicOptions<-readRDS("int/scenicOptions after threshold set.Rds")
#use shiny to adjust the threshold
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)
#save the adjusted thresholds
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions after threshold set.Rds")

print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("MITF", "TFEC", "RUNX2","STAT3")],], plots="Expression")

runSCENIC_4_aucell_binarize(scenicOptions)

tsneAUC(scenicOptions, aucType="AUC") # choose settings
getSettings(scenicOptions,"defaultTsne/dims")
getSettings(scenicOptions,"aucell/tsne_perpl")
getSettings(scenicOptions,"seed")
scenicOptions@settings$seed <- 123

#create cell activity values
nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
pdf("int/AUC_tsne.pdf",width=10,height=10)
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
dev.off()




# if you have clustering information from Seurat try this
cellInfo<- data.frame(CellType=Idents(seuratobj))
nGene <- seuratobj@meta.data$nCount_RNA
nGene <- as.data.frame(nGene)
nUMI <- seuratobj@meta.data$nFeature_RNA
nUMI <- as.data.frame(nUMI)
cellInfo["nGene"] <- nGene
cellInfo["nUMI"] <- nUMI

plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

#visualize the SCENIC result in Seurat 
scRNA<-seuratobj
##inout original regulonAUC matrix
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)

AUCmatrix$colname<-rownames(AUCmatrix)
cellInfo$colname<-rownames(cellInfo)
AUCmatrix_label<-inner_join(AUCmatrix,cellInfo,by="colname")
head(BINmatrix_label)
rownames(AUCmatrix_label)<-AUCmatrix_label$colname
AUCmatrix_order<-AUCmatrix_label[order(AUCmatrix_label$CellType),]
AUCmatrix_map<-AUCmatrix_order[,-(325:327)]

RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
#scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

##inout biniary regulonAUC matrix
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
#manipulate data for pheatmap
library(dplyr)
BINmatrix$colname<-rownames(BINmatrix)
head(BINmatrix)
cellInfo$colname<-rownames(cellInfo)
BINmatrix_label<-inner_join(BINmatrix,cellInfo,by="colname")
head(BINmatrix_label)
rownames(BINmatrix_label)<-BINmatrix_label$colname
BINmatrix_order<-BINmatrix_label[order(BINmatrix_label$CellType),]
head(BINmatrix_order)
table(BINmatrix_order$CellType)
celltype = subset(BINmatrix_order,select = 'CellType')  
nGene=subset(BINmatrix_order,select = 'nGene') 
nUMI=subset(BINmatrix_order,select = 'nUMI')
cellInfo=subset(BINmatrix_order,select = c('CellType','nGene','nUMI'))
BINmatrix_map<-BINmatrix_order[,-(217:220)]
saveRDS(BINmatrix_order,"int/4.1_binaryRegulonActivity_ordered_bycelltype.Rds")
#BINmatrix_order<-readRDS("int/4.1_binaryRegulonActivity_ordered_bycelltype.Rds")

RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
#scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

## visualize AUC using Seurat
dir.create('scenic_seurat')
#FeaturePlot   
library(ggplot2)
p1 = FeaturePlot(scRNAauc, features='KDM5B_extended_1754g', label=F, cols=c("lightgray","red"),reduction = 'umap')
p2 = FeaturePlot(scRNAbin, features='KDM5B_extended_1754g', label=F, reduction = 'umap')
#p3 = DimPlot(scRNA, reduction = 'umap', group.by = "new.ident", label=T)
plotc =p1|p2  
ggsave('KDM5B_extended_1754g.pdf', plotc, width=10 ,height=4)

#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "MITF_244g", group.by="new.ident") + 
  theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features = "RUNX2_extended_77g", pt.size = 0, group.by="new.ident") + 
  theme(legend.position='none')
plotc = p2
ggsave('scenic_seurat/Ridge-Vln_RUNX2_extended_77g.png', plotc, width=10, height=8)

#heatmap of the cell activity
library(pheatmap)
install.packages("pheatmap")
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'CellType')

AUCmatrix_map_t <- data.frame(t(AUCmatrix_map))
BINmatrix_map_t <- data.frame(t(BINmatrix_map))
#select interest regulons
my.regulons <- c("CREB3L2_extended (87g)",	
                 "MAF (66g)",	
                 "MAFB (64g)",	
                 "NFIA (24g)",	
                 "HIF1A (11g)",	
                 "ETV5 (17g)",	
                 "MYC (240g)",	
                 "CEBPZ_extended (876g)",	
                 "YBX1_extended (25g)",	
                 "STAT3 (35g)",	
                 "ETS1 (22g)",	
                 "RUNX3_extended (39g)",	
                 "VDR (20g)",	
                 "STAT1 (192g)",	
                 "BCL11A_extended (519g)",	
                 "ARID5B (24g)",	
                 "NR1H3_extended (74g)",	
                 "MITF (244g)",	
                 "TFEC (163g)",	
                 "RUNX2 (53g)",	
                 "IRF1 (250g)",	
                 "ELF4 (156g)",	
                 "IKZF1 (100g)",	
                 "POU2F2_extended (89g)",	
                 "KLF13 (96g)"
)

myAUCmatrix <- AUCmatrix_map_t[rownames(AUCmatrix_map_t)%in%my.regulons,]
myBINmatrix <- BINmatrix_map_t[rownames(BINmatrix_map_t)%in%my.regulons,]
myAUCmatrix2<-as.data.frame(sapply(myAUCmatrix, as.numeric))#convert to numeric
rownames(myAUCmatrix2)<-rownames(myAUCmatrix)
#use the original AUC value to do heatmap, not a good choice, binmap is much better
rownames(celltype) <- colnames(AUCmatrix_map_t)
pheatmap(myAUCmatrix2, show_colnames=F, annotation_col=celltype,cluster_cols = F,cluster_rows =T,
         filename = 'scenic_seurat/myAUCmatrix_heatmap.png',
         width =15, height = 20)

#图use binary AUC value to do heatmap
rownames(celltype) <- colnames(BINmatrix_map_t)
rownames(nGene)<- colnames(BINmatrix_map_t)
rownames(nUMI)<- colnames(BINmatrix_map_t)
pheatmap(BINmatrix_map_t, show_colnames=F, annotation_col=celltype, cluster_cols =F,cluster_rows =T,
         filename = 'scenic_seurat/BINmatrix_heatmap_cor.pdf',
         color = colorRampPalette(colors = c("white","black"))(100),
         width =20, height =20)
myBINmatrix$TFs<-rownames(myBINmatrix)
myBINmatrix$TFs<-factor(myBINmatrix$TFs,levels=c("CREB3L2_extended (87g)",	
                                                 "MAF (66g)",	
                                                 "MAFB (64g)",	
                                                 "NFIA (24g)",	
                                                 "HIF1A (11g)",	
                                                 "ETV5 (17g)",	
                                                 "MYC (240g)",	
                                                 "CEBPZ_extended (876g)",	
                                                 "YBX1_extended (25g)",	
                                                 "STAT3 (35g)",	
                                                 "ETS1 (22g)",	
                                                 "RUNX3_extended (39g)",	
                                                 "VDR (20g)",	
                                                 "STAT1 (192g)",	
                                                 "BCL11A_extended (519g)",	
                                                 "ARID5B (24g)",	
                                                 "NR1H3_extended (74g)",	
                                                 "MITF (244g)",	
                                                 "TFEC (163g)",	
                                                 "RUNX2 (53g)",	
                                                 "IRF1 (250g)",	
                                                 "ELF4 (156g)",	
                                                 "IKZF1 (100g)",	
                                                 "POU2F2_extended (89g)",	
                                                 "KLF13 (96g)"
))
myBINmatrix_order<-myBINmatrix[order(myBINmatrix$TFs),]
myBINmatrix_order<-myBINmatrix_order[,-3969]









