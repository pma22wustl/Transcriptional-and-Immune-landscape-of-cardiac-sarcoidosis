#cellcall packages
if (!requireNamespace ("BiocManager", quietly = TRUE)) 
  install.packages ("BiocManager") 
BiocManager::install ("clusterProfiler") 
BiocManager::install("enrichplot",force = TRUE)
BiocManager::install("DOSE",force = TRUE)
library(DOSE)
library(devtools)
devtools::install_github("ShellyCoder/cellcall")
library(cellcall)
install.packages('Seurat')
library(Seurat)

seuratobj <- readRDS("renamed t cells 5group.rds")
seuratobj<- RenameIdents(seuratobj, `CD4_Treg`="CD4",`CD4_naive`="CD4",`CD4_trm`="CD4")
seuratobj1 <- readRDS("Renamed_Myeloid_6group.rds")
seuratobj1<- RenameIdents(seuratobj1, `Mac_HLA-DR`="Mac inf",`Mac_VCAN`="Mac inf",`Mac_res`="Mac res")

DimPlot(seuratobj1, dims=c(1,2), reduction = "umap",label=TRUE) #split.by = "group"
seuratobj.combined <- merge(seuratobj, y = seuratobj1, add.cell.ids = c("T", "myeloid"), project = "combined")
seuratobj2<-seuratobj
seuratobj<-seuratobj.combined

#Part I: Data input & processing and initialization of CellChat object
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)

data.input  <- data.frame(seuratobj@assays$RNA@counts)
celltype = data.frame(group =seuratobj@active.ident,row.names = names(seuratobj@active.ident))
colnames(data.input) <- paste(1:length(celltype$group), celltype$group, sep = "_")
saveRDS(data.input,"input data counts with cellinfo sarcoid combined.rds")
mt<-readRDS("input data counts with cellinfo sarcoid combined.rds")
#Create object

mt <- CreateNichConObject(data=data.input, min.feature = 0,
                          names.field = 2,
                          names.delim = "_",
                          source ="TPM",# "fullLength",#"UMI",, or .
                          scale.factor = 10^6,
                          Org ="Homo sapiens", # mouse "Mus musculus"
                          project = "Microenvironment")
#Infer the cell-cell communication score
#The communication score of an L-R interaction between cell types 
#is evaluated by integrating the L2- norm of the L-R interaction and the activity score of the downstream TFs
mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.9,#default=0.05
                        CorValue = 0.05, #default=0.1
                        topTargetCor=1,
                        p.adjust = 0.9,#default=0.05
                        use.type="median",
                        probs = 0.7, #default=0.9
                        method="weighted",
                        IS_core = TRUE,
                        Org ='Homo sapiens')  #"Mus musculus"

#Pathway activity analysis
n<-mt@data$expr_l_r_log2_scale
pathway.hyper.list<-lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")#"Mus musculus"
  return(tmp)
})#sender cellA and receiver cellB, user can set cella_cellb="A-B".

myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))#between sender cellA and receiver cellB: cella_cellb="A-B"
p <- plotBubble(myPub.df)
p

#Circle plot is adopted to present the global cell-cell communications between cell types.
#Setting the color and name of each cell type:
cell_color <- data.frame(color=c('#e31a1c','#1f78b4',"forestgreen","darkorange","magenta4",
                                 '#e78ac3','black',"yellow"), stringsAsFactors = FALSE)
rownames(cell_color) <- c( "CD8", "CD4", "NK","Mac res", "Mac inf","Giant cell","DCs","Mono")

#present the global cell-cell communications
order.vector<-c("CD8", "CD4", "NK","Mac res", "Mac inf","Giant cell","DCs","Mono")
ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 2,
                order.vector=order.vector,
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)

#Pheatmap plot is adopted to present the detailed communication scores
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 12,angle_col = "90", 	
             main="score")

#Sankey plot is adopted to present the detailed L-R-TF axis
mt <- LR2TF(object = mt, sender_cell=c("Giant cell"), recevier_cell="Mac inf",
            slot="expr_l_r_log2_scale", org="Homo sapiens") #)
mt <- LR2TF(object = mt, sender_cell="Cxcl9 Gbp2b Mac", recevier_cell="Cxcl9 Gbp2b Mac",
            slot="expr_l_r_log2_scale", org="Mus musculus") #"Homo sapiens")
mt <- LR2TF(object = mt, sender_cell="Cxcl9 Gbp2b Mac", recevier_cell="CD8 T",
            slot="expr_l_r_log2_scale", org="Mus musculus") #"Homo sapiens")
mt <- LR2TF(object = mt, sender_cell=c("CD4","CD8"), recevier_cell="Giant cell",
            slot="expr_l_r_log2_scale", org="Homo sapiens") #"Homo sapiens")
head(mt@reductions$sankey)
p<-mt@reductions$sankey
p
if(!require(networkD3)){
  BiocManager::install("networkD3")
}

sank <- LRT.Dimplot(mt, fontSize = 8, nodeWidth = 30, height = NULL, width = 1200, 		 
                    sinksRight=FALSE, DIY.color = FALSE)
networkD3::saveNetwork(sank, "~/ST-SSC_full.html")

library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0) ## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  ## transform weight1 and weight2 to one value (weight)
head(tmp.df)

## set the color of node in sankey graph
mycol.vector = c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
elments.num <-  length(unique(tmp.df$Ligand))
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], 
             isGrandSon = TRUE, nudge_x = nudge_x, font.size = 2, boder.col="white", 			
             set_alpha = 0.8)
ssc.tf <- names(mt@data[["gsea.list"]][["Giant cell"]]@geneSets)
ssc.tf

## gsea object
egmt <- mt@data[["gsea.list"]][["Giant cell"]]@result
## filter TF
egmt.df <- data.frame(egmt)
head(egmt.df[,1:6])
flag.index <- which(egmt.df$p.adjust < 0.05)




getGSEAplot(gsea.list=mt@data[["gsea.list"]], geneSetID=c("HIF1A", "PARP1", "NFKBIA","STAT2"), 
            myCelltype=mt@data[["gsea.list"]][["Giant cell"]], fc.list=mt@data[["fc.list"]],  
            selectedGeneID = mt@data[["gsea.list"]][["Giant cell"]]@geneSets[["PARP1"]][1:10],
            mycol = NULL)

ridgeplot.DIY(x=egmt, fill="pvalue", showCategory=flag.index, core_enrichment = T,
              orderBy = "NES", decreasing = FALSE)

#cell-cell interatcion
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap",force = TRUE)
devtools::install_github("jokergoo/circlize")
library(circlize)
#install.packages('NMF')
library(NMF)
#remotes::install_github("sqjin/CellChat",force = TRUE)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(patchwork)
library(tidyverse)

seuratobj <- readRDS("renamed t cells 5group.rds")
seuratobj1 <- readRDS("Renamed_Myeloid_6group.rds")
DimPlot(seuratobj1, dims=c(1,2), reduction = "umap",label=TRUE) #split.by = "group"
seuratobj<- RenameIdents(seuratobj, `CD4_trm`="CD4",`CD4_naive`="CD4",`CD4_Treg`="CD4")
seuratobj1<- RenameIdents(seuratobj1, `Mac_VCAN`="Mac_HLA-DR")
#remove or subset clusters
seuratobj <- subset(x= seuratobj, ident= c("myeloid","NK&T"), invert = F)#T=remove ident
#Idents(seuratobj) <- seuratobj$integrated_snn_res.0.9

#Part I: Data input & processing and initialization of CellChat object
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)
seuratobj.combined <- merge(seuratobj, y = seuratobj1, add.cell.ids = c("T", "myeloid"), project = "combined")
seuratobj<-seuratobj.combined

data.input  <- seuratobj@assays$RNA@data
identity = data.frame(group =seuratobj@active.ident,row.names = names(seuratobj@active.ident)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat <- createCellChat(data.input)
summary(cellchat)

#add metadata to CellChat object
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


#CellChatDB.human,CellChatDB.mouse import the L-R pairs database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
CellChatDB.use <- CellChatDB # use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)  

#Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Extract the inferred cellular communication network as a data frame
#access the inferred cell-cell communications of interest
df.net <- subsetCommunication(cellchat)#slot.name = "netP", to access the inferred communications at the level of signaling pathways
write.csv(df.net,"cell-cell communications myeloid&T with less celltype.csv")
levels(cellchat@idents)
df.net <- subsetCommunication(cellchat, sources.use = c(4,5,6), targets.use = 1)#get cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5
df.net <- subsetCommunication(cellchat, signaling = c("CD45", "CD6")) #cell-cell communications mediated by signaling WNT and TGFb.


#inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.
cellchat <- computeCommunProbPathway(cellchat)#Infer the cell-cell communication at a signaling pathway level
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)#sources.use= and targets.use=: calculate the aggregated network among a subset of cell groups

# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)#creating a simple multi-paneled plot
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#examine the signaling sent from each cell group. control the parameter edge.weight.max to compare edge weights between differet networks.
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Part III: Visualization of cell-cell communication network
#Hierarchy plot: define vertex.receiver, the target cell groups,solid and open circles represent source and target, respectively
#Chord diagram: netVisual_chord_cell(cell groups), netVisual_chord_gene(pathway, ligand or receptor)
netVisual_chord_cell("TGFb")
saveRDS(cellchat,"myeloid&T with less celltype.rds")
cellchat<-readRDS("myeloid&T with less celltype.rds")
#check the result
cellchat@netP$pathways
LRsig<-cellchat@LR$LRsig
write.csv(LRsig,"sig L-R pairs.csv")
#check the cell clusters
levels(cellchat@idents) 
#visualize the signal pathway
pathways.show <- cellchat@netP$pathways
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# Chord diagram
group.cellType <- c(rep("T", 4),rep("myeloid", 6),"T") # grouping cell clusters into T and myeloid cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

#contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
#extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR<- extractEnrichedLR(cellchat, signaling =pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[26,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

saveRDS(cellchat,"cellchat myeloid&T combined after calculation.rds")
cellchat<-readRDS("cellchat myeloid&T combined after calculation.rds")
#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#show (L-R pairs) from some cell groups to other cell groups using netVisual_bubble
netVisual_bubble(cellchat, sources.use =4, targets.use =c(1:8), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:8), signaling = c("TGFb","PTPRM"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MHC-II","PTPRM"))
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2), targets.use = c(4:6), slot.name = "netP", legend.pos.x = 10)
netVisual_chord_gene(cellchat, sources.use = c(5:10), targets.use =c(1,2,3,4) , slot.name = "netP", legend.pos.x = 10)

#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "PTPRM")#enriched.only = FALSE show all genes
plotGeneExpression(cellchat, signaling ="TGFb")
# extract the signaling genes related to the inferred L-R pairs or signaling pathway
gene_expression<-extractEnrichedLR(cellchat,signaling = pathways.show.all)
gene_expression

#Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(4,2))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTPRM","LAMININ","CD99"))
gg1 + gg2
gg2

# which signals contributing most to outgoing or incoming signaling of certain cell groups.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1
ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PTPRM","LAMININ","CD99"))
ht
#Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
#run selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4 #Cophenetic values begin to drop when the number of incoming patterns is 4.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
install.packages("Miniconda")
install.packages("Python")
installed.packages('umap-learn')
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

saveRDS(cellchat,"cellchat myeloid&T combined after calculation.rds")
sessionInfo()


