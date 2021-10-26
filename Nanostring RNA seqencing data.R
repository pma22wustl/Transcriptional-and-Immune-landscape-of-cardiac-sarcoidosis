
#box plot for QC data
library(ggplot2)
library(tidyverse)
data <- read.csv("nanostring 3 group QC for box plot.csv")
p<-ggplot(data, aes(x=ScanLabel,y=TrimmedReads,color=ScanLabel))+ 
  geom_boxplot()+
  geom_point()+
  theme_classic()
p


#PCA plot
rm(list = ls())
library(dplyr)
library(ggplot2)
df <-read.csv("pca plot sarc vs. GCM&LM CD68+.csv")
df1 <- df[5:6]
pca_res <- prcomp(df1, scale. = TRUE)
p<- autoplot(pca_res, data = df, colour = 'Slide',frame = TRUE, frame.type = 'norm')
p<- p+geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2)+
  theme_minimal() 
p


#Venn plot
library(gplots)
venn_dat1<-read.csv("Sarc vs lymphocytic DE genes 3 grou_w.csv")
venn_up_SarcvsLymph<-venn_dat1$symbol[venn_dat1$change =="up"]
venn_dat2<-read.csv("Sarc vs GCM DE genes 3 grou_w.csv")
venn_up_SarcvsGCM<-venn_dat2$symbol[venn_dat2$change =="up"]
interup_Sarc<-intersect(venn_up_SarcvsLymph,venn_up_SarcvsGCM)
write.csv(interup_Sarc,"inter_up_in_Sarc_3 grou_w.csv")
venn(data = list(venn_up_SarcvsGCM,venn_up_SarcvsLymph))
interdown_Sarc<-venn_dat1$symbol[venn_dat1$change =="down"]
# two dimentional plot
install.packages("VennDiagram")
library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    "Sarc vs GCM"=venn_up_SarcvsGCM,
    "Sarc vs Lymph"=venn_up_SarcvsLymph
  ),
  filename = "Sarc_up.png", imagetype = "png",
  lwd = 2,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.6,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("cornflowerblue", "darkorchid1"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 20)
)
ggsave(venn.plot,filename = "venn_up_Sarc.png",path = "3 group download from web/")

#Venn plot
UpvsGCM<-read.csv("2groupDE_up Sarcoid vs GCM.csv")
UpvsLM<-read.csv("2groupDE_up Sarcoid vs LM.csv")
inner_up<-intersect(UpvsGCM$symbol,UpvsLM$symbol)
write.csv(inner_up,"inner_up genes in sarcoid vs GCM&LM.csv")
Venn_data <- list(UpvsGCM$symbol,UpvsLM$symbol)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
# Default plot
ggVennDiagram(Venn_data)
# Remove labels background color,Change category names,Change the gradient fill color
ggVennDiagram(
  Venn_data, label_alpha = 0,
  category.names = c("Sarcoidosis vs GCM","Sarcoidosis vs LM"))+
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")




#Volcano plot using DE genes
rm(list = ls())
library(dplyr)
library(ggplot2)
dat1 <-read.csv("sarco vs GCM/2group DE genes Sarc vs GCM.csv")
dat2<-read.csv("sarco vs lympho/2groupDEgenes GCMvsLM.csv")
dat<-full_join(dat1,dat2,by="Target.name")
write.csv(dat,"all DEgenes in 3group 2compare.csv")
dat<-read.csv("all DEgenes in 3group 2compare.csv")
colnames(dat)[1]="symbol"
colnames(dat)[2]="avg_logFC"
colnames(dat)[3]="p_val_adj"
dat<-dat[order(dat$avg_logFC),]
logFC_t=1
P.Value_t = 0.05
n1 = (dat$p_val_adj < P.Value_t)&(dat$avg_logFC < -logFC_t)
n2 = (dat$p_val_adj < P.Value_t)&(dat$avg_logFC > logFC_t)
change = ifelse(n1,"down",ifelse(n2,"up","stable"))
dat <- mutate(dat,change)
table(dat$change)
deg_down<-dat[dat$change=="down",]
write.csv(deg_down,"2groupDE_down Sarcoid vs GCM.csv")
write.csv(deg_down,"2groupDE_down Sarcoid vs LM.csv")
deg_up<-dat[dat$change=="up",]
write.csv(deg_up,"2groupDE_up Sarcoid vs GCM.csv")
write.csv(deg_up,"2groupDE_up Sarcoid vs LM.csv")
#selected genes to label
colnames(dat)
for_label <- dat%>% 
  filter(symbol %in% c("HLA-DRA", "STAT1","MMP9","LYZ", "LAMP3",
                       "CHIT1","LYZ","CHI3L1","FBP1","HLA-DRA","HLA-DPB1","FTH1","CXCL9",
                       "HLA-DPA1","STAT1","MMP9","HLA-DQB1","HLA-DQA2","CCL5","VIM","HLA-DRB1",
                       "LAMP1","S100A11","ITGAX","GPNMB","FKBP5")) 
#label head(10) genes
x1 = head(dat[dat$change=="down",] ,10)
#write.csv(x1,"sarc vs lympho down.genes.csv")
x2 = tail(dat[dat$change=="up",],15)
#write.csv(x2,"sarc vs lympho up.genes.csv")
for_label = rbind(x1,x2)
#data visualization with ggplot2
library(ggplot2)
p <- ggplot(data = dat, 
            aes(x = avg_logFC, 
                y = -log10(p_val_adj))) +
  geom_point(alpha=0.7, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("red", "grey","blue"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
#plot with gene label
volcano_plot <- p +
  geom_point(size = 3.0, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",max.overlaps = Inf
  )
volcano_plot
volcano_plot+xlim(-2,4)
ggsave("sarcoid vs lympho volcano top 10.png",width = 10,height = 6, path = "sarco vs lympho/")
ggsave("sarcoid vs lympho volcano filtered gene.png",width = 10,height = 6, path = "sarco vs lympho/")




#heatmap for nanostring data
rm(list = ls())
exp<-read.csv("heatmap/CD3 CS vs. LM raw data.csv",header =FALSE)
#remove first col and row
exp4<-exp[-1,]
exp2<-exp[,-1]
exp3<-exp2[-1,]
#transform dataframe into matrix
ex<-as.matrix(exp3)
#transform string into numeric
ex=apply(ex,2,as.numeric)
ex=log2(ex)#log2 transform normalization
#set colnames and rownames
colnames(ex)=exp2[1,]
rownames(ex)=exp4[,1]
#remove all 0 valuse
ex[ex==0]<-NA
ex=na.omit(ex)##reomve rows with "NA"
head(ex)

# produce grouplist
pl<-read.csv("heatmap/pl 2 group cd3.csv",header = FALSE)
colnames(pl)=c("sample","group")
group_list<-pl$group
print(group_list)
#set the levels
group_list = factor(group_list,
                    levels = c("Sarcoid","Lymphocytic"))
#heatmap for the first 1000 DEGs 
cg=names(tail(sort(apply(ex,1,sd)),1000))#查看前1000个sd值最大（组间差异最大）的基因
n=ex[cg,]
head(n)
#generate heatmap
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
#install.packages('pheatmap')
library(pheatmap)
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         cluster_cols = F,
         annotation_col=annotation_col,
         scale = "row")

library(ggplot2)
ggsave("Heatmap_Sarc_vs_LM.png",width = 15,height =20,path = "heatmap/")

#read the selest genes you want to present
cg<-read.csv("heatmap/heatmap for proposal cd3.csv",header=F)
class(cg)
#check these genes
cg1<-cg[,1]
print(cg1)
class(cg1)
#order the grouplist 
group_list = group_list[order(group_list)]
# gene label for heatmap
n<-ex[cg1,]
print(cg1)
dim(n)
head(n)

#generate heatmap
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
#install.packages("ggplotify")
library(ggplotify)
library(RColorBrewer)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                "RdBu")))(100),
                                   show_rownames = T,
                                   scale = "row",
                                   cluster_cols = F,
                                   cluster_rows = F,
                                   annotation_col=annotation_col)) 
heatmap_plot
ggsave(" Sarc_vs LM_selectheatmap.png",width=10,height=15,path = "heatmap/")


#Gene enrichment analysis and KEGG pathway analysis using clusterProfiler
library(ggplot2)
library(clusterProfiler)
install.packages("BiocManager")
library(BiocManager)
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
install.packages("clusterProfiler")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(bitr)
s2e <- bitr(cg, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

s2e_all <- bitr(deg_Lymphotic_vs_Sarc$symbol, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)#人类
#other species:  http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

deg_Lymphotic_vs_Sarc <- inner_join(deg_Lymphotic_vs_Sarc,s2e_all,by=c("symbol"="SYMBOL"))
write.csv(deg_inner,"deg_inner_3group_w.csv")
write.csv(deg_innner,"DEgens for KeGG.csv")
library(clusterProfiler)
library(dplyr)
library(ggplot2)
source("kegg_plot_function.R")

#1.GO database analysis ----

#(1)data input
deg<-deg_innner
gene_up=s2e$ENTREZID
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#(2)Go analysis
if(T){
  #CC
  ego_CCup <- enrichGO(gene = gene_up, #gene_diff,  or gene_down
                       OrgDb= org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.01,
                       readable = TRUE)
  #BP
  ego_BPup <- enrichGO(gene = gene_up, #gene_diff,  or gene_down
                       OrgDb= org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.01,
                       readable = TRUE)
  #MF：
  ego_MFup <- enrichGO(gene = gene_up, #gene_diff,  or gene_down
                       OrgDb= org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.01,
                       readable = TRUE)
  save(ego_CCup,ego_BPup,ego_MFup,file = "GO analysis intergenes_3group_w.Rdata")
}
load(file = "GO analysis intergenes_3group_w.Rdata")

#(3)data visualization
#bar plot
barplot(ego_MFup)
barplot
ggsave("barplot_up_BPinSarc.png",path = "3 group download from web/")
#bubble plot
dotplot(ego_MFup)
p1<-dotplot(ego_CCup,x = "GeneRatio")
ggsave("dotplotCCupinCcl17high.png",width = 8,height = 8,path = "3 group/" )
p2<-dotplot(ego_BPup,x = "GeneRatio")
ggsave("dotplotBPup.png",width = 8,height = 8,path = "3 group/" )
p3<-dotplot(ego_MFup,x = "GeneRatio")
ggsave("dotplotMFupinSarc.png",width = 8,height = 8,path = "3 group/" )
library(patchwork)
(p1 + p2 +p3)+ plot_annotation(tag_levels = "A")
ggsave(filename = "Combined cc_bp-mf.png",width = 30,height = 10,path = "3 group download from web/")


#2.KEGG pathway analysis
#data inout
library(dplyr)
logFC_t=1
P.Value_t = 0.05
k1 = (deg_inner$logFC < -logFC_t)
k2 = (deg_inner$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg_inner <- mutate(deg_inner,change)
table(deg_inner$change)
deg_inner<-deg_inner[order(deg$logFC),]
save(deg_inner,file="deg_inner.Rdata")
gene_up = deg_inner[deg_inner$change == 'up','ENTREZID'] 
gene_down = deg_inner[deg_inner$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg_Lymphotic_vs_Sarc[,'ENTREZID']

#pathway enrichment using DEGs
if(T){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9)
  save(kk.diff,kk.down,kk.up,file = "kegg_inner_3group_w.Rdata")
}
load("Ccl17kegg.Rdata")
#(3)subset the result from database
kegg_up_dt <- kk.up@result

#(4)screen by p value
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% 
  mutate(group=-1) 
up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)
#(5)visualization
g_kegg <- kegg_plot(down_kegg,up_kegg)
g_kegg 
#+ scale_y_continuous(labels = c(15,10,5,0,5,10))
ggsave(g_kegg,filename = 'INNER_kegg_up_down.png',width=20,height=20,path = "3 group download from web/")

#Enrichment analysis using gsea 
#(1)check the example data
data(geneList, package="DOSE")
#(2)format transformation
geneList=deg_inner$logFC
names(geneList)=deg_inner$ENTREZID
head(geneList)
geneList=sort(geneList,decreasing = T)
#(3)enrichment analysis
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
save(down_kegg,up_kegg,file="GSEA Ccl17.Rdata")
#(4)visualization
gse_d<-kegg_plot(up_kegg,down_kegg)
gse_d
gse_d+ scale_y_continuous(labels = c(3,2,1,0,1,2,3))
ggsave(filename = 'GSEA.png',width=15,height=10,path = "3 group/")





