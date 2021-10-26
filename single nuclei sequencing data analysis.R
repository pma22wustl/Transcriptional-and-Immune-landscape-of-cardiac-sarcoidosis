

#TFs wgene regulation network production
library(tidyverse)
library(igraph)
library(ggraph)
data <- read_csv("Step2regulonTargetsInfo.csv")
TFs <- c('MITF', 'TFEC', 'ARID5B', 'NR1H3',"KDM5B")
marker<-read_csv("giant cell markers_fc2_0.5.csv")
marker<-marker$symbol
data2 <- data %>%
  filter(TF %in% TFs,
         gene %in% marker)

# number of packages to display
number_of_genes <- 200
top_gene <- data2 %>%
  group_by(gene)%>%
  summarize(total_nMotifs = sum(nMotifs)) %>%
  arrange(desc(total_nMotifs)) %>%
  select(gene, total_nMotifs)
# all popular languages per package
top_TF_per_gene <- data %>%
  filter(
    gene %in% top_gene$gene,
    TF %in% TFs) %>%
  arrange(gene, desc(nMotifs)) %>%
  group_by(gene) %>%
  mutate(
    main = row_number() == 1, # main language of package should be opaque
    total_nMotifs = sum(nMotifs)
  ) %>%
  ungroup() %>%
  select(TF, gene, nMotifs, total_nMotifs, main)

# only following languages found in given packages
(top_TFs <- top_TF_per_gene %>%
    pull(TF) %>%
    unique %>%
    sort)
top_TF_colors <- c(
  "forestgreen",
  "darkorange",     
  "hotpink",    
  "red3",
  "darkblue",
  "#000000")
names(top_TF_colors) <- c(
  "ARID5B", "KDM5B",  "MITF", "NR1H3","TFEC",  
  'All')

edges1 <- top_TF_per_gene %>%
  transmute(from = TF, to = gene, total_nMotifs = nMotifs, main)

edges2 <- top_TF_per_gene %>%
  count(TF, wt = nMotifs, name = 'total_nMotifs') %>%
  transmute(
    from = '',
    to = TF,
    total_nMotifs,
    main = TRUE )
edges <- bind_rows(edges1, edges2)

vertices1 <- top_TF_per_gene %>%
  filter(main)%>%
  transmute(
    node = gene, TF, total_nMotifs, level = 1)

vertices2 <- edges2 %>%
  transmute(
    node = to, TF = to, total_nMotifs, level = 2)

vertices3 <- tibble(
  node = '', TF = NA, total_nMotifs = 0, level = 3)
vertices <- bind_rows(vertices1, vertices2, vertices3) %>%
  mutate(
    radius = total_nMotifs**(1.8), # scaling circles
    TF = factor(TF, names(top_TF_colors))
  ) %>%
  arrange(level, TF, node)

graph <- graph_from_data_frame(edges, vertices = vertices)
# create custom layout by updating existing circle layout
layout <- create_layout(graph, layout = 'circle')
outer_circle <- layout %>%
  filter(level == 1) %>%
  mutate(TF = factor(TF, names(top_TF_colors))) %>%
  arrange(TF, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / number_of_genes * 2 * pi),
    y = sin((row_number() - 1) / number_of_pkgs * 2 * pi))

# positioning circle centers manually by specifying polar coords
angles <- c(3, 43,255, 350, 190, 0)
radii <- c(0.8, 0.5, 0.6, 0.6, 0.38, 0)
centers <- tibble(
  x = radii * cos(angles / 180 * pi),
  y = radii * sin(angles / 180 * pi)
)
inner_circle <- bind_cols(centers, select(filter(layout, level != 1), -x, -y))

layout[] <- bind_rows(outer_circle, inner_circle) %>%
  arrange(.ggraph.index)

ggraph(layout) +
  geom_edge_diagonal(
    aes(edge_color = node1.TF, edge_alpha = as.factor(main)),
    edge_width = 0.3, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius, color = TF),
    alpha = 0.6, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% top_TFs)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = name,
      filter = name %in% top_TFs
    ),
    size = 6, hjust = 0.5, family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y - 0.045,
      label = ifelse(
        total_nMotifs > 1000,
        format(total_nMotifs, big.mark = ','),
        total_nMotifs
      ),
      filter = name %in% top_TFs
    ),
    size = 3, hjust = 0.5, family = 'Oswald'
  ) +
  scale_edge_color_manual(values = top_TF_colors) +
  scale_color_manual(values = top_TF_colors) +
  scale_size_area(max_size = 150) +
  scale_edge_alpha_manual(values = c(0.15, 1)) +
  coord_fixed() +
  labs(
    title = 'TFs of giant cells',
    subtitle = '',
    caption = '#tidytuesday 46|2019 spren9er'
  ) +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 20, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )

ggsave(
  'images/tidytuesday_201946_cran_packages.png',
  width = 12, height = 12.5, dpi = 300)
  
  
  
#Combined Violin Plots for marker genes in each celltype
  modify_vlnplot<- function(obj, 
                            feature, 
                            pt.size = 0, 
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {
    p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
      xlab("") + ylab(feature) + ggtitle("") + 
      theme(legend.position = "none", 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.y = element_text(size = rel(1), angle = 0), 
            axis.text.y = element_text(size = rel(1)), 
            plot.margin = plot.margin,
            axis.text=element_text(size=15, angle = 45)) 
    return(p)
  }
  
  ## extract the max value of the y axis
  extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  
  ## main function
  StackedVlnPlot<- function(obj, features,
                            pt.size =5, 
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {
    
    plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
    
    # Add back x-axis title to bottom plot. patchwork is going to support this?
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
      theme(axis.text.x=element_text(), axis.ticks.x = element_line())
    
    # change the y-axis tick to only max value 
    ymaxs<- purrr::map_dbl(plot_list, extract_max)
    plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                              scale_y_continuous(breaks = c(y)) + 
                              expand_limits(y = y))
    
    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
  }
  
  features<- c('COL4A3','F13A1','CD247','RGS5',
               'MYH11','MYOM2','VWF','NRG1','KIT','ADIPOQ','NRXN1',
               'PIEZO2')
  png(file = "Violin Plots modified.png", 
      width = 10*100,
      height = 15*100,
      res = 100)
  StackedVlnPlot(obj = seuratobj, features = features)
  dev.off()
  

  
  
#single cell data analysis by Seurat packages
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  seuratobj <- readRDS("cardiac sarcoidosis snRNA seq.rds")
  seuratobj <- subset(x= seuratobj, ident= c("0","2","5","8"), invert = F)#remove ident
  saveRDS(seuratobj, "myeloid_UMAP.ident.rds")
  
  seuratobj <- readRDS("Renamed_Myeloid_UMAP.new.ident.rds")
  DefaultAssay(seuratobj) <- "RNA"
  seuratobj <- NormalizeData(seuratobj)#,normalization.method = "LogNormalize", scale.factor = 10000
  seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seuratobj), 10)
  plot1 <- VariableFeaturePlot(seuratobj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  #scaling data
  all.genes <- rownames(seuratobj)
  seuratobj <- ScaleData(seuratobj, features = all.genes)
  DefaultAssay(seuratobj) <- "integrated"
  seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object =seuratobj))
  # Examine and visualize PCA results a few different ways
  print(seuratobj[["pca"]], dims = 1:15, nfeatures = 5)
  VizDimLoadings(seuratobj, dims = 1:5, reduction = "pca")
  DimPlot(seuratobj, reduction = "pca")
  DimHeatmap(seuratobj, dims = 1:10, cells = 500, balanced = TRUE)
  
  #Determine the ‘dimensionality’ of the dataset
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
  seuratobj <- JackStraw(seuratobj, num.replicate = 100)
  seuratobj <- ScoreJackStraw(seuratobj, dims = 1:20)
  JackStrawPlot(seuratobj, dims = 1:20)
  ElbowPlot(seuratobj)
  #Cluster the cells
  seuratobj <- FindNeighbors(seuratobj, dims = 1:30)
  seuratobj <- FindClusters(seuratobj, resolution = c(0.25,0.3,0.35,0.4,0.45,0.5,0.6))
  # Look at cluster IDs of the first 5 cells
  head(Idents(seuratobj), 5)
  #Run non-linear dimensional reduction (UMAP/tSNE)
  DefaultAssay(seuratobj) <- "integrated"
  seuratobj <- RunUMAP(seuratobj, dims = 1:30)
  DimPlot(seuratobj, reduction = "umap")
  Idents(seuratobj)<-seuratobj$integrated_snn_res.0.35
  DimPlot(seuratobj, reduction = "umap")
  
  
  #for integration data
  # Run the standard workflow for visualization and clustering
  seuratobj <- readRDS("myeloid_UMAP.ident.rds")
  DefaultAssay(seuratobj) <- "integrated"
  seuratobj <- ScaleData(seuratobj, verbose = FALSE)
  seuratobj <- RunPCA(seuratobj, npcs = 30, verbose = FALSE)
  seuratobj <- RunUMAP(seuratobj, reduction = "pca", dims = 1:30)
  seuratobj <- FindNeighbors(seuratobj, reduction = "pca", dims = 1:30)
  seuratobj <- FindClusters(seuratobj, resolution = c(0.2,0.3,0.32,0.35,0.38,0.4,0.5,0.6))
  Idents(seuratobj) <- seuratobj$integrated_snn_res.0.3
  DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 1.0)
  ggsave("res0.8_pt1.png")
  #rename each cellytype
  seuratobj<- RenameIdents(seuratobj, `0` = "Mac_res", `1` = "Mac_res", `2` = "Mac_VCAN", 
                           `3` = "Mac_HLA-DR",`4` = "Giant cell_GPNMB",  `5` = "DCs_CD1C",`6` = "Mono_CD16")

  DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE) #split.by = "group"
  saveRDS(seuratobj,"renamed sarcoid only_UMAP.rds")
  VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          group.by = "orig.ident" ,ncol = 3)+ #split.by="orig.ident"
    theme(axis.text.x = element_text(angle = 90))
  
  
  #for UMAPs
  DefaultAssay(seuratobj) <- "integrated"
  Idents(seuratobj) <- seuratobj$integrated_snn_res.0.5
  DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 0.5)
  
  # remove cluster
  #seuratobj <- SubsetData(object = seuratobj, ident.remove = c("8","9"))
  #DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 2.0)
  library(ggplot2)
  ggsave("res0.4 split.png")
  DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 1.5)
  ggsave("res0.4.png")
  
  
  seuratobj <- subset(x= seuratobj, ident= c("Mac_HLA-DR","Giant cell_GPNMB"), invert = F)#remove ident
  DefaultAssay(seuratobj) <- "RNA"
  seuratobj <- NormalizeData(seuratobj)
  seuratobj <- ScaleData(seuratobj)
  cluster.markers <- FindAllMarkers(seuratobj, only.pos =T, min.pct = 0.25, logfc.threshold = 0.25)
  top5<-cluster.markers %>%
    group_by(cluster) %>%
    top_n(n =30, wt = avg_log2FC)
  DoHeatmap(seuratobj, features = top5$gene) + NoLegend()
  write.csv(top5, "top20 genes sarcoid_all by celltype.csv")
  write.csv(cluster.markers,"all marker genes fibro0.3.csv")
  write.csv(cluster.markers,"DEGs SYTL3 vs Giant cells neg+pos.csv")
  saveRDS(cluster.markers,"all marker genes sarcoid_all by celltype.rds")
  
  DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=T,pt.size = 0.8)
  saveRDS(seuratobj,"Renamed_Myeloid_UMAP.new.ident.rds")
  library(ggplot2)
  ggsave("res0.4 new total with name.png")
  DimPlot(seuratobj, dims=c(1,2), reduction = "umap",split.by="group",label=FALSE,pt.size = 1.0)
  ggsave("res0.4 new split with name pt1.0-2 .png",path = "exported figures/")
  
  #cluster cell population in each group
  library(scales)
  table(Idents(seuratobj))
  table(Idents(seuratobj),seuratobj$group)
  prop.table(table(Idents(seuratobj)))
  cell.prop<-as.data.frame(prop.table(table(Idents(seuratobj), seuratobj$group)))
  colnames(cell.prop)<-c("cluster","group","proportion")
  ggplot(cell.prop,aes(group,proportion,fill=cluster))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("")+
    scale_y_continuous(labels = percent)+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'))+
    
    guides(fill=guide_legend(title=NULL))
  ggsave("percent.png")
  
  # cluster composition
  table(Idents(seuratobj),seuratobj$group)
  cell.count<- as.data.frame(table(Idents(seuratobj),seuratobj$group))
  colnames(cell.count)<-c("cluster","group","count")
  write.csv(cell.count,"cell.count.csv")
  
  DF_Count2<-cell.count %>% group_by(cluster)%>%
    mutate(freq = count / sum(count)*100)
  
  png(file = "Myeloid0.4_Cell Composition Percent.png", 
      width = 15*100,
      height = 10*100,
      res = 100)
  ggplot(DF_Count2, aes(x =cluster, y = freq, fill = group))+
    geom_col() +
    labs(x="Groups",y="Percentage",fill="Clusters")+
    #geom_text(aes(label = paste(round(freq, 2),"%")),position = position_stack(vjust = 0.5))+
    theme_classic()
  dev.off()
  
  #switch assay to RNA Expression analysis
  
  DefaultAssay(seuratobj) <- "RNA"
  seuratobj <- NormalizeData(seuratobj)
  seuratobj <- ScaleData(seuratobj)
  
  #plot with original data
  #vlnplot
  VlnPlot(seuratobj, features = c( "GPNMB","FBP1","CHIT1","CHI3L1", "MMP9","TPRG1"),split.by="group",split.plot=TRUE,pt.size = 0,ncol = 2)
  library(ggplot2)
  ggsave(filename ="genes vlnplot3.png",width = 15,height = 20)
  
  #Featureplot 
  FeaturePlot(seuratobj, features = c("GPNMB","FBP1","CHIT1","CHI3L1", "MMP9","TPRG1"),pt.size = 0.2,label = FALSE,
              blend = FALSE,cols = c("lightgrey","red"))
  FeaturePlot(seuratobj, features = c("CD68", "MRC1", "F13A1", "C1QA"), pt.size = 0.2,label = FALSE,
              blend = FALSE,cols = c("lightgrey","red"))
  ggsave("Ccl7.png")
  #Featureplot two genes in one plot
  FeaturePlot(seuratobj, features = c("S100A6","CD14"),pt.size = 0.5,label = F,
              blend = TRUE )

  
  #density plot with Nebulosa
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Nebulosa")
  library("Nebulosa")
  plot_density(seuratobj, c("MITF","KDM5B","TFEC","NR1H3"))
  plot_density(seuratobj, c("PECAM1"))
  #p3 + plot_layout(ncol = 1)#show the figure in one colomn
  
  
  #boxplot entropy and pseudotime 
  data<-read.csv("myeloid_palantir_meta_data.csv")
  data<-data[,7:8]
  data$ClusterName<-as.factor(data$ClusterNam)
  
  p<-ggplot(data, aes(x=ClusterName, y=entropy, fill=ClusterName)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, position = position_jitterdodge())+
    theme_classic()
   
  p<-ggplot(data, aes(x=ClusterName, y=entropy, fill=ClusterName)) +
    geom_boxplot()+
    ylim(0, 1.5)+
    theme_classic()
  p
  
  #waterfall plot for pseudotime
  data<-read.csv("myeloid_palantir_meta_data -new.csv")
  data<-read.csv("myeloid_palantir_meta_data.csv")
  data<-data[,c(1,6:8)]
  data$ClusterName<-as.factor(data$ClusterName)
  data1<-data[data$ClusterName!=1,]
  data1<-data1[data1$ClusterName!=5,]
  p<-ggplot(data, aes(x=pseudotime, y=entropy, color=ClusterName)) +
    geom_point(size=1.5)+
    theme_classic()
  p
  p+scale_color_manual(values=c("#F8766D","#7CAE00","#00BFC4","#C77CFF"))
  p+scale_color_manual(values=c("lightgray","lightgray","#00BFC4","#C77CFF"))
  p+scale_color_manual(values=c("#F8766D","lightgray","lightgray","#C77CFF"))
  p+scale_color_manual(values=c("lightgray","#7CAE00","lightgray","#C77CFF"))
  
  #show default colors in ggplot2
  library(scales)
  show_col(hue_pal()(4))
  hue_pal()(6)
  
  
  #Plotting multiple genes in the same FeaturePlot
  object<-seuratobj
  marker_gene_list <- list(c("CD68", "MRC1","C1QA"))
  marker_gene_list <- list(c("F13A1","CD163L1","SCN9A","IL2RA","IL12RB2"))
  object <- AddModuleScore(object, features = marker_gene_list, name = "Macrophage")
  object <- AddModuleScore(object, features = marker_gene_list, name = "Mac_res")
  head(object$Macrophage1)
  FeaturePlot(object = object, features = "Macrophage1",pt.size = 1.5,
              max.cutoff = 1.0,min.cutoff = 0,
              cols=c("darkblue","cyan","springgreen","yellow","orange","red"),
              label = FALSE)+
    labs(title="Macrophage CD68,MRC1,C1QA,F13A1")+
    theme(plot.title = element_text(size=10))
  
  marker_gene_list2 <- list(c("VCAN","SYTL3","THBS1","AREG","VDR","SORL1"))
  object <- AddModuleScore(object, features = marker_gene_list2, name = "Mac_VCAN")
  head(object$Mac_VCAN1)
  FeaturePlot(object = object, features = "Mac_VCAN1",pt.size = 1.5,
              max.cutoff = 1.4,min.cutoff =0.4,
              cols=c("darkblue","cyan","springgreen","yellow","orange","red"),
              label = FALSE)+
    labs(title="VCAN,SYTL3,THBS1")+
    theme(plot.title = element_text(size=10))
  
  marker_gene_list3 <- list(c("RUNX2","HLA-DRA","ATF7IP2","CIITA"))
  object <- AddModuleScore(object, features = marker_gene_list3, name = "Mac_HLA")
  head(object$Mac_HLA1)
  FeaturePlot(object = object, features = "Mac_HLA1",pt.size = 1.5,
              max.cutoff =0.8,min.cutoff =0.3,
              cols=c("darkblue","cyan","springgreen","yellow","orange","red"),
              label = FALSE)+ labs(title="HLA-DRA,HLA-DRB1,HLA-DQB1,RUNX2,ATF7IP2,CIITA")+
    theme(plot.title = element_text(size=10))
  
  marker_gene_list4 <- list(c( "GPNMB","SNTB1","TPRG1","PLA2G7"))
  object <- AddModuleScore(object, features = marker_gene_list4, name = "Giant")
  head(object$Giant1)
  FeaturePlot(object = object, features = "Giant1",pt.size = 1.5,
              max.cutoff =1.1,min.cutoff =0.3,
              cols=c("darkblue","cyan","springgreen","yellow","orange","red"),
              label = FALSE)+ labs(title="GPNMB,SNTB1,TPRG1,PLA2G7")+
    theme(plot.title = element_text(size=10))
  
  marker_gene_list5 <- list(c( "FLT3","CD1C","ZBTB46","IL1R2","IL18R1","ADAM19"))
  object <- AddModuleScore(object, features = marker_gene_list5, name = "DCs")
  FeaturePlot(object = object, features = "DCs1",pt.size = 1.5,
              max.cutoff =1.0,min.cutoff =0.3,
              cols=c("darkblue","cyan","springgreen","yellow","orange","red"),
              label = FALSE)+ labs(title="FLT3,CD1C,ZBTB46,IL1R2,IL18R1,ADAM19")+
    theme(plot.title = element_text(size=10))
  
  marker_gene_list6 <- list(c( "FCGR3A","LILRB1","ECE1"))
  object <- AddModuleScore(object, features = marker_gene_list6, name = "Mono")
  FeaturePlot(object = object, features = "Mono1",pt.size = 1.5,
              max.cutoff =1,min.cutoff =0.4,
              cols=c("darkblue","cyan","springgreen","yellow","orange","red"),
              label = FALSE)+ labs(title="FCGR3A,LILRB1,ECE1,NOCT,PRAM1,LILRA5")+
    theme(plot.title = element_text(size=10))

  
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  # Dotplot
  Idents(seuratobj) <- factor(Idents(seuratobj))
  markers.to.plot <- c("F13A1","CD163L1","SCN9A",
                       "VCAN","SYTL3","THBS1",
                       "HLA-DRB1","HLA-DRA","RGS1",
                       "GPNMB","SNTB1","TPRG1",
                       "FLT3","CD1C","ZBTB46",
                       "FCGR3A","LILRB1","ECE1")
  Idents(seuratobj) <- factor(Idents(seuratobj), levels = c("Mac_res","Mac_VCAN","Mac_HLA-DR","Giant cell_GPNMB","DCs_CD1C","Mono_CD16"))
  DotPlot(seuratobj, features = markers.to.plot, dot.scale =5)+ #, split.by = "group"
    #scale_color_gradientn(colours = rainbow(5))+
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 90))+
    scale_colour_gradient2(low = "blue", mid = "white", high = "red")
  
  
  
  #violinplot for sarcoid celltype
  plots <- VlnPlot(seuratobj, features = c("SLC19A3","ADIPOQ","TRHDE","RP11-286B14.1","PTGER3"),
                   pt.size = 0)
  plots
  ggsave(plots, "violinplot sarcoid adipocyte.pdf")
  VlnPlot(seuratobj, features = c("CFAP61",	"MYH6",	"MIR133A1HG",	"MYO18B",	"FHOD3"),pt.size = 0)
  ggsave("violinplot sarcoid cardiomyocyte.png")
  
  VlnPlot(seuratobj, features = c("NRG1",	"PKHD1L1",	"BMX",	"BMP6",	"NFATC1"),pt.size = 0)
  ggsave("violinplot sarcoid endocardium.png")
  
  VlnPlot(seuratobj, features = c("ADGRD1",	"MEG3",	"PDGFD",	"COL4A3",	
                                  "ABI3BP",	"FGF14"),pt.size = 0)
  ggsave("violinplot sarcoid fibroblast.png",height = 20, width = 20)
  
  VlnPlot(seuratobj, features = c("GRIK3",	"CHL1",	"PCSK2",	"XKR4",	"NRXN1"),pt.size = 0)
  ggsave("violinplot sarcoid glia.png",height = 20, width = 20)
  
  VlnPlot(seuratobj, features = c("CPA3",	"AC004791.2",	"HPGDS",	"HDC",	"SLC18A2","KIT"),pt.size = 0)
  ggsave("violinplot sarcoid mast.png",height = 20, width = 20)
  
  VlnPlot(seuratobj, features = c("CD163",	"F13A1","RBM47","MERTK",
                                  "CTSB",	"STAB1",	"FCGR2A",	"FGD2"),pt.size = 0)
  ggsave("violinplot sarcoid myeloid2.png",height = 20, width = 20)
  
  VlnPlot(seuratobj, features = c("PRKCQ",	"THEMIS",	"CD247"),pt.size = 0)
  ggsave("violinplot sarcoid NK&T2.png",height = 20, width = 20)
  
  VlnPlot(seuratobj, features = c("DACH1",	"GUCY1A2",	"RGS5",	"KCNAB1",	"NR2F2-AS1",	"FRMD3"),pt.size = 0)
  ggsave("violinplot sarcoid pericytes.png",height = 20, width = 20)
  
  VlnPlot(seuratobj, features = c("CARMN",	"NTRK3",	"MYH11",	"SLIT3",	
                                  "RCAN2",	"KCNAB1",	"PLA2G5"),pt.size =0.1)
  ggsave("violinplot sarcoid SMC.png",height = 20, width = 20)
  

  # Find marker gene one cellype vs. other one or more celltypes
  Mac_res.markers <- FindMarkers(seuratobj, ident.1 = "Mac_res", min.pct = 0.25)
  head(Mac_res.markers, n = 50)
  write.csv(Mac_res.markers,"Mac_res.markers.csv") 
  
  Giant_cell.markers <- FindMarkers(seuratobj, ident.1 = "Giant cell_GPNMB", min.pct = 0.25)
  head(Giant_cell.markers, n = 50)
  write.csv(Giant_cell.markers,"DEG Giant_cell VS HLA-DR.csv") 
  
  Mac_giant_vs_inf.markers <- FindMarkers(seuratobj, ident.1 = "Giant cell_GPNMB",ident.2 =c("Mac_VCAN","Mac_HLA-DR"), min.pct = 0.25)
  head(Mac_giant_vs_inf.markers, n = 50)
  write.csv(Mac_giant_vs_inf.markers,"Mac_giant_vs_inf.markers.csv") 
  
  Mac_res_vs_Mac_inf.markers <- FindMarkers(seuratobj, ident.1 = "Mac_res",ident.2 =c("Mac_VCAN","Mac_HLA-DR"), min.pct = 0.25)
  head(Mac_res_vs_Mac_inf.markers, n = 50)
  write.csv(Mac_res_vs_Mac_inf.markers,"Mac_res_vs_Mac_inf.markers.csv") 
  
  
  # Find all marker genes
  myeloid.markers0.8 <- FindAllMarkers(seuratobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # Heatmap 
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  myeloid.markers0.8 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  top10 <- myeloid.markers0.8 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  DoHeatmap(seuratobj, features = top10$gene) + NoLegend()
  ggsave("Heatmap105reso.8.png",width = 15,height = 10)
  #using select genes to do heatmap
  for_heatmap <- myeloid.markers%>% 
    filter(gene %in% c("F13A1","CD163L1","SCN9A",
                       "VCAN","SYTL3","THBS1",
                       "HLA-DRB1","HLA-DRA","RGS1",
                       "GPNMB","SNTB1","TPRG1",
                       "FLT3","CD1C","ZBTB46",
                       "FCGR3A","LILRB1","ECE1")) 
  DoHeatmap(seuratobj, features = for_heatmap$gene)
  ggsave("Heatmap selected gene.png")
  
  # DE genes in Mac_res marker genes(res vs inf)
  rm(list = ls())
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  dat <-read.csv("Mac_res.markers(pos &neg).csv")
  colnames(dat)[1]="symbol"
  colnames(dat)[2]="avg_log2FC"
  colnames(dat)[3]="p_val_adj"
  dat<-dat %>%
    arrange(desc(avg_log2FC))
  logFC_t=0.5
  P.Value_t = 0.05
  n1 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC > logFC_t)
  n2 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC < -logFC_t)
  change = ifelse(n1,"up",ifelse(n2,"down","stable"))
  dat <- mutate(dat,change)
  table(dat$change)
  deg_up<-dat[dat$change=="up",]
  write.csv(deg_up,"DE_up Sarcoid res vs inf.csv")
  deg_down<-dat[dat$change=="down",]
  write.csv(deg_down,"DE_down Sarcoid vs lympho.csv")
  
  # DE genes in giant cells vs Mac_HLA-DT
  rm(list = ls())
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  dat <-read.csv("DEG Giant_cell VS HLA-DR.csv")
  colnames(dat)[1]="symbol"
  colnames(dat)[3]="avg_log2FC"
  colnames(dat)[6]="p_val_adj"
  dat<-dat %>%
    arrange(desc(avg_log2FC))
  logFC_t=0.5
  P.Value_t = 0.05
  n1 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC > logFC_t)
  n2 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC < -logFC_t)
  change = ifelse(n1,"up",ifelse(n2,"down","stable"))
  dat <- mutate(dat,change)
  table(dat$change)
  deg_up<-dat[dat$change=="up",]
  write.csv(deg_up,"DE_up Sarcoid res vs inf.csv")
  deg_down<-dat[dat$change=="down",]
  write.csv(deg_down,"DE_down Sarcoid vs lympho.csv")
  
  
  # DE genes in Giant vs inf marker genes
  rm(list = ls())
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  dat <-read.csv("giant cell vs others/Mac_giant_vs_inf.markers.csv")
  colnames(dat)[1]="symbol"
  colnames(dat)[2]="avg_log2FC"
  colnames(dat)[3]="p_val_adj"
  dat<-dat %>%
    arrange(desc(avg_log2FC))
  logFC_t=0.5
  P.Value_t = 0.05
  n1 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC > logFC_t)
  n2 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC < -logFC_t)
  change = ifelse(n1,"up",ifelse(n2,"down","stable"))
  dat <- mutate(dat,change)
  table(dat$change)
  deg_up<-dat[dat$change=="up",]
  write.csv(deg_up,"DE_up Sarcoid Giant_vs_mac_inf.csv")
  deg_down<-dat[dat$change=="down",]
  write.csv(deg_down,"DE_down Sarcoid Giant_vs_mac_inf.csv")

  
  # DE genes in res vs MAC_inf marker genes
  rm(list = ls())
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  dat <-read.csv("Mac_res_vs_Mac_inf.markers.csv")
  colnames(dat)[1]="symbol"
  colnames(dat)[2]="avg_log2FC"
  colnames(dat)[3]="p_val_adj"
  dat<-dat %>%
    arrange(desc(avg_log2FC))
  logFC_t=0.5
  P.Value_t = 0.05
  n1 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC > logFC_t)
  n2 = (dat$p_val_adj < P.Value_t)&(dat$avg_log2FC < -logFC_t)
  change = ifelse(n1,"up",ifelse(n2,"down","stable"))
  dat <- mutate(dat,change)
  table(dat$change)
  deg_up<-dat[dat$change=="up",]
  write.csv(deg_up,"DE_up Mac_res_vs_Mac_inf.csv")
  deg_down<-dat[dat$change=="down",]
  write.csv(deg_down,"DE_down Mac_res_vs_Mac_inf.csv")
  
  #label head(10) genes
  x1 = head(dat[dat$change=="up",] ,15)
  #write.csv(x1,"sarc vs lympho down.genes.csv")
  x2 = tail(dat[dat$change=="down",],15)
  #write.csv(x2,"sarc vs lympho up.genes.csv")
  for_label = rbind(x1,x2)
  
  library(ggplot2)
  p <- ggplot(data = dat, 
              aes(x = avg_log2FC, 
                  y = -log10(p_val_adj))) +
    geom_point(alpha=0.6, size=6, 
               aes(color=change)) +
    ylab("-log10(Pvalue)")+
    scale_color_manual(values=c("red", "grey","blue"))+
    geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
    theme_bw()
  p
  volcano_plot <- p +
    geom_point(size = 6.0, shape = 1, data = for_label) +
    ggrepel::geom_label_repel(
      aes(label = symbol),data = for_label, size=3,color="black",max.overlaps = Inf)
  volcano_plot+xlim(-2.5,2.5)
  
  
  #KEGG analysis by cell types
  rm(list=ls())
  library(Seurat)
  library(gplots)
  library(ggplot2)
  myeloid_marker<-read.csv("myeloid_marker correct.csv") 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  ids=bitr(myeloid_marker$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
  myeloid_marker=merge(myeloid_marker,ids,by.x='gene',by.y='SYMBOL')
  View(myeloid_marker)
  gcSample=split(myeloid_marker$ENTREZID, myeloid_marker$cluster) 
  ## KEGG
  kegg <- compareCluster(gcSample,
                         fun = "enrichKEGG",
                         organism = "hsa", pvalueCutoff = 0.05
  )
  write.csv(kegg,"kegg cluters combined.csv")
  #select kegg for plot in csv
  kegg_selected<-read.csv("kegg cluters combined selected.csv")
  
  kegg@compareClusterResult<-kegg_selected
  
  
  p <- dotplot(kegg,showCategory=50)
  p + theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
  
  
  
  