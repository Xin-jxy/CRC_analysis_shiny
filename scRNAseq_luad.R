library(data.table)
library(Seurat)
library(ggplot2)
library(dplyr)

#load the data that we use after the selection of lung tumor, how to select-->readme
counts <- data.table::fread("/data/data_cremer_lab/scRNAseq/data/GSE131907_LungTumor_log2TPM.txt.gz")
df=data.frame(counts)
list_gene=data.table::fread("/data/data_cremer_lab/scRNAseq/data/list_gene.txt")
rownames(df)<-list_gene$Index

# get DEGs
rt <- read.csv("/data/data_cremer_lab/scRNAseq/data/list_gene_function.csv",header=T,sep=",",row.names = 1)
rownames(rt)=gsub("\\.", "-", rownames(rt))
up=rownames(rt[which(rt$regulated=="up"),])
down=rownames(rt[which(rt$regulated=="down"),])
up_sc=na.omit(df[up,])
down_sc=na.omit(df[down,])

##### combine the up and down genes with total by calculating the mean of up/down ####
df=rbind(df,colSums(up_sc)/length(up_sc))
df=rbind(df,colSums(down_sc)/length(down_sc))
row.names(df)[29635]<-"up"
row.names(df)[29636]<-"down"

#annotation
a <- read.table("/data/data_cremer_lab/scRNAseq/data/GSE131907_Lung_Cancer_cell_annotation.txt",header=T,sep="\t",row.names = 1)
rownames(a)=gsub("_", "-", rownames(a))
colnames(df)=gsub("_", "-", colnames(df))
rownames(df)=gsub("_", "-", rownames(df))
Obj_luad=CreateSeuratObject(counts = df, project = "luad-tumor",meta.data=a)
all.genes <- rownames(Obj_luad)
#construct our object
pbmc <- ScaleData(Obj_luad, features = all.genes)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc) )
Umap <- RunUMAP(pbmc, dims = 1:10)

#give annotation, the annotation could be Cell_subtype or Cell_type
Idents(Umap)=Obj_luad@meta.data$Cell_subtype
DimPlot(Umap, reduction = "umap",label=TRUE)

#the different plots with all up and down genes
FeaturePlot(Umap, features = rownames(rt)[1:227])
RidgePlot(Umap, features = rownames(rt)[1:227], ncol = 3)
VlnPlot(Umap, features = rownames(rt)[1:227])
DotPlot(Umap, features = rownames(rt)[1:227]) + RotatedAxis()

FeaturePlot(Umap, features = c("up","down"),cols=c("#cccccc","#ff6600","#660000"))
FeaturePlot(Umap, features = c("up","down"))+ scale_colour_gradientn(colours = c("#cccccc","#ff6600","#660000"))
RidgePlot(Umap, features = c("up","down"), ncol = 3)
VlnPlot(Umap, features = c("up","down"))
DotPlot(Umap, features = c("up","down") )+ RotatedAxis()

#get genes in clusters of umap
levels(Umap)
MALTB=subset(x = Umap, idents = c("MALT B cells"))
gene_related=MALTB@assays$RNA@var.features

#specific one type of cell-sous type
Markers_maltB=FindMarkers(Umap, ident.1 ="MALT B cells", min.pct = 0.25)
#write.table(Markers_maltB[1:50,],"../scRNAseq/Markergenes_top50_by_maltB.csv",sep=",")

#get the data of all markers in terms of the global comparaison 
MaltB.markers_allcluster <- FindAllMarkers(Umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Top 50 by each cluster
clustering_soustype=MaltB.markers_allcluster %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

#write.table(clustering_soustype,"../scRNAseq/Markergenes_top50_by_clustering_cellsoustype.csv",sep=",")
