source("~/Documents/analysis_script_RNAseq/script/Universal_Function.R")

packages=c("org.Hs.eg.db","clusterProfiler","dplyr","ggplot2","ggrepel","GSVA","limma","tidyr",
           "enrichplot","ggplot2","fgsea","msigdbr","ggnewscale","openxlsx","tximport","DESeq2","GOplot")
fun_packages(packages)
load("~/Documents/data_GOenrichissement.RData")

########################## Attention!!! all the method by using the GESA method the data will be the result of Deseq2/edgeR analysis####
######## Other normal result like GO and Kegg, use just the data of significant DEGs ######################

#deseq2
args <- commandArgs(trailingOnly = FALSE)
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))

rt=read.table(paste0(basedir,"data/DGE_statRes_Deseq2.csv"),sep=",",header=T,row.names = 1) 
rownames(rt)=gsub("\\.", "-", rownames(rt))
################## Get gene function ##################

file_function=Fun_annot_genefunc(rt)
write.table(file_function[,6:8],"~/list_gene_function.csv",sep=",")

#edgeR
#rt=read.table(paste0(basedir,"data/DGE_statRes_edgeR.csv"),sep=",",header=T,row.names = 1) 

############## choose the regulated gene if you want to limit the type of genes ###################
input_value <- readline(prompt = "Enter the regulation that you want to applied (up,down or all): ")
up_or_down <- function(input_value) {
  if (input_value=="all") {
    df=rt
  }
  else
  {df=rt[rt$regulated==input_value,]}
  return(df)}

df=up_or_down(input_value) 

################### WGCNA hubgene ########################
rt=read.table("~/hubgene_KME_yellow_significant.csv",sep=",",header=T)
rownames(rt$gene)=gsub("\\.", "-", rownames(rt$gene))
df=rt$gene
go <- enrichGO(df, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05,
               keyType = 'SYMBOL')
barplot(go,title="GO enrichment of pathway in WGCNA")
go@result=go@result[which(go@result$ONTOLOGY=="BP"),]
dotplot(go, showCategory=20,split="ONTOLOGY",title="Enrichissement d'analyse de GO de process biologique pathway en WGCNA ")


kegg_names<- bitr(df, 
                  fromType="SYMBOL", 
                  toType="ENTREZID",
                  OrgDb="org.Hs.eg.db")

kegg <- enrichKEGG(kegg_names$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH')
kegg_cat=dim(kegg)[1]
barplot(kegg,showCategory=kegg_cat,title="KEGG enrichment of pathway")


############################### GO #####################################
#OVA 
go <- enrichGO(rownames(rt), 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05,
               keyType = 'SYMBOL')
cat=dim(go)[1]

#barplot
barplot(go,showCategory=cat,drop=T,title=paste("GO enrichment of pathway of DEGs",sep=" "))

#geneList prparation
genelist=rt$log2FoldChange
names(genelist)=rownames(rt)
edox=pairwise_termsim(go, showCategory = cat)
options(ggrepel.max.overlaps = Inf)

######################### choose the pathway interested for showing with plot, this step can be skipped ##################
y <- as.data.frame(edox)
pathway_select=y[-c(5,6,11,25,27,28,42,53,54,58,60,63,65,66),]
data_ordered <- pathway_select %>% 
  group_by(ONTOLOGY) %>% # group by condition
  arrange(ONTOLOGY, GeneRatio) %>% # arrange by Condition and GeneRatio
  mutate(GOs_ordered = factor (row_number(), labels = ONTOLOGY)) 
go_down@result = go_down@result[go_down@result$ID %in% data_ordered$ID, ]
###########################################################################################################################

#net work presentation
p1_go <- cnetplot(edox, node_label="category", showCategory=cat,categorySize="pvalue",
                    cex_label_category = 1.2) 
p2_go <- cnetplot(edox, node_label="gene", showCategory=cat,categorySize="pvalue",
                    cex_label_gene = 0.8) 
cowplot::plot_grid(p2_go,p1_go,  labels=LETTERS[1:2])
cnetplot(edox, foldChange=genelist, circular = TRUE, showCategory=cat,categorySize="p.adjust",colorEdge = T) 
#heatplot
heatplot(edox, foldChange=genelist, showCategory=cat)

#treeplot, clustering method could be changed
treeplot(edox, hclust_method = "average")

#barplot in the terms of ontology
barplot(go,showCategory=dim(go)[1], split="ONTOLOGY",title="GO enrichment of edgeR Outliers split with ontology")+ facet_grid(ONTOLOGY~.,scale="free")
#show with ggplot, can be skipped
#ggplot(go, showCategory=cat,aes(Description, Count, fill=ONTOLOGY), split='ONTOLOGY') + geom_col()  + 
  #theme(axis.text.x=element_text(angle=-90, hjust=0)) + facet_grid(.~ONTOLOGY, scale="free_x")

#dotplot separated with MP, BP, CC
dotplot(go,showCategory=dim(go)[1],split="ONTOLOGY",title=paste("GO enrichment of pathway of DEGs up regulated",sep=" "))+
  facet_grid(ONTOLOGY~.,scale="free_y")+scale_y_discrete(labels=function(x) str_wrap(x, width=80))+geom_count()+ scale_size_area(max_size = 15)+
  theme(axis.text.x=element_text(size=rel(2)),axis.text.y=element_text(size=rel(2)),plot.title = element_text(size = 20),legend.title = element_text(size=15),legend.text = element_text(size=15))


###################################################### chose only interested pathway ###########################################
go=go_down
cat=dim(go)[1]
edox=pairwise_termsim(go, showCategory = cat)
y <- as.data.frame(edox)
pathway_select=y[c(1,2,3,4,5,8,9,10,12,13,14,16,17,18,19,21,22,23,24),]
data_ordered <- pathway_select %>% 
  group_by(ONTOLOGY) %>% # group by condition
  arrange(ONTOLOGY, GeneRatio) %>% # arrange by Condition and GeneRatio
  mutate(GOs_ordered = factor (row_number(), labels = ONTOLOGY)) 
go@result = go@result[go@result$ID %in% data_ordered$ID, ]
edox=pairwise_termsim(go, showCategory = dim(go)[1])

dotplot(go_down,showCategory=dim(go_down)[1],split="ONTOLOGY",title=paste("GO enrichment of pathway of DEGs down regulated",sep=" "))+
  facet_grid(ONTOLOGY~.,scale="free_y")+scale_y_discrete(labels=function(x) str_wrap(x, width=80))+geom_count()+ scale_size_area(max_size = 15)+
  theme(axis.text.x=element_text(size=rel(2)),axis.text.y=element_text(size=rel(2)),plot.title = element_text(size = 20),legend.title = element_text(size=15),legend.text = element_text(size=15))
########################################################################################################################################

##### circle plot to show the up/down genes in each pathway
test_go=data.frame(go)
go_select=subset(test_go,select=c("ONTOLOGY","ID","Description","p.adjust","geneID"))
rownames(go_select)=c(1:length(rownames(go_select)))
go_select$ID=c(1:length(rownames(go_select)))
colnames(go_select)=c("Category","ID","Term","adj_pval","Genes")
go_select$Genes <- gsub('/', ',', go_select$Genes)
###################################### selection of pathway interested ############################################
ids=c("1","2","3","4","5","8","9","10","12","13","14","16","17","18","19","21","22","23","24")
go_select=go_select[ids,]
######################################################################################################

#load the file of data (not just list of gene but logFC..., the rename of all col is important)
rt=read.table(paste0(basedir,"data/DGE_statRes_Deseq2.csv"),sep=",",header=T,row.names = 1)
rt$ID=rownames(rt)
rownames(rt)=c(1:length(rownames(rt)))
rt=rename(rt,c("logFC"="log2FoldChange","AveExpr"="baseMean","P.Value"="pvalue","adj.P.Val"="padj"))
##objet of circle plot
circ <- circle_dat(go_select, rt)
GOCircle(circ,nsub = length(rownames(go_select)),lfc.col = c('#E66100','#5D3A9B'))
GOBubble(circ,title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  

#only bp of ontology
go@result=go@result[which(go@result$ONTOLOGY=="BP"),]
dotplot(go,showCategory=dim(go)[1],title=paste("GO enrichment of pathway of Biological Process of DEGs down regulated",sep=" "))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))+geom_count()+ scale_size_area(max_size = 15)+
  theme(axis.text.x=element_text(size=rel(2)),axis.text.y=element_text(size=rel(2)),plot.title = element_text(size = 20),legend.title = element_text(size=15),legend.text = element_text(size=15))

################################################################################################################
##method gesa with GO database,Attention!!! we need to load the database of result in Deseq2 (not significant DEGs)
###############################################################################################################
load("~/Documents//DEGs_multivariates_DesEdg_withoutCov7.RData")
results.DESeq2=read.csv("~/Documents/analysis_script_RNAseq/data/StatRes_Deseq2.csv",header=T,row.names=1)
eg_gsea<- bitr(rownames(results.DESeq2), 
               fromType="ENSEMBL", 
               toType="ENTREZID",
               OrgDb="org.Hs.eg.db")
gene_all=data.frame(results.DESeq2[eg_gsea$ENSEMBL,])
gene_all$ENSEMBL=eg_gsea$ENSEMBL
merge_gene=merge(eg_gsea,gene_all,by="ENSEMBL")
geneList <- merge_gene$log2FoldChange
names(geneList) <- merge_gene$ENTREZID
genelist_sort=sort(geneList,decreasing = T)
kk<-gseGO(geneList = genelist_sort,
          OrgDb = "org.Hs.eg.db",
          pvalueCutoff = 0.05,
          eps=0,
          nPermSimple = 10000,
          verbose=F)
kkcat=dim(kk)[1]
dotplot(kk,showCategory=kkcat,split=".sign")+facet_wrap(~.sign,scales="free")

gseaplot(kk,by="all",title=kk$Description[1],geneSetID=1)
gseaplot(kk,geneSetID=1:3)

a=length(unique(kk$Description))
p <- list()
a=10
for(i in 1:a){
  p[[i]]=gseaplot(kk,by="all",title=kk$Description[i],geneSetID=i)
}
p_all <- grid.arrange(grobs=p,ncol=5,top="GSEA")
ridgeplot(kk,showCategory=kkcat)

##########################     KEGG     ##############################

# KEGG pathway over-representation analysis
kegg_names<- bitr(rownames(rt), 
               fromType="SYMBOL", 
               toType="ENTREZID",
               OrgDb="org.Hs.eg.db")

kegg <- enrichKEGG(kegg_names$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH')
kegg_cat=dim(kegg)[1]
barplot(kegg,showCategory=kegg_cat,title="KEGG enrichment of pathway")

over_kegg=pairwise_termsim(kegg, showCategory = kegg_cat)
edoOver_kegg <- setReadable(over_kegg, 'org.Hs.eg.db', 'ENTREZID')
p1_over <- cnetplot(edoOver_kegg, node_label="category", 
               cex_label_category = 1.2) 
p2_over <- cnetplot(edoOver_kegg, node_label="gene", 
               cex_label_gene = 0.8) 
cowplot::plot_grid(p1_over, p2_over, labels=LETTERS[1:2])

######################################### KEGG database of gene set enrichment analysis using GESA method #######################
results.DESeq2=read.csv("~/Documents/analysis_script_RNAseq/data/StatRes_Deseq2.csv",header=T,row.names=1)
eg_gsea<- bitr(rownames(results.DESeq2), 
               fromType="ENSEMBL", 
               toType="ENTREZID",
               OrgDb="org.Hs.eg.db")
gene_all=data.frame(results.DESeq2[eg_gsea$ENSEMBL,])
gene_all$ENSEMBL=eg_gsea$ENSEMBL
merge_gene=merge(eg_gsea,gene_all,by="ENSEMBL")
#merge_gene=merge_gene %>% rename(log2FoldChange = logFC,padj = FDR)
geneList <- merge_gene$log2FoldChange
names(geneList) <- merge_gene$ENTREZID
genelist_sort=sort(geneList,decreasing = T)


kk<-gseKEGG(geneList = genelist_sort,
            organism="hsa",
            nPerm=1000,
            minGSSize = 120,
            pvalueCutoff = 0.05,
            verbose=F)
kkcat=dim(kk)[1]
dotplot(kk,showCategory=kkcat,title=paste("KEGG enrichment of pathway of DEGs",sep=" "))


edo_kegg=pairwise_termsim(kk, showCategory = kegg_cat)
edox_kegg <- setReadable(edo_kegg, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox_kegg, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox_kegg, node_label="gene", 
               cex_label_gene = 0.8) 
cowplot::plot_grid(p1, p2, labels=LETTERS[1:2])
emapplot(edox)

heatplot(edox_kegg, foldChange=genelist_sort, showCategory=kegg_cat)
browseKEGG(kegg, "4489")
treeplot(edox_kegg, hclust_method = "average")

#############################################   GSEA method with C7 database ##################################
args <- commandArgs(trailingOnly = FALSE)
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
results.DESeq2=read.table(paste0(basedir,"data/StatRes_Deseq2.csv"),sep=",",header=T,row.names = 1) 
df_symbol=normTransSymbol(results.DESeq2)
df_symbol=df_symbol %>% 
  rename(
    log2FoldChange = logFC,
    padj = FDR)
geneList_gesa <- df_symbol$log2FoldChange
names(geneList_gesa) <- rownames(df_symbol)
head(geneList_gesa)
geneList_gesa_sort=sort(geneList_gesa,decreasing = T)
#load pathway
kegmt<-read.gmt("~/Documents/analysis_script_RNAseq/data/c7.all.v7.5.1.symbols.gmt")
res_gsea = GSEA(geneList_gesa_sort, TERM2GENE=kegmt,eps=0)

ridgeplot(res_gsea)

cat_gesa=dim(res_gsea@result)[1]
dotplot(res_gsea,showCategory=cat_gesa,split=".sign")+facet_wrap(~.sign,scales="free")
p1_gesa <- cnetplot(res_gsea, node_label="category", 
               cex_label_category = 1.2) 
p2_gesa <- cnetplot(res_gsea, node_label="gene", 
               cex_label_gene = 0.8) 
cowplot::plot_grid(p1_gesa, p2_gesa, labels=LETTERS[1:2])
gseaplot2(res_gsea,4500,color="red",pvalue_table = T)

edo_gesa=pairwise_termsim(res_gsea, showCategory = cat_gesa)
heatplot(edo_gesa, foldChange=geneList_gesa_sort, showCategory=cat_gesa)
treeplot(edo_gesa, hclust_method = "average")
