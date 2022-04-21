source("~/Documents/analysis_script_RNAseq/script/Universal_Function.R")

packages=c("org.Hs.eg.db","clusterProfiler","dplyr","ggplot2","ggrepel","GSVA","limma","tidyr",
           "enrichplot","ggplot2","fgsea","msigdbr","ggnewscale","openxlsx","tximport","DESeq2")
fun_packages(packages)
load("~/data_GOenrichissement.RData")

########################## Attention!!! all the method by using the GESA method the data will be the result of Deseq2/edgeR analysis####
######## Other normal result like GO and Kegg, use just the data of DEGs ######################

#deseq2
args <- commandArgs(trailingOnly = FALSE)
basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))

rt=read.table(paste0(basedir,"data/DGE_statRes_Deseq2.csv"),sep=",",header=T,row.names = 1) 

#edgeR
rt=read.table("~/Documents/analysis_script_RNAseq/data/DGE_statRes_edgeR.csv",sep=",",header=T,row.names = 1) 

#regulated

input_value <- readline(prompt = "Enter the regulation that you want to applied (up,down or all): ")
up_or_down <- function(input_value) {
  if (input_value=="all") {
    df=rt
  }
  else
  {df=rt[rt$regulated==input_value,]}
  return(df)}

df=up_or_down(input_value) 

############################### GO #####################################
#OVA 
go <- enrichGO(rownames(df), 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05,
               keyType = 'SYMBOL')
cat=dim(go)[1]

#barplot
barplot(go,showCategory=cat,drop=T,title=paste("GO enrichment of pathway of DEGs",input_value,"regulated",sep=" "))

#geneList prparation
genelist=df$log2FoldChange
names(genelist)=rownames(df)
edox=pairwise_termsim(go, showCategory = cat)
options(ggrepel.max.overlaps = Inf)
#net work presentation
p1_go <- cnetplot(edox, node_label="category", 
                    cex_label_category = 1.2) 
p2_go <- cnetplot(edox, node_label="gene", 
                    cex_label_gene = 0.8) 
cowplot::plot_grid(p1_go, p2_go, labels=LETTERS[1:2])
#heatplot
heatplot(edox, foldChange=genelist, showCategory=cat)
#treeplot
treeplot(edox, hclust_method = "average")
#barplot in the terms of ontology
barplot(go,showCategory=cat, split="ONTOLOGY",title="GO enrichment of edgeR Outliers split with ontology")+ facet_grid(ONTOLOGY~.,scale="free")
ggplot(go, showCategory=cat,aes(Description, Count, fill=ONTOLOGY), split='ONTOLOGY') + geom_col()  + 
  theme(axis.text.x=element_text(angle=-90, hjust=0)) + facet_grid(.~ONTOLOGY, scale="free_x")
#dotplot
dotplot(go,showCategory=cat,split="ONTOLOGY",title=paste("GO enrichment of pathway of DEGs",input_value,"regulated",sep=" "))+
  facet_grid(ONTOLOGY~.,scale="free")+scale_y_discrete(labels=function(x) str_wrap(x, width=80))

##########################     KEGG     ##############################
load("~/Documents//DEGs_multivariates_DesEdg_withoutCov7.RData")

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

#KEGG pathway gene set enrichment analysis using GESAs
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

#########################   GSEA  ##################################

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

cat_gesa=dim(res_gsea@result)[1]
dotplot(res_gsea,split=".sign")+facet_wrap(~.sign,scales="free")
p1_gesa <- cnetplot(res_gsea, node_label="category", 
               cex_label_category = 1.2) 
p2_gesa <- cnetplot(res_gsea, node_label="gene", 
               cex_label_gene = 0.8) 
cowplot::plot_grid(p1_gesa, p2_gesa, labels=LETTERS[1:2])
gseaplot2(res_gsea,4500,color="red",pvalue_table = T)

edo_gesa=pairwise_termsim(res_gsea, showCategory = cat_gesa)
heatplot(edo_gesa, foldChange=geneList_gesa_sort, showCategory=cat_gesa)
treeplot(edo_gesa, hclust_method = "average")

#######################   GSVA   #######################
####using the normalization data or rawdata
setwd("~/Documents/RSEM_res")
dir="~/Documents/RSEM_res"
filelist <- list.files(dir)
#load all data 
ele1=loadData(filelist)
txi.counts=ele1[[1]]
TPMnorm=normData(txi.counts$counts)
gene_symbol=normTransSymbol(TPMnorm)
res_gsva=gsva(data.matrix(gene_symbol),kegmt,min.sz>1,method="gsva",verbose=T,kcdf="Poisson")
