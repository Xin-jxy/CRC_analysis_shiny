source("~/Documents/analysis_script_RNAseq/script/Universal_Function.R")

#Liste des packages necessaires:
packages=c("plyr", "dplyr", "tm", "readxl", "wordcloud", "SnowballC", "stringdist", "tidytext",
                "rmarkdown", "knitr", "quanteda", "reshape", "stringr", "RecordLinkage", "plotly",
                "data.table", "rvest",  "shiny", "shinydashboard", "shinyWidgets", "DT","shinythemes","InteractiveComplexHeatmap",
                "ggplot2" ,"clusterProfiler","pheatmap" ,"edgeR" , "statmod","DESeq2","NMF","ggbeeswarm",
                "genefilter","pheatmap","ade4","viridis","tidyverse","dplyr","tximeta","tximport","openxlsx",
                'EnhancedVolcano',"ggrepel","biomaRt","reshape","RColorBrewer","VennDiagram","ComplexHeatmap",
                "crosstalk","ggvenn","DT","curl","ConsensusTME","org.Hs.eg.db","GSVA","limma","tidyr",
                "enrichplot","fgsea","msigdbr","ggnewscale","devtools","MCPcounter","gridExtra","ggsignif") 
#Lancer la fonction de chargement/installation: 
fun_packages(packages)
#s'il y a des nécessité des packages avec bioconductor 
#BiocManager::install(c("tximeta","DESeq2","tximport","clusterProfiler","edgeR",'EnhancedVolcano',"biomaRt" ))

#load the RSEM rawdata
setwd("~/Documents/RSEM_res")
dir="~/Documents/RSEM_res"
filelist <- list.files(dir)
txi.rsem <- tximport(filelist, type = "rsem", txIn = FALSE, txOut = FALSE)
colnames(txi.rsem$counts)=substring(filelist,1,5)
colnames(txi.rsem$abundance)=substring(filelist,1,5)
colnames(txi.rsem$length)=substring(filelist,1,5)
#function is in Universal_Function.R
ele=loadData(filelist)
txi.rsem=ele[[1]]

#verify our data
head(txi.rsem$counts)

#write.table(txi.rsem$counts, file="~/matrixCount.csv", sep=",", quote=F)

##############################  Set up our model  ################################

colnames(txi.rsem$counts)=gsub("Pos",'Neg',colnames(txi.rsem$counts),fixed=T)
Groups=substring(colnames(txi.rsem$counts),1,3)

#the histology group without COV7
histology_type=c(rep("adénocarcinome",4),rep("Other",2),rep("adénocarcinome",4),"Other", 
            rep("adénocarcinome",3),"Other",rep("adénocarcinome",3),rep("Other",3),"adénocarcinome")

#the histology group within COV7
#histology_type=c(rep("adenocarcinome",4),rep("Other",3),rep("adenocarcinome",4),"Other", 
                 #rep("adenocarcinome",3),"Other",rep("adenocarcinome",3),rep("Other",3),"adenocarcinome")



# In our data set, we want to identify the genes differentially between the different groups
#take account of only one principale factors
#sampleTable <- data.frame(condition = factor(Groups))

#take account of multiple factors
sampleTable <- data.frame(condition = factor(Groups), Histological_type=factor(histology_type))
rownames(sampleTable) <- colnames(txi.rsem$counts)

# TODO : chose a reference group (first-level factor)
# Here any group can do, if there was a control or time-0 group, you would chose this one.
#the reference group is less "severe" or important one
reference_group <- "Neg"
sampleTable$condition<- relevel(sampleTable$condition, reference_group)
sampleTable$condition

reference_histo<-"Other"
sampleTable$Histological_type<- relevel(sampleTable$Histological_type, reference_histo)
sampleTable$Histological_type

######################################################################
#                             Deseq2                                 #
######################################################################

#get DESeqDataSetFromTximport functionnal
txi.rsem$abundance=txi.rsem$abundance[apply(txi.rsem$length,1,
                function(row) all(row !=0 )),]

txi.rsem$counts <-txi.rsem$counts[apply(txi.rsem$length,1,
                             function(row) all(row !=0 )),]

txi.rsem$length <-
  txi.rsem$length[apply(txi.rsem$length,1,
                             function(row) all(row !=0 )),]


#sampleTable$sample=rownames(sampleTable)

# Replicates must have the same group identifier, so that the package can use them to estimate dispersion
# generate the DESeqDataSet : here we test for 2 samples only with the design = ~group argument 

#Only one main factor
#DESeq.ds <- DESeqDataSetFromTximport(txi.rsem,sampleTable , ~ condition)

#two or more factors
DESeq.ds <- DESeqDataSetFromTximport(txi.rsem,sampleTable , ~ condition+Histological_type)

# investigate  library  sizes
colSums(counts(DESeq.ds))

################################# Visualization of non-supervised analysis #############################

#### rlog transformations for exploratory analysis and visualization
# PCA and clustering should be done on normalized and preferably transformed read counts, so
# that the high variability of read counts does not occlude potentially informative data trends

# Here I used the regularized log-transformation of the DESeq2 package that reduce variance 
# of low read counts by using the dispestion-mean trend seen for the entire data set as a reference
rld <- rlog(DESeq.ds, blind = FALSE)
rlog.norm.counts <- assay(rld)

# Read count correlations
# we compute 1 minus the correlation of read counts to obtain a distance
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))

# we use the distance to perform a hierarchical clustering and plot the resulting tree
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
     main = "rlog  transformed  read  counts\ndistance: Pearson  correlation")
# TODO : 
# -check that replicates are clustered in the tree
# -see if some branches are unexpectedly long 
###### probl?mes : deux r?plicats sont ?loign?s ？
# PCA : the PCA will inform on the structure of the signal. It helps to verify that the overall experiment produced a signal.
plotPCA(rld,"condition")
pcaData=plotPCA(rld,"condition",returnData=T)
pcaData_sort <- pcaData[order(pcaData$condition,decreasing=F),]

#PCA without Cov7 samples
ggplot(data = pcaData, aes(x = pcaData[,1], y = pcaData[,2],label=rownames(pcaData))) +scale_color_manual(values=c("#FFC20A","#0C7BDC"))+
  theme_minimal()+geom_point(aes(color = condition))+geom_text_repel(hjust=0, vjust=0,size=3)+labs(title="PCA with all patient samples",x="PC1:23% variance",y="PC2:18% variance")+
  theme(plot.title = element_text(size = 15,hjust=0.5),axis.text.x = element_text( hjust = -0.25,size=15),
        axis.text.y=element_text(size=15),axis.title.x= element_text(size=15),axis.title.y= element_text(size=15),
        legend.key.size = unit(2,"line"),legend.title = element_text(size = 15),legend.text = element_text(size = 15))

#all samples
ggplot(data = pcaData, aes(x = pcaData[,1], y = pcaData[,2])) +scale_color_manual(values=c("#FFC20A","#0C7BDC"))+
  theme_minimal()+geom_point(aes(color = condition,size=20))+labs(title="PCA of all samples",x="PC1:23% variance",y="PC2:19% variance")+
  theme(plot.title = element_text(size = 15,hjust=0.5),axis.text.x = element_text( hjust = -0.25,size=15),
        axis.text.y=element_text(size=15),axis.title.x= element_text(size=15),axis.title.y= element_text(size=15),
        legend.key.size = unit(2,"line"),legend.title = element_text(size = 15),legend.text = element_text(size = 15))

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(rld))))
pheatmap(sampleDistMatrix,main = "Heatmap of correlation between each samples")

#-----------------------------------------------------------------------
#============   Differential Expression with DESeq2     ================
#-----------------------------------------------------------------------

# Warning : The DGE analysis has to be performed on the raw read counts (untransformed, not normalized
# for sequencing depth)

#============== Differential Expression Gene Analysis ==================

# The DESeq function will performed three successive steps of analysis :
# sequencing depth normalization between the samples
# gene-wise dispersion estimation across all samples
# fits a negative  binomial  GLM  and applies  Wald  statistics  to each  gene

DESeq_new <- DESeq(DESeq.ds)
plotDispEsts(DESeq_new)
resultsNames(DESeq_new)

#Don't cutoff the outlier
#DESeq_new_outliers <- DESeq(DESeq.ds,minReplicatesForReplace = Inf )

#shrink when there is many dispersion of expression of gene, donesn't have impact for DGE
#DESeq_shrink <- lfcShrink(DESeq_new, coef="condition_COV_vs_NEG",type="apeglm")

############if need, we will performe the group to compare, if not, the first level group is reference#########
###########!!!!!!control group must in the last position!!!!!################
# Building the results table 
# TODO : chose two groups to compare. The group name must match one present in the group vector defined line 93 in SampleInfo
#group1 <- "COV"
#group2 <- "NEG"
# here the method of adjustement of p-values is Benjamini and Hochberg (BH), also called FDR (false discovery rate)
#DGE.results <- results(DESeq_new, c("condition", treated, untreated), pAdjustMethod = "BH")
#DGE.results

#Don't cut off the outlier data
#DGE_outliers_res=results(DESeq_new_outliers,pAdjustMethod = "BH",cooksCutoff = FALSE)
#DGE.results=DGE_outliers_res

#Only one factor is present 
#DGE.results=results(DESeq_new,pAdjustMethod = "BH")

#calculate with a condition that we want(main factors)
DGE.results=results(DESeq_new,pAdjustMethod = "BH",contrast=c("condition","COV","Neg"))
summary(DGE.results)

#filter the data that we don't want (p-value is NA, p-value ajuste is NA, basemean is 0)
filter_DGE_res=DGE.results[which(! is.na(DGE.results$pvalue)&
                            ! is.na(DGE.results$padj)&
                            DGE.results$baseMean>0),]

mcols(filter_DGE_res, use.names = TRUE)

results.DESeq2 <- filter_DGE_res

# Histogram of p-values. Mainly useful to see if there is a low statistical power.
hist(DGE.results$pvalue, col = "grey", border = "white", xlab = "p-value of all", ylab = "frequency",
     main = "frequencies  of p-values")
# NB: This plot is best formed by excluding genes with very small counts, 
# which otherwise generate spikes in the histogram.
hist(DGE.results$pvalue[DGE.results$baseMean > 1], col = "grey", border = "white", xlab = "p-value of basemean>1", ylab = "frequency",
     main = "frequencies  of p-values")

# MA plot:  provides a global view of the relationship between the expression change between conditions
# (log-ratios M), the average espression strength of the genes(average mean, A) and the abilitiy of 
# the algorithm to detect differential gene expression
# TODO : chose an appropriate alpha risk. Do the plot with different alpharisk to see the differences in statistical power.
plotMA(results.DESeq2 , alpha = 0.0001,  main = "MA plot: group1 vs group2 (padj < 0.0001)", ylim = c(-4,4))
abline(h = c(-1,1), col = "blue", lty = 2)
mtext(c(paste("-2 fold"), paste("+ 2 fold")), side = 4, at = c(-1, 1),
      cex = 0.8, line = 0.5, col = "blue")


#============== Differential Expression Extraction ==================

# Read counts of the top gene and plot, mainly useful that there was no mistake when setting the constrasts (which groups to compare)
topGene <- rownames(results.DESeq2)[which.min(results.DESeq2$padj)]
plotCounts(DESeq.ds, gene = topGene, intgroup="condition")

#sort with the p-value adjust
DGE.results.sorted  <- results.DESeq2[order(results.DESeq2$padj), ]
#transform ENSEMBL genenames to SYMBOL,function is at Universal_Function.R
DGE_sorted=normTransSymbol(DGE.results.sorted)
#we think the p-value ajuste <0.05 is signifiance
DGE_subset=subset(DGE_sorted,padj < 0.05)
#DGE_subset$ENSEMBL<-genenames_DGE$ENSEMBL

# identify  genes  with  the  desired  adjusted p-value cut -off
DGEgenes=DGE_sorted[DGE_sorted$padj<0.05,]
head(DGEgenes,10)

#up and downregulated
up_deseq=DGE_sorted[which(DGE_sorted$padj<0.05 & DGE_sorted$log2FoldChange>=1),]
up_deseq$regulated="up"
down_deseq=DGE_sorted[which(DGE_sorted$padj<0.05 & DGE_sorted$log2FoldChange<=(-1)),]
down_deseq$regulated="down"
total_deseq=rbind(up_deseq, down_deseq)
total_deseq_sort=total_deseq[order(total_deseq$padj), ]

#write.table(total_deseq, file="~/DGE_statRes_Deseq2.csv", sep=",",quote=F,row.names = T,col.names=NA)

#-----------------------------------------------------------------------
#============      Visualization of data of Deseq2      ================
#-----------------------------------------------------------------------

################ the table of DEGs ################
#COV_vs_NEG, so log2FC>0-->COV

grp=c("UP","DOWN")
grp.col <- c("#568875", "#73FAFC")

DESeq_result_df <- DGE_sorted %>% data.frame() %>% arrange(log2FoldChange, padj)

log2fc_threshold <- with(DESeq_result_df, mean(abs(log2FoldChange)) + 1.5 * sd(abs(log2FoldChange)))
message(paste("Threshold of log2Foldchange [Mean+1.5(SD)] is", log2fc_threshold))
pval <- 0.05
DESeq_result_df[which(DESeq_result_df$log2FoldChange >= log2fc_threshold & DESeq_result_df$padj < pval), "Enrichment"] <- grp[2]
DESeq_result_df[which(DESeq_result_df$log2FoldChange <= -log2fc_threshold & DESeq_result_df$padj < pval), "Enrichment"] <- grp[1]
DESeq_result_df[which(abs(DESeq_result_df$log2FoldChange) < log2fc_threshold | DESeq_result_df$padj >= pval), "Enrichment"] <- "Nonsignf"
table(DESeq_result_df$Enrichment)

############volcano plot############
#function in Universal_Function.R
DGE_sorted[which(total_deseq$regulated=="up"),"Enrichment"]<-"Up"
DGE_sorted[which(total_deseq$regulated=="down"),"Enrichment"]<-"Down"
DGE_sorted[is.na(DGE_sorted$Enrichment),"Enrichment"]<-"Nonsignf"
table(DGE_sorted$Enrichment)
keyvals.colour <- ifelse(
  DGE_sorted$padj<0.05 & DGE_sorted$log2FoldChange>=1, '#E66100',
  ifelse(DGE_sorted$padj<0.05 & DGE_sorted$log2FoldChange<=(-1), '#5D3A9B',
         'grey'))
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == '#5D3A9B'] <- 'Down'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Nonsignf'
names(keyvals.colour)[keyvals.colour == '#E66100'] <- 'Up'

volcano=EnhancedVolcano(DGE_sorted,
                lab = rownames(DGE_sorted),
                selectLab = rownames(DGE_sorted)[which(names(keyvals.colour) %in% c('Up', 'Down'))],
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals.colour,
                #ylim = c(0, 5.5),
                pCutoff = 0.05,
                FCcutoff = 1,
                labSize = 5.0,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1/2,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 4.0,
                colConnectors = 'black',
                xlab = bquote(~~Log[2]~ 'fold change'),
                maxoverlapsConnectors=Inf,
                drawConnectors = F,
                widthConnectors = 0.75,
                title = 'Volcano plot of DEGs colored with COV and Neg patients',
                legendLabels = c("NS", expression(Log[2] ~ FC), "adjusted p-value", expression(p - adj ~ and
                                                                                               ~ log[2] ~ FC)),
                ylab = bquote(~-Log[10] ~ italic(Padj)),
                axisLabSize = 18) 
library("ggallin")
volcano+scale_y_continuous(trans=pseudolog10_trans,limits=c(1,5.5))

############### Heatmap #############

#load("~/Documents/Data_Deseq2_edgeR_multivariate_withoutPOS_COV7.RData")
#normalization of Deseq2 for the visualization of heatmap
DESeq.ds <- estimateSizeFactors(DESeq.ds)
normalized_counts <- counts(DESeq.ds, normalized=TRUE)

#write.table(normalized_counts, file="~/normalized_medianOfRatios_withoutCov7.csv", sep=",", quote=F,col.names = NA)

#rld_heatmap <- vst(DESeq_new, blind=FALSE)
#de<- rownames(DGEgenes)
#de_mat <- assay(rld_heatmap)[de,]

#estabilish the annotation col or row
### prepare pheatmap
histology=c(rep("adenocarcinoma",4),"adenosquamous carcinoma","squamous cell carcinoma",
            rep("adenocarcinoma",4),"squamous cell carcinoma", rep("adenocarcinoma",3),"neuroendocrine carcinoma",rep("adenocarcinoma",3),rep("squamous cell carcinoma",3),"adenocarcinoma")
age=c("40-60",rep("60-80",2),"40-60",rep("60-80",2),"40-60",rep("60-80",4),rep("40-60",2),rep("60-80",9))
sexe=c(rep("M",2),rep("F",2),rep("M",3),rep("F",2),rep("M",4),"F",rep("M",7),"F")
sample_info_commun = data.frame(
  Sample_Type=Groups,
  Histology_type=histology,
  Age=age,
  Sex=sexe)

#For DEGs
#annot_colnames, from ENSEMBL to SYMBOL
norm_data=normTransSymbol(normalized_counts)
#rename according our samples
rownames(sample_info_commun)=colnames(norm_data) 

#row(gene) annotation
row_sample_deseq=subset(total_deseq,select=regulated)
degs_norm_deseq <- norm_data[rownames(row_sample_deseq),]
degs_norm_deseq=na.omit(degs_norm_deseq)

des="Deseq2"
pheatmap(degs_norm_deseq, scale="row",cluster_cols = F, cluster_rows = T,clustering_method="mcquitty",
         show_rownames = T, annotation_row = row_sample_deseq,clustering_distance_rows="correlation",
         border=FALSE,display_numbers = F,annotation_col = sample_info_commun,
         main = paste("Heatmap of", des,"for DEGs genes",sep=" "),gaps_col =5,color= colorRampPalette(c("navy", "white", "firebrick3"))(10))


######################  ComplexHeatmap   ####################

logTPM <- function(x) {return(log2(x+1))}
df_norm=degs_norm_deseq %>% mutate_if(is.numeric, logTPM)

colours_col <- list("Histology_type"=c("adenocarcinoma"="#6699CC",
                         "adenosquamous carcinoma"="#40B0A6",
                         "squamous cell carcinoma"="#D35FB7",
                         "neuroendocrine carcinoma"="#FEFE62"),
                "Age"=c("40-60"="#99DDFF","60-80"="#CCDDAA"),
                "Sex"=c("M"="#EE8866","F"="#AA4499"),
                "Sample_Type"= c("COV" = "#0C7BDC" ,"Neg" ="#FFC20A"))
colours_row <-list("regulated"=c("up"="#E66100","down"="#5D3A9B"))

row_sample_deseq=total_deseq[rownames(degs_norm_deseq),]
row_sample_deseq=subset(row_sample_deseq,select=regulated)

scale_data_deseq=degs_norm_deseq%>% pheatmap:::scale_rows()
annot_degs=HeatmapAnnotation(df=sample_info_commun,which="col",col=colours_col,boxplot=anno_boxplot(df_norm))
annot_row_degs_deseq2=rowAnnotation(df=row_sample_deseq,col=colours_row,bar=anno_barplot(rowMeans(df_norm)))
Heatmap(scale_data_deseq,name = "Z-score",cluster_rows=F,show_row_names = FALSE,show_column_names = FALSE,top_annotation=annot_degs,left_annotation=annot_row_degs_deseq2,
        column_title = "Heatmap of Deseq2 for DEGs analysis", #column_split = sample_info_commun$Sample_Type,
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))

#######################################################################
#                              edgeR                                 #
######################################################################

#============== Data Normalisation for library sample size ==================s
cts <- txi.rsem$counts
normMat <- txi.rsem$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

#filter NA data wih dataframe
df_normCounts=as.data.frame(normCts)
df_normCounts=rapply(df_normCounts, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
#df_normCounts=df_normCounts[apply(df_normCounts,1,function(row) all(row !=0 )),]

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(df_normCounts) * colSums(df_normCounts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

df_normMat=data.frame(normMat)
df_no_nan=df_normMat[is.finite(rowSums(df_normMat)),]
matrix_normM=as.matrix(df_no_nan)

cols <- intersect(rownames(df_normCounts), rownames(df_no_nan))
df_normCounts_fn<- df_normCounts %>% filter(row.names(df_normCounts) %in% cols)
normCts_final=as.matrix(df_normCounts_fn)

#normCts_output=normTransSymbol(normCts_final)
#write.table(normCts_output, file="~/normalized_TMM_edgeR.csv", sep=",", quote=F,row.names = T,col.names = NA)

y <- DGEList(normCts_final)
list_calculate <- scaleOffset(y, matrix_normM)

#construction of model
Des_Group <- factor(paste(sampleTable$condition,sampleTable$Histological_type,sep="."))
desigN<-model.matrix(~0+Des_Group,data=sampleTable)
cbind(sampleTable,Group=Des_Group)
colnames(desigN) <- levels(Des_Group)

keep <- filterByExpr(list_calculate,desigN)
list_pret <- list_calculate[keep, ]

#estimate the dispersion of the data
#estimeDisp <- estimateDisp(list_pret, design=desigN,robust=T)
estimeDisp <- estimateGLMRobustDisp(list_pret, design=desigN)
plotBCV(estimeDisp)
plotMDS(estimeDisp)

#fit <- glmQLFit(estimeDisp, contrast=desigN,robust=T)
fit <- glmFit(estimeDisp, desigN,robust=T)

colnames(fit)
plotQLDisp(fit)

#construction of contrasts.matrix
contrasts.matrix<-makeContrasts(CovVSNegOther=COV.Other-Neg.Other,
                                CovVSNegAde=COV.adenocarcinome-Neg.adenocarcinome,levels=desigN)

#test qlf
#qlf<-glmLRT(fit,contrast=contrasts.matrix[,"CovVSNegOther"])

qlf<-glmLRT(fit,contrast=contrasts.matrix[,"CovVSNegAde"])
#glmQLFTest
qlf$table$FDR <- p.adjust(qlf$table$PValue,"BH")
qlf_res=topTags(qlf,n = nrow( qlf$table ))$table
deg_edger_qlf=qlf_res[qlf_res$FDR<0.05,]

FDR_sorted=qlf$table[order(qlf$table$FDR), ]
FDR_sorted[1:10,]
#decision test to know the significance gene or not
is.de <- decideTests(qlf,adjust.method="BH", p.adjust=0.05)
summary(is.de)

#get symbole
idnames_sorted_edger=normTransSymbol(FDR_sorted)
class(idnames_sorted_edger)
idnames_sorted_edger= idnames_sorted_edger %>% rename("log2FoldChange"="logFC" ,"padj"="FDR" )

#up and downregulated
up=idnames_sorted_edger[which(idnames_sorted_edger$padj<0.05 & idnames_sorted_edger$log2FoldChange>=1),]
up$regulated="up"
down=idnames_sorted_edger[which(idnames_sorted_edger$padj<0.05 & idnames_sorted_edger$log2FoldChange<=(-1)),]
down$regulated="down"
total_edger=rbind(up, down)
total_edger_sort=total_edger[order(total_edger$padj), ]

write.table(total_edger, file="~/DGE_statRes_edgeR.csv", sep=",", quote=F,row.names = T,col.names=NA)

#-----------------------------------------------------------------------
#============      Visualization of data of EdgeR       ================
#-----------------------------------------------------------------------

#COV_vs_NEG, so log2FC>0-->COV
##########volcano plot#########
Fun_Volcano(idnames_sorted_edger)

keyvals.colour <- ifelse(
  idnames_sorted_edger$padj<0.05 & idnames_sorted_edger$log2FoldChange>=1, '#E66100',
  ifelse(idnames_sorted_edger$padj<0.05 & idnames_sorted_edger$log2FoldChange<=(-1), '#5D3A9B',
         'grey'))
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == '#5D3A9B'] <- 'Down'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Nonsignf'
names(keyvals.colour)[keyvals.colour == '#E66100'] <- 'Up'

volcano=EnhancedVolcano(idnames_sorted_edger,
                        lab = rownames(idnames_sorted_edger),
                        selectLab = rownames(idnames_sorted_edger)[which(names(keyvals.colour) %in% c('Up', 'Down'))],
                        x = 'log2FoldChange',
                        y = 'padj',
                        colCustom = keyvals.colour,
                        #ylim = c(0, 5.5),
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        labSize = 5.0,
                        labCol = 'black',
                        labFace = 'bold',
                        colAlpha = 1/2,
                        legendPosition = 'right',
                        legendLabSize = 15,
                        legendIconSize = 4.0,
                        colConnectors = 'black',
                        xlab = bquote(~~Log[2]~ 'fold change'),
                        maxoverlapsConnectors=Inf,
                        drawConnectors = F,
                        widthConnectors = 0.75,
                        title = 'Volcano plot of DEGs colored with COV and Neg patients',
                        legendLabels = c("NS", expression(Log[2] ~ FC), "adjusted p-value", expression(p - adj ~ and
                                                                                                       ~ log[2] ~ FC)),
                        ylab = bquote(~-Log[10] ~ italic(Padj)),
                        axisLabSize = 18) 
library("ggallin")
volcano+scale_y_continuous(trans=pseudolog10_trans,limits=c(0,4.5))

##########################heatmap###############################

### prepare pheatmap

rownames(sample_info_commun)=colnames(list_pret$counts)
up_down_edge=c(rep("up",length(rownames(data.frame(up)))),rep("down",length(rownames(data.frame(down)))))

#For DEGs
norm_edger=normTransSymbol(normCts_final)
row_sample_edger=subset(total_edger,select=regulated)
degs_norm_edger <- norm_edger[rownames(row_sample_edger),]


#for histo
#colnames(degs_norm_edger)=histology_cov
#degs_norm_edger=degs_norm_edger[,order(colnames(degs_norm_edger))]

edg="EdgeR"
pheatmap(degs_norm_edger,scale="row", cluster_cols = F, cluster_rows = T,clustering_method="mcquitty",
         show_rownames = T, annotation_row = row_sample,clustering_distance_rows="correlation",
         border=FALSE,display_numbers = F,main = paste("Heatmap of",edg,"for DEGs genes",sep=" "),
         gaps_col =6,color= colorRampPalette(c("navy", "white", "firebrick3"))(10))
#annotation_col = sample_info_commun,

######################  ComplexHeatmap   ####################

logTPM <- function(x) {return(log2(x+1))}
df_norm=degs_norm_edger %>% mutate_if(is.numeric, logTPM)

scale_data_edger=degs_norm_edger%>% pheatmap:::scale_rows()
annot_degs_edger=HeatmapAnnotation(df=sample_info_commun,col=colours_col,boxplot=anno_boxplot(df_norm))
annot_row_degs_edger=rowAnnotation(df=row_sample_edger,col=colours_row,bar = anno_barplot(rowMeans(df_norm)))
Heatmap(scale_data_edger,name = "Z-score",cluster_rows=F,show_row_names = FALSE,show_column_names = FALSE,top_annotation=annot_degs_edger,left_annotation=annot_row_degs_edger,
        column_title = "Heatmap of EdgeR for DEGs analysis",column_title_gp = gpar(fontsize = 20, fontface = "bold"))
#column_split = sample_info_commun$Sample_Type,

######################################################################
#              comparaison of two toolkits                           #
######################################################################

load("~/Documents/DEGs_multivariates_DesEdg_withoutCov7.RData")
#test if the result of two packages is correlated
combined_res_des=merge(data.frame(DGE_sorted),data.frame(idnames_sorted_edger),by="row.names",all=F)
summary(combined_res_des)

intersect_res=intersect(rownames(total_edger_sort),rownames(total_deseq_sort))
intersect(rownames(total_edger_sort[1:10,]),rownames(total_deseq_sort[1:10,]))
#write.table(intersect_res, file="~/res_intersect.txt", quote=F,row.names = F)
edger_eliminate=total_edger_sort[-which(rownames(total_edger_sort)%in% intersect_res),]
deseq_eliminate=total_deseq_sort[-which(rownames(total_deseq_sort)%in% intersect_res),]
combined_res_des <- combined_res_des %>% 
  mutate(Groups =
           case_when(Row.names %in% intersect_res~"Intersect",
                     Row.names %in% rownames(edger_eliminate)~"edgeR",
                     Row.names %in% rownames(deseq_eliminate)~"DESeq2"
                     ))
combined_res_des$Groups<-combined_res_des$Groups %>% replace_na("NoSigf")
ggplot(data = combined_res_des, mapping = aes(x = log2FoldChange.x, y = log2FoldChange.y,color=Groups,alpha = Groups)) +scale_colour_brewer(palette = "Set2")+ geom_point()+
  scale_alpha_discrete(range = c(1,0.1,1,1))+theme_minimal()+
  labs(title="Correlation between analysis DESeq2 and edgeR",x="Deseq2LFC",y="EdgeRLFC")

x = list(
  EdgeR = rownames(total_edger_sort),
  DESeq2 = rownames(total_deseq_sort))

overlaps <- calculate.overlap(x)

names(x) <- c("EdgeR","Deseq2")
ggvenn(x, 
       fill_color = c("#FC8D62","#66C2A5"),
       stroke_size = 0.5, set_name_size = 4)
brewer.pal(5, "Set2")
#bscols(ggplotly(ggvenn(x, 
 # fill_color = c("#0073C2FF", "#EFC000FF"),
  #stroke_size = 0.5, set_name_size = 4)),datatable(data.frame(overlaps$a3)))

