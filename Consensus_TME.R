source("~/Documents/analysis_script_RNAseq/Universal_Function.R")
devtools::install_github("cansysbio/ConsensusTME")
packages=c("ConsensusTME","tximport","org.Hs.eg.db","clusterProfiler","openxlsx","reshape2","ggplot2","curl")
fun_packages(packages)

Fun_ConsensusTME<-function (dir){
  setwd(dir)
  filelist <- list.files(dir) 
  # get names of dataset
  element=loadData(filelist)
  DataList=element[[2]]
  setwd("~")
  
  #get normalized data here is TPM, could be FPKM as well
  data_norm=normData(DataList,"TPM")
  # the data will be SYMBOL setted
  data_pret=normTransSymbol(data_norm)
  colnames(data_pret)=gsub("Pos",'Neg',colnames(data_pret),fixed=T)
  
  #write.table(data_pret, file="~/normalized_TPM_Symbol_withCOV7.csv", sep=",", quote=F)
  #load directly the data
  #data_pret=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withoutCov7.csv",header=T,row.names=1)
  bulkExpMatrix <- as.matrix(data_pret)
  return (bulkExpMatrix)}

bulkExpMatrix=Fun_ConsensusTME("~/Documents/RSEM_res")
consen=ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "LUAD", statMethod = "singScore")

#write.csv(consen, file="~/TMEConsensus_withoutCov7.csv")
################################################## fin of calculate score of consensusTME ###################################
################## not nessaire, if you don't want to #################################################################
#personlizer our own signature to see the score is significant or not
#Covid, DEGs

rt_degs=read.table("~/Documents/analysis_script_RNAseq/data/DGE_statRes_Deseq2.csv",sep=",",header=T,row.names = 1) 
rt=read.xlsx("~/Documents/COVID19_signature.xlsx") 
#gene interferons
geneList=read.xlsx("~/Documents/Genes_list.xlsx")
row_annot=melt(t(geneList))
colnames(row_annot)=c("GenesNature","id","GeneName")
row_annot=na.omit(row_annot)
#gene TLS
TLS=read.xlsx("~/Documents/signature_TLS.xlsx") 
row_TLS=melt(t(TLS))
colnames(row_TLS)=c("GenesNature","id","GeneName")
row_TLS=na.omit(row_TLS)
gene_design=list()
gene_design[["Covid"]]=rt$Genes
gene_design[["DEGs"]]=rownames(rt_degs)
gene_design[["INF"]]=row_annot$GeneName
gene_design[["VirusSensor"]]=geneList$Virus_sensor
gene_design[["IFN_receptor"]]=geneList$typeI_IFN_receptor
gene_design[["ISG"]]=geneList$ISG
gene_design[["Immunocheck_act"]]=geneList$Immune_checkpoint_activator
gene_design[["Immunocheck_inhibitor"]]=geneList$Immune_checkpoint_inhibitor
gene_design[["cytokines"]]=geneList$Cytokines
gene_design[["TLS"]]=row_TLS$GeneName
gene_design[["TLSCancer"]]=TLS$Signature.TLS.Cancer
gene_design[["TLSMEL"]]=TLS$Signature.TLS.MEL
gene_design[["TLSCRC"]]=TLS$Signature.TLS.CRC

consen=ConsensusTME::geneSetEnrichment(bulkExpMatrix, gene_design)
colnames(consen)=c(rep("COV",6),rep("NEG",16))

res_COV7=consen[,order(colnames(consen),decreasing = T)]

######################################The part in commun that need to be executed for violin plot#############################################
#within COV7
colnames(consen$Scores)=c(rep("COV",6),rep("NEG",16))
res_COV7=consen$Scores[,order(colnames(consen$Scores),decreasing = T)]
df_cov7=data.frame(reshape2::melt(res_COV7))

#arrange the data
colnames(df_cov7)=c("Cell_Type","Sample","value")

########## violin plot ###########
Fun_violin(df_cov7,"COV","NEG")

############# pheatmap ########### it's not nessairy

Info_gene_CT=genes$`Cell population`
names(Info_gene_CT)=genes$`HUGO symbols`
df_gene_CT=data.frame(Info_gene_CT)
df_CT=data_pret[rownames(df_gene_CT),]
CT_gene=merge(df_CT,df_gene_CT,by="row.names")
rownames(CT_gene)=CT_gene$Row.names

df_CT <-CT_gene[order(CT_gene$Info_gene_CT),]
CT_heatmap=df_CT[2:(max(ncol(df_CT))-1)]
Info=df_CT$Info_gene_CT
row_sample_Info=data.frame(CellType=Info)
rownames(row_sample_Info)=rownames(CT_heatmap)

#rename according our samples and annotation withCOV7 or withoutCOV7
sample_info_commun=Fun_annot("withoutCOV7")
rownames(sample_info_commun)=colnames(CT_heatmap)
substring=substring(rownames(sample_info_commun),1,3)

COV_amount=length(grep("COV", substring))
pheatmap(CT_heatmap, scale="row",cluster_cols = F, cluster_rows = F,clustering_method="ward.D2",
         show_rownames = T, border=FALSE,display_numbers = F,annotation_row=row_sample_Info,
         annotation_col = sample_info_commun,
         main ="Heatmap of CT for all genes",gaps_col =c(COV_amount,COV_amount+10),color= colorRampPalette(c("navy", "white", "firebrick3"))(10))

mcpGeneSet <- ConsensusTME::methodSignatures$MCP.Counter
genesetenrichssement=ConsensusTME::geneSetEnrichment(bulkExpMatrix, mcpGeneSet)

