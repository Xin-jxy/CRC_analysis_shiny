source("~/Documents/analysis_script_RNAseq/Universal_Function.R")

packages=c("devtools","MCPcounter","org.Hs.eg.db","clusterProfiler","dplyr","reshape","ggplot2","gridExtra","ggsignif","pheatmap")
fun_packages(packages)

devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")

###################################### MCPcounter use the data normalized by log(TPM+1) ##########################

################ Load the data normalized  ################

#load orignaly data
setwd("/home/xin/Documents/RSEM_res")
dir="/home/xin/Documents/RSEM_res"
filelist <- list.files(dir) 

# get names of dataset
ele=loadData(filelist)
DataList=ele[[2]]

#get normalized data here is TPM, could be FPKM as well, check Universal_Function.R fille out
data_norm=normData(DataList)
# the data will be SYMBOL setted
data_pret=normTransSymbol(data_norm)
#change colnames pos to neg
colnames(data_pret)=gsub("Pos",'Neg',colnames(data_pret),fixed=T)
#write.table(data_pret, file="~/Normalized_TPM_withinCov7.csv", sep=",",quote=F,row.names = T,col.names=NA)
#logTPM
logTPM <- function(x) {return(log2(x+1))}
df_norm=data_pret %>% mutate_if(is.numeric, logTPM)


#within COV7
res_COV7=MCPcounter.estimate(df_norm,featuresType=c("HUGO_symbols")[1],
                    probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                    genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
)


#write.csv(res_COV7, file="~/MCPCounter_withinCov7.csv")

#rearrange the result
colnames(res_COV7)=substring(colnames(df_norm),1,3)
res_COV7=res_COV7[,order(colnames(res_COV7),decreasing = T)]
df_cov7=data.frame(melt(res_COV7))

#arrange the data
colnames(df_cov7)=c("Cell_Type","Sample","value")
split_cellType=df_cov7 %>% group_split(Cell_Type)

########## violin plot ###########
a=length(unique(df_cov7$Cell_Type))
p <- list()

for(i in 1:a){
  p[[i]]=ggplot(split_cellType[[i]],aes(x=Sample,y=value,fill=Sample))+geom_violin(trim=F)+scale_fill_brewer(palette = "Accent")+
    geom_signif(comparisons = list(c("COV","Neg")),test = "wilcox.test",map_signif_level = TRUE,y_position=c(max(split_cellType[[i]]$value)-0.2,max(split_cellType[[i]]$value),max(split_cellType[[i]]$value)+0.3))+
    theme(legend.position = "none",axis.title.x = element_blank())+
    geom_boxplot(width=0.1,fill="white")+ggtitle(split_cellType[[i]]$Cell_Type)
}
p_all <- grid.arrange(grobs=p,ncol=5,top="Violin plots of Cell Types")


#preparation of data
genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
Info_gene_CT=genes$`Cell population`
names(Info_gene_CT)=genes$`HUGO symbols`
df_gene_CT=data.frame(Info_gene_CT)
df_CT=df_norm[rownames(df_gene_CT),]
CT_gene=merge(df_CT,df_gene_CT,by="row.names")
rownames(CT_gene)=CT_gene$Row.names


############# pheatmap ###########

df_CT <-CT_gene[order(CT_gene$Info_gene_CT),]
CT_heatmap=df_CT[2:(max(ncol(df_CT))-1)]
Info=df_CT$Info_gene_CT
row_sample_Info=data.frame(CellType=Info)
rownames(row_sample_Info)=rownames(CT_heatmap)
#within COV7
histology=c(rep("adénocarcinome",4),"carcinome adénosquameux","carcinome épidermoïde","carcinome neuroendocrine",
            rep("adénocarcinome",4),"X", rep("adénocarcinome",3),"carcinome neuroendocrine",rep("adénocarcinome",3),rep("carcinome épidermoïde",3),"adénocarcinome")
age=c("40-60",rep("60-80",2),"40-60",rep("60-80",2),"30","40-60",rep("60-80",3),"X",rep("40-60",2),rep("60-80",9))
sexe=c(rep("M",2),rep("F",2),rep("M",2),"F","M",rep("F",2),"M","X",rep("M",2),"F",rep("M",7),"F")
sample_info_commun = data.frame(
  Histology_type=histology,
  Age=age,
  Sexe=sexe)

#without COV7
histology=c(rep("adénocarcinome",4),"carcinome adénosquameux","carcinome épidermoïde",
            rep("adénocarcinome",4),"X", rep("adénocarcinome",3),"carcinome neuroendocrine",rep("adénocarcinome",3),rep("carcinome épidermoïde",3),"adénocarcinome")
age=c("40-60",rep("60-80",2),"40-60",rep("60-80",2),"40-60",rep("60-80",3),"X",rep("40-60",2),rep("60-80",9))
sexe=c(rep("M",2),rep("F",2),rep("M",3),rep("F",2),"M","X",rep("M",2),"F",rep("M",7),"F")
sample_info_commun = data.frame(
  Histology_type=histology,
  Age=age,
  Sexe=sexe)

#rename according our samples
rownames(sample_info_commun)=colnames(CT_heatmap)
substring=substring(rownames(sample_info_commun),1,3)

COV_amount=length(grep("COV", substring))
pheatmap(CT_heatmap, scale="row",cluster_cols = F, cluster_rows = F,clustering_method="ward.D2",
         show_rownames = T, border=FALSE,display_numbers = F,annotation_row=row_sample_Info,
         annotation_col = sample_info_commun,
         main ="Heatmap of CT for all genes",gaps_col =c(COV_amount,COV_amount+10),color= colorRampPalette(c("navy", "white", "firebrick3"))(10))


scale_data_CT=CT_heatmap%>% pheatmap:::scale_rows()
annot_row_CT=rowAnnotation(CellType=Info,Means = anno_barplot(rowMeans(CT_heatmap)))
annot_col=HeatmapAnnotation(df=sample_info_commun,boxplot=anno_boxplot(CT_heatmap))
Heatmap(scale_data_CT,name = "Z-score",cluster_rows = F,
        left_annotation=annot_row_CT,top_annotation=annot_col,row_title = NULL,
        column_split = substring,row_split = Info,
        column_title = "Heatmap of Cell type for samples without COV7")


