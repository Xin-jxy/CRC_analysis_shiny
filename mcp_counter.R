source("~/Documents/analysis_script_RNAseq/Universal_Function.R")

packages=c("devtools","MCPcounter","org.Hs.eg.db","clusterProfiler","dplyr","reshape","ggplot2","gridExtra","ggsignif","pheatmap")
fun_packages(packages)

#devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")

###################################### MCPcounter use the data normalized by log(TPM+1) ##########################

################ Load the data normalized  ################

Fun_MCPcounter<-function (dir){
  setwd(dir)
  #load orignaly data
  filelist <- list.files(dir) 
  
  # get names of dataset
  ele=loadData(filelist)
  DataList=ele[[2]]
  
  #get normalized data here is TPM, could be FPKM as well, check Universal_Function.R fille out
  data_norm=normData(DataList,"TPM")
  # the data will be SYMBOL setted
  data_pret=normTransSymbol(data_norm)
  #change colnames pos to neg
  colnames(data_pret)=gsub("Pos",'Neg',colnames(data_pret),fixed=T)
  #write.table(data_pret, file="~/Normalized_TPM_withinCov7.csv", sep=",",quote=F,row.names = T,col.names=NA)
  #logTPM
  logTPM <- function(x) {return(log2(x+1))}
  df_norm=data_pret %>% mutate_if(is.numeric, logTPM)
  setwd("~")
  
  #within COV7
  res_COV7=MCPcounter.estimate(df_norm,featuresType=c("HUGO_symbols")[1],
                      probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                      genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
  )
  return (res_COV7)
}

res_COV7=Fun_MCPcounter("~/Documents/RSEM_res")
#write.csv(res_COV7, file="~/MCPCounter_withoutCov7.csv")

#rearrange the result
colnames(res_COV7)=substring(colnames(df_norm),1,3)
res_COV7=res_COV7[,order(colnames(res_COV7),decreasing = T)]
df_cov7=data.frame(melt(res_COV7))

#arrange the data
colnames(df_cov7)=c("Cell_Type","Sample","value")

########## violin plot ###########
Fun_violin(df_cov7,"COV","Neg")

