
#### Function for loading the packages
fun_packages <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      #Si le package ne peut pas être chargé, on l'installe
      install.packages( i , dependencies = TRUE )
      #Charger le package après l'installation
      require( i , character.only = TRUE )
    }
    if (!requireNamespace("BiocManager", quietly = TRUE))
    {install.packages("BiocManager")}
  }
}

###### Function for the load matrix count ######
## Function which contains a list of data from RSEM prepared for analysis for DESEQ2 and edgeR :ele[[1]]
# And the pre-prepared characteres for other normalizaed data will use with function normData: ele[[2]]
loadData<-function(sampleList){
  message("Attention !! choose to input your file, if the number of input samples dones't correspond to what you want, please check out your folder and make sure there is just your raw data file")
  data <- readline(prompt = "Do you want sample POS? (Yes or No)")
  df_fl=as.data.frame(sampleList)
  df_withoutbis=df_fl[-8,]
  df_withoutbis=df_withoutbis[-6]
  if (data=="No"){
      print("This is a dataset without POS and COV7 sample")
      df_withoutCov7=df_withoutbis[-7]
      data_selected=df_withoutCov7[1:16]
    }
  
  else if (data=="Yes"){
    data1 <- readline(prompt = "Do you want sample COV7? (Yes or No)")
    if(data1=="No"){
      print("This is a dataset without COV7 samples")
      df_withoutCov7=df_withoutbis[-7]
      data_selected=df_withoutCov7
    }
    else if(data1=="Yes"){
      print("This is a dataset only without bis sample")
      data_selected=df_withoutbis
    }
    else{
      message("Please give an argument between Yes or No")
    }}
  else{
    message("Please input some arguments")
  }
  
  txi.rsem <- tximport(data_selected, type = "rsem", txIn = FALSE, txOut = FALSE)
  colnames(txi.rsem$counts)=substring(data_selected,1,8)
  colnames(txi.rsem$abundance)=substring(data_selected,1,8)
  colnames(txi.rsem$length)=substring(data_selected,1,8)
  return (list(txi.rsem,data_selected))
  }

###### This function allowed to get the normalization data, it could be FPKM as well if needed #####
# change the value of function : new_data[[sampleNames[i]]]=list_data[[sampleNames[i]]]$FPKM
normData<-function(charctere_pret){
  #rearrange the data to get other type of normalization from RSEM
  sampleNames <- substring(charctere_pret,1,8)
  list_data=list()
  new_data=list()
  for (i in 1:length(charctere_pret)){
    list_data[[sampleNames[i]]]=read.table(file.path(dir, charctere_pret[i]), header = TRUE)
    rownames(list_data[[sampleNames[i]]])=list_data[[sampleNames[i]]]$gene_id
    new_data[[sampleNames[i]]]=list_data[[sampleNames[i]]]$TPM
    names(new_data[[sampleNames[i]]])=rownames(list_data[[sampleNames[i]]])
  }
  
  df_all=lapply(1:length(charctere_pret), function(i) {
    new_data[i] %>% 
      data.frame()
  })
  
  df_norm=data.frame(Reduce(cbind, df_all))
  return(df_norm)
  
}

##### This function can change the ENSEMBL ID to SYMBOL ##########
normTransSymbol<-function(ENSEMBLfile){
#transfer the ensembl id to symbol
rank_names <- bitr(rownames(ENSEMBLfile), 
                   fromType="ENSEMBL", 
                   toType=c("SYMBOL"),
                   OrgDb="org.Hs.eg.db")
expression_filter=data.frame(ENSEMBLfile[rank_names$ENSEMBL,])
expression_filter$SYMBOL=rank_names$SYMBOL
merge_data=merge(rank_names,expression_filter,by="SYMBOL",all=F)
rownames(merge_data)=make.names(merge_data$SYMBOL,unique=T)
data_pret=merge_data[ , -which(colnames(merge_data) %in% c("SYMBOL","ENSEMBL"))]

return(data_pret)
}

####### this function is a function of anntotation for heatmap (with or without Cov7)  ########
Fun_annot<-function(sampleNames){
  if (sampleNames=="withCOV7"){
    #within COV7
    histology=c(rep("adénocarcinome",4),"carcinome adénosquameux","carcinome épidermoïde","carcinome neuroendocrine",
                rep("adénocarcinome",4),"X", rep("adénocarcinome",3),"carcinome neuroendocrine",rep("adénocarcinome",3),rep("carcinome épidermoïde",3),"adénocarcinome")
    age=c("40-60",rep("60-80",2),"40-60",rep("60-80",2),"30","40-60",rep("60-80",3),"X",rep("40-60",2),rep("60-80",9))
    sexe=c(rep("M",2),rep("F",2),rep("M",2),"F","M",rep("F",2),"M","X",rep("M",2),"F",rep("M",7),"F")
    Groups=c(rep("COV",7),rep("Neg",16))
    sample_info_commun = data.frame(
      Sample_Type=Groups,
      Histology_type=histology,
      Age=age,
      Sexe=sexe)
  }
  else if(sampleNames=="withoutCOV7"){
    #without COV7
    histology=c(rep("adénocarcinome",4),"carcinome adénosquameux","carcinome épidermoïde",
                rep("adénocarcinome",4),"X", rep("adénocarcinome",3),"carcinome neuroendocrine",rep("adénocarcinome",3),rep("carcinome épidermoïde",3),"adénocarcinome")
    age=c("40-60",rep("60-80",2),"40-60",rep("60-80",2),"40-60",rep("60-80",3),"X",rep("40-60",2),rep("60-80",9))
    sexe=c(rep("M",2),rep("F",2),rep("M",3),rep("F",2),"M","X",rep("M",2),"F",rep("M",7),"F")
    Groups=c(rep("COV",6),rep("Neg",16))
    sample_info_commun = data.frame(
      Sample_Type=Groups,
      Histology_type=histology,
      Age=age,
      Sexe=sexe)
  }
  else{message("Error, please input correct sample")}
  
  return(sample_info_commun)
  
}
############### This function used for the loading of normalized data (part of shiny, normal data analyse doesn't involve) #############
Fun_loadNormFile<-function(dataNorm){
  args <- commandArgs(trailingOnly = FALSE)
  basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
  #load normalization data
    if (dataNorm=="Deseq2"){
      normalized_counts=read.csv(paste0(basedir,"../data/normalized_medianOfRatios_Deseq2.csv"),header=T,row.names=1)
    }
    else if (dataNorm=="edgeR"){
      normalized_counts=read.csv(paste0(basedir,"../data/normalized_TMM_edgeR.csv"),header=T,row.names=1)
    }
    else{
      message("error of input the normalized data")
    } 
  return (normalized_counts)
}

########### This function is applied for prepare and generation of volcanoplot (for analysis) ##############
Fun_Volcano<-function(DGE_sorted){
VC=EnhancedVolcano(DGE_sorted,
                lab = rownames(DGE_sorted),
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0, 8),
                pCutoff = 0.05,
                FCcutoff = 1,
                labSize = 5.0,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1/2,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                colConnectors = 'black',
                title = "Volcano plot",
                legendLabels = c("NS", expression(Log[2] ~ FC), "adjusted p-value", expression(p - adj ~ and
                                                                                               ~ log[2] ~ FC)),
                ylab = bquote(~-Log[10] ~ italic(Padj)),
                axisLabSize = 18) 
return (VC)}

########### This function is applied for prepare and generation of heatmap (for shiny) ##############
Fun_heatmap<-function(data_DEGs,normalized_counts,sample_info_commun){
#add annotation for up and downregulated
  up_deseq=data_DEGs[which(data_DEGs$padj<0.05 & data_DEGs$log2FoldChange>=1),]
  up_deseq$regulated="up"
  down_deseq=data_DEGs[which(data_DEGs$padj<0.05 & data_DEGs$log2FoldChange<=(-1)),]
  down_deseq$regulated="down"
  total_deseq=rbind(up_deseq, down_deseq)
  
  #rename according our samples
  rownames(sample_info_commun)=colnames(normalized_counts) 
  
  #row(gene) annotation
  row_sample_deseq=subset(total_deseq,select=regulated)
  degs_norm_deseq <- normalized_counts[rownames(row_sample_deseq),]
  logTPM <- function(x) {return(log2(x+1))}
  df_norm=degs_norm_deseq %>% mutate_if(is.numeric, logTPM)
  
  scale_data_deseq=degs_norm_deseq%>% pheatmap:::scale_rows()
  annot_degs=HeatmapAnnotation(df=sample_info_commun,boxplot=anno_boxplot(df_norm))
  annot_row_degs_deseq2=rowAnnotation(df=row_sample_deseq,bar=anno_barplot(rowMeans(df_norm)))
  Heatmap(scale_data_deseq,name = "Z-score",cluster_rows=F,show_row_names = FALSE,top_annotation=annot_degs,left_annotation=annot_row_degs_deseq2,
          column_split = sample_info_commun$Sample_Type,column_title = "Heatmap of DEGs genes")
  
}

########################### enrichment analysis ##################################
Fun_enrichmment<-function(data_DEGs){
  
# Run GO enrichment analysis 
go <- enrichGO(rownames(data_DEGs), 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05,
               keyType = 'SYMBOL')

return(go)
}

######################## Tis function is for the cnetplot ######################################
Fun_cnetplot<-function(data_DEGs,go){
  cat=dim(go)[1]
  genelist=data_DEGs$log2FoldChange
  names(genelist)=rownames(data_DEGs)
  edox=pairwise_termsim(go, showCategory = cat)
  options(ggrepel.max.overlaps = Inf)
  p1_go <-cnetplot(go,categorySize="padj", foldChange=genelist,showCategory = cat,  node_label="category", 
                   cex_label_gene = 0.8,cex_label_category = 0.8)
  
  p2_go <-cnetplot(edox,categorySize="padj", foldChange=genelist,showCategory = cat, node_label="gene",
                   circular = TRUE, ggrepel.max.overlaps = Inf,cex_label_gene = 0.8,cex_label_category = 0.8)
  
  whole_plot<-cowplot::plot_grid(p1_go, p2_go, labels=LETTERS[1:2])
  return(whole_plot)
}
###################### This function is for the preparation of genelist for all analysis based on GESA ###############################
Fun_genelist<-function(dataNames){
  args <- commandArgs(trailingOnly = FALSE)
  basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
  if (dataNames=="Deseq2"){
    results.DESeq2=read.csv(paste0(basedir,"../data/StatRes_Deseq2.csv"),header=T,row.names=1)
  }
  else if (dataNames=="edgeR"){
    results.DESeq2=read.csv(paste0(basedir,"../data/StatRes_EdgeR.csv"),header=T,row.names=1)
  }
  else{
    message("error of input the normalized data")
  } 
  
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
  return(list(genelist_sort,results.DESeq2))
}
###################### GESA analysis of GO ###############################
Fun_gseGO<-function(genelist_sort){
  
  kk<-gseGO(geneList = genelist_sort,
            OrgDb = "org.Hs.eg.db",
            pvalueCutoff = 0.05,
            eps=0,
            nPermSimple = 10000,
            verbose=F)
  return(kk)
}
###################### GESA analysis of KEGG ###############################
Fun_gesa<-function(geneList_gesa_sort){
  args <- commandArgs(trailingOnly = FALSE)
  basedir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
  kegmt<-read.gmt(paste0(basedir,"../data/c7.all.v7.5.1.symbols.gmt"))
  res_gsea = GSEA(geneList_gesa_sort, TERM2GENE=kegmt,eps=0)
  return(res_gsea)
}

######################This function is for MCPcounter TPM normalisation data, get TPKM normalization, 
#see the function normData()  and the start of mcp_counter script #################

Fun_mcpcounter<-function(df_normalized){
  logTPM <- function(x) {return(log2(x+1))}
  df_norm=df_normalized %>% mutate_if(is.numeric, logTPM)
  genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
  res=MCPcounter.estimate(df_norm,featuresType=c("HUGO_symbols")[1],
                               probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                               genes=genes)
  colnames(res)=substring(colnames(df_norm),1,3)
  df_melt=data.frame(melt(res))
  colnames(df_melt)=c("Cell_Type","Sample","value")
  #violin plot
  
  return(list(df_melt,df_norm,genes))
}


############################ This function ois for TMEConsensus, same principal with MCPcounter ###########
Fun_TME<-function(df_normalized){
  bulkExpMatrix <- as.matrix(df_normalized)
  consen=ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "LUAD", statMethod = "singScore")
  #within COV7
  colnames(consen$Scores)=substring(colnames(df_normalized),1,3)
  df_cov7=data.frame(reshape2::melt(consen$Scores))
  colnames(df_cov7)=c("Cell_Type","Sample","value")
  return(df_cov7)
}

   ###########violin plot for immunology analysis ############
Fun_violin<-function(df_cov7){
  split_cellType=df_cov7 %>% group_split(Cell_Type)
  
  a=length(unique(df_cov7$Cell_Type))
  p <- list()
  for(i in 1:a){
    p[[i]]=ggplot(split_cellType[[i]],aes(x=Sample,y=value,fill=Sample))+geom_violin(trim=F)+scale_fill_brewer(palette = "Accent")+
      geom_signif(comparisons = list(c("COV","Neg")),test = "wilcox.test",map_signif_level = TRUE,y_position=c(max(split_cellType[[i]]$value)-0.2,max(split_cellType[[i]]$value),max(split_cellType[[i]]$value)+0.3))+
      theme(legend.position = "none",axis.title.x = element_blank())+
      geom_boxplot(width=0.1,fill="white")+ggtitle(split_cellType[[i]]$Cell_Type)
  }
  p_all <- grid.arrange(grobs=p,ncol=5)#,top="Violin plots of Cell Types without COV7"
  return(p_all)
}

   ############ heatmap for immunology analysis  ############

Fun_heatmcp<-function(df_norm,genes,sample_info_commun){
  Info_gene_CT=genes$`Cell population`
  names(Info_gene_CT)=genes$`HUGO symbols`
  df_gene_CT=data.frame(Info_gene_CT)
  df_CT=df_norm[rownames(df_gene_CT),]
  CT_gene=merge(df_CT,df_gene_CT,by="row.names")
  
  #rownames
  rownames(CT_gene)=CT_gene$Row.names
  df_CT <-CT_gene[order(CT_gene$Info_gene_CT),]
  CT_heatmap=df_CT[2:(max(ncol(df_CT))-1)]
  Info=df_CT$Info_gene_CT
  row_sample_Info=data.frame(CellType=Info)
  rownames(row_sample_Info)=rownames(CT_heatmap)
  
  #rename according our samples
  rownames(sample_info_commun)=colnames(CT_heatmap)
  substring=substring(rownames(sample_info_commun),1,3)
  
  scale_data_CT=CT_heatmap%>% pheatmap:::scale_rows()
  annot_row_CT=rowAnnotation(CellType=Info,Means = anno_barplot(rowMeans(CT_heatmap)))
  annot_col=HeatmapAnnotation(df=sample_info_commun,boxplot=anno_boxplot(CT_heatmap))
  Heatmap(scale_data_CT,name = "Z-score",cluster_rows = F,
          left_annotation=annot_row_CT,top_annotation=annot_col,row_title = NULL,
          column_split = substring,row_split = Info)#,column_title = "Heatmap of Cell type for samples without COV7"
  
}


