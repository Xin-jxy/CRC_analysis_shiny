source("~/Documents/analysis_script_RNAseq/Universal_Function.R")
packages=c("openxlsx","reshape2","ComplexHeatmap",'devtools',"survminer","tidyverse","survival","devtools",
           "cgdsr","DT","ggpubr")

fun_packages(packages)


#######################################################################################
############################# analyse de survie ##########################################
#######################################################################################
mycgds <- CGDS("http://www.cbioportal.org/")

# Test the CGDS endpoint URL using a few simple API tests
test(mycgds) 
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study  
mycancerstudy = getCancerStudies(mycgds)[163,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]


##clinical data##
myclinicaldata <- getClinicalData(mycgds,mycaselist) #230

dat <- myclinicaldata %>% filter(OS_MONTHS > 0)
table(dat$OS_STATUS)
#density of living and death of this dataset luad
ggplot(dat,       
       aes(x = OS_MONTHS, group = OS_STATUS,colour = OS_STATUS,           
           fill = OS_STATUS
       )) + geom_density(alpha = 0.5)

#### with sex
dat$STATUS= ifelse(dat$OS_STATUS == '1:DECEASED', 1, 0)
my.surv <- Surv(dat$OS_MONTHS,dat$STATUS) 
kmfit1 <- survfit(my.surv~dat$SEX) 
ggsurvplot(kmfit1,data = dat)
#with age
my.surv <- Surv(dat$OS_MONTHS,dat$OS_STATUS=='DECEASED') 
kmfit2 <- survfit(my.surv~dat$AGE) 
ggsurvplot(kmfit2,data = dat)

##get RNA_seq analysis##
rt <- read.table("~/Documents/analysis_script_RNAseq/data/DGE_statRes_Deseq2.csv",header=T,sep=",",row.names = 1)
up=rownames(rt[which(rt$regulated=="up"),])
down=rownames(rt[which(rt$regulated=="down"),])

mycaselist <- getCaseLists(mycgds,mycancerstudy)[20,1] 
getGeneticProfiles(mycgds,mycancerstudy) #see RNA_seq data
mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy)[12,1]

################ chose the gene (Up or Down) #####################
expr <- getProfileData(mycgds,down,
                       mygeneticprofile,mycaselist)
exprfilter=data.frame(t(na.omit(t(expr))))

#getProfileData to get the data analysis

data=merge(myclinicaldata[,c('OS_STATUS','OS_MONTHS')],exprfilter[rownames(myclinicaldata),],by="row.names",all.x=TRUE )
data=data[data$OS_MONTHS > 0,]
data=data[!is.na(data$OS_STATUS),]
data$OS_STATUS=as.character(data$OS_STATUS)
my_gene_list=colnames(data)[-(1:4)]

df_select=data[,-2]


############# test de survie ##############
#observer one specific gene, for example here is AMPD1 and TTC16
dat$ampd1_group=ifelse(df_select$AMPD1 > median(df_select$AMPD1),'high','low')
dat$TTC16_group=ifelse(df_select$TTC16 > median(df_select$TTC16),'high','low')
attach(dat)
table(ampd1_group,TTC16_group)


my.surv <- Surv(dat$OS_MONTHS,dat$STATUS)
kmfit1 <- survfit(my.surv~TTC16_group,data = dat) 

ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE)

################## cox #######################
myclinicaldata$STATUS= ifelse(myclinicaldata$OS_STATUS == '1:DECEASED', 1, 0)
data=merge(myclinicaldata[,c('STATUS','OS_MONTHS')],exprfilter[rownames(myclinicaldata),],by="row.names",all.x=TRUE )
gene_list <- gsub(colnames(data), pattern = '-', replacement = '_')
data=data[data$OS_MONTHS > 0,]
data=data[!is.na(data$STATUS),]
formula <- formula(paste('Surv(OS_MONTHS, STATUS)~', gene_list[-(1:4)]))
summary(coxph(formula, data = data))


uni_cox <- function(single_gene){
  formula <- as.formula(paste0('Surv(OS_MONTHS, STATUS)~', single_gene))
  surv_uni_cox <- summary(coxph(formula, data = data))
  ph_hypothesis_p <- cox.zph(coxph(formula, data = data))$table[1,3]
  if (surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){  #get the pvalue
    single_cox_report <- data.frame('uni_cox_sig_genes'=single_gene,
                                    'beta'=surv_uni_cox$coefficients[,1],
                                    'Hazard_Ratio'=exp(surv_uni_cox$coefficients[,1]),
                                    'z_pvalue'=surv_uni_cox$coefficients[,5],
                                    'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
                                    'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
    single_cox_report
  }
 }
uni_cox_list <- lapply(gene_list[-(1:4)], uni_cox)
candidate_gene=do.call(rbind, uni_cox_list)


#########################  modele construit ##################
#perform the multi-variates cox regression using qualified genes.
formula_for_multivariate <- as.formula(paste0('Surv(OS_MONTHS, STATUS)~', paste(candidate_gene$uni_cox_sig_genes, sep = '', collapse = '+')))
multi_variate_cox <- coxph(formula_for_multivariate, data = data)
#check if variances are supported by PH hypothesis.
ph_hypo_multi <- cox.zph(multi_variate_cox)
#The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#Remove variances not supported by ph hypothesis and perform the 2nd regression.
formula_for_multivariate <- as.formula(paste0('Surv(OS_MONTHS, STATUS)~', paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05], sep = '', collapse = '+')))
multi_variate_cox_2 <- coxph(formula_for_multivariate, data = data)

################ forest plot #####################
correlation <- cor(data[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], method = 'pearson')
library('GGally')
ggpairs(data[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]], 
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())
library('rms')
vif <- rms::vif(multi_variate_cox_2)
#Some people said if the square root of VIF >2, they might be co-linear.
sqrt(vif) < 2

ggforest(model = multi_variate_cox_2, data = data, main = 'Hazard ratios of candidate genes', fontsize = 1)

########################### heatmap gene interesting (autres analyses, not very necessary) ##################################
load("~/Documents/DEGs_multivariates_DesEdg_withoutCov7.RData")
#choose one
geneList=read.xlsx("~/Documents/Genes_list.xlsx")
geneList=read.xlsx("~/Documents/signature_TLS.xlsx") 
row_annot=melt(t(geneList))
colnames(row_annot)=c("GenesNature","id","GeneName")
row_annot=na.omit(row_annot)
list_gene=unique(row_annot$GeneName)

sample_info_commun=Fun_annot("withoutCOV7")

norm_data=normalized_counts[list_gene,]
norm_data=na.omit(norm_data)
#INF
norm_data=norm_data[-which(rownames(norm_data)=="IFNA11P"),]
#rename according our samples
rownames(sample_info_commun)=colnames(norm_data) 

#row(gene) annotation
row_sample=row_annot[rownames(norm_data),]
row_sample=subset(row_sample,select=GenesNature)

colours_col <- list("Histology_type"=c("adenocarcinoma"="#6699CC",
                                       "adenosquamous carcinoma"="#40B0A6",
                                       "squamous cell carcinoma"="#D35FB7",
                                       "neuroendocrine carcinoma"="#FEFE62","X"="orange"),
                    "Age"=c("40-60"="#99DDFF","60-80"="#CCDDAA","X"="#E1BE6A"),
                    "Sex"=c("M"="#EE8866","F"="#AA4499","X"="slategrey"),
                    "Sample_Type"= c("COV" = "#0C7BDC" ,"Neg" ="#FFC20A"))


scale_data_deseq=norm_data%>% pheatmap:::scale_rows()
annot_degs=HeatmapAnnotation(df=sample_info_commun,which="col",col=colours_col,boxplot=anno_boxplot(norm_data))
annot_row_degs_deseq2=rowAnnotation(df=row_sample,bar=anno_barplot(rowMeans(norm_data)))
Heatmap(scale_data_deseq,name = "Z-score",cluster_rows=F,show_row_names = T,show_column_names = FALSE,top_annotation=annot_degs,left_annotation=annot_row_degs_deseq2,
        column_split = sample_info_commun$Sample_Type,column_title = "Heatmap of immunology process gene avec different patients", 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))