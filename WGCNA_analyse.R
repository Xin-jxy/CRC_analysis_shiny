source("~/Documents/analysis_script_RNAseq/script/Universal_Function.R")

#Liste des packages necessaires:
packages=c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "igraph",
           "fastcluster", "dynamicTreeCut", "survival","NMF","ggbeeswarm","clusterProfiler",
           "GO.db", "preprocessCore", "impute","WGCNA","reshape2","stringr","dplyr","tximport") 
#Lancer la fonction de chargement/installation: 
fun_packages(packages)
#BiocManager::install(c("GO.db", "preprocessCore", "impute","WGCNA"))

Fun_data_normlized_FPKM<-function (dir){
### load the data ####
setwd(dir)
filelist <- list.files(dir) 
# get names of dataset
element=loadData(filelist)
DataList=element[[2]]

setwd("~")
#get normalized data here is FPKM 
data_norm=normData(DataList,"FPKM")

# the data will be SYMBOL setted

colnames(data_norm)=gsub("Pos",'Neg',colnames(data_norm),fixed=T)
data_norm=data_norm[which(rowSums(data_norm)>10),]
logFPKM <- function(x) {return(log2(x+1))}
df_norm=data_norm %>% mutate_if(is.numeric, logFPKM)

#transfert_to_symbol_names
df_norm=normTransSymbol(df_norm)
return(df_norm)}

df_norm=Fun_data_normlized_FPKM("~/Documents/RSEM_res")
###########################################
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# !!!! Attention, the type of network correlation, recommand "signed" or "signed hybrid", but notice it will lose the data negatively correlated
type = "unsigned"

# calcule of correlation
# recommand biweight mid-correlation & bicor
# corType: pearson or bicor
corType="bicor"
corFnc = ifelse(corType=="bicor", cor, bicor)
# need to be setted specially when multivariates settled  
maxPOutliers = ifelse(corType=="bicor",1,0.05)

# multivariates
robustY = ifelse(corType=="bicor",T,F)

######################### filtering  ######################

#filter the variance 75% of gene, at least MAD>0.01
#or don't filter, make MAD>0,
#the amount of calcul will be reduced after filtering and lose some informations

m.mad <- apply(df_norm,1,mad)
dataExprVar <- df_norm[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## gene in row and sample in col
dataExpr <- as.data.frame(t(dataExprVar))

## verification of leak data
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
head(dataExpr)[,1:8]

#see if there is outlier samples
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

#begin of analysis
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# R-square=0.85
abline(h=0.85,col="red")

# Soft threshold and connexetion
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")

power = sft$powerEstimate

#there is some settle down that we could use when there is no power satisfaited

################## construction of network ######################
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "~/normlog2FPKM_WGCNA_0.25",
                       verbose = 3)
table(net$colors)

#load directly the data of net
load("~/Documents/0.25_log2FPKM_WGCNA_symbol.RData")

## grey is non-categorized
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# if don't satisfait，recutBlockwiseTrees，save times
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# module eigengene
MEs = net$MEs

### rename the colnames for avoiding the recalculating
### 
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# correlation heatmap with different module 
# marDendro/marHeatmap for setting the distance of edges
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)

#chose hub gene(the apparition of hub gene is just one principal)
hubs = chooseTopHubInEachModule(dataExpr, colorh=moduleColors, power=power, type=type)
hubs
list_hub_module=data.frame(hubs)
#write.csv(list_hub_module,"~/list_hub_module.csv")

############################################################################################################
############################### Only for the generation of heatmap for each gene, time consuming, not really necessary#######
############################################################################################################
TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
#load if you don't calculate the network but load
load("~/Documents/norm_log2(FPKM+1)_0.25_symbol.tom-block.1.RData", verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# Attention!!!clustering plot but time consuming, not adapt for large volume of data
TOMplot(plotTOM, net$dendrograms, moduleColors,
        main = "Network heatmap plot, all genes")

probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)


############################################################################################################
################################## module trait related to their phenotype #################################
############################################################################################################
datTraits=factor(c(rep("COV",6),rep("Neg",16)))
table(datTraits)

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
design=model.matrix(~0+ datTraits)
colnames(design)=levels(datTraits)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) ##value of ME with different colors (samples vs modules)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

#################### Display the correlation values within a heatmap plot and their phenotype #########################
png("~/Module-trait-relationships_0.25.png",width = 800,height = 1200,res = 120)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships 0.25"))
dev.off()
plot(rnorm(50), rnorm(50))

## the genes correlated to module with their phenotype
### calculate the correlation of module and phenotype

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
## each module and gene with their pearson correlation value
## MEs is the value of each module in each sample
## datExpr is the each quantification of sample in each gene
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

## only for the phenotype limite calculation
Cov = as.data.frame(design[,1])
names(Cov) = "Cov"
geneTraitSignificance = as.data.frame(cor(dataExpr, Cov, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Cov), sep="")
names(GSPvalue) = paste("p.GS.", names(Cov), sep="")

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.


##################### all the plot with the realtion of GS and module membership ##########################
colorlevels=unique(moduleColors)
pdf("GS vs. Module Membership.pdf",width = 14,height = 8)
par(mfrow=c(4,5))
par(mar = c(5,5,3,3))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (moduleColors==whichmodule)
  column = match(whichmodule, modNames)
  verboseScatterplot(abs(geneModuleMembership[restrict1, column]),
                     abs(geneTraitSignificance[restrict1, 1]), 
                     col=moduleColors[restrict1],
                     main=whichmodule, 
                     xlab = "Module Membership", ylab = "Gene Significance", abline = TRUE)
}
dev.off()

##################### one interested module and its hubgenes ###########################
module = "skyblue"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Covid",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


########################  get infos of all hub genes of desired module ################################
probes = colnames(dataExpr)
geneInfo0 = data.frame(Genes = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
##Order modules by their significance for trait_x
modOrder = order(-abs(cor(MEs, design, use = "p")));
##Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Cov))
geneInfo = geneInfo0[geneOrder, ] 

Module_Related_genes=geneInfo[which(geneInfo$moduleColor=="yellow"),]
write.csv(geneInfo, "module_genes.csv")

################################ Relation intra-module ##########################
################################ GS et connectivity #############################
ADJ1=abs(cor(dataExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)
colorlevels=unique(moduleColors)


pdf("GS vs. Connectivity.pdf",width = 14,height = 8)
par(mfrow=c(4,5))
par(mar = c(5,5,3,3))
for (i in c(1:length(colorlevels))){
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     geneTraitSignificance[restrict1, 1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
dev.off()

#################################  geneModuleMembership$MMyellow)>0.8 & abs(geneTraitSignificance)>0.2, hub genes #########
############ Attention!!! change MM.yellow name to MM.lightcyan or or other ########################
datKME=signedKME(dataExpr, MEs, outputColumnName="MM.")
hub_selected<- abs(datKME$MM.yellow)>0.8 & abs(geneTraitSignificance)>0.2
h=data.frame(hub_selected)
h$gene=rownames(h)
hub_selected=h[which(h$GS.Cov=="TRUE"),]
##verify the intersect of all genes and the significant gene intra-module
intersect(hub_selected$gene,Module_Related_genes$Genes)
write.csv(hub_selected, "hubgene_KME_yellow_significant.csv")

############################################################################################################################
############################ calculate to see the position of our own covid sample with other modules #######################
############################################################################################################################

adjacency = adjacency(dataExpr, power = power)

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
#claculate the geneTree
geneTree = hclust(as.dist(dissTOM), method = "average")

### melt the modules 
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

################## assembly the similarity module ####################
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged
#save(MEs,moduleLabels,moduleColors,geneTree,file="~/data_module_color_gene.RData")


################################ add our phenotype, here covid to observe their position with other modules ############

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, Cov))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)
par(cex = 0.9)
png("Eigengene-dendrogram.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# Plot the dendrogram
sizeGrWindow(6,6)
par(cex = 1.0)
## only the thee
png("Eigengene-dendrogram-hclust.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)

## only heatmap
png("Eigengene-adjacency-heatmap.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

############################################################################################################
############################### generation of cytoscape if needed, for the network visualization ##########################################
############################################################################################################

# Export the network into edge and node list files Cytoscape can read
# threshold is 0.5 for default, but we can adjust in cytoscape
#cyt = exportNetworkToCytoscape(TOM,
#edgeFile = paste("0.25Symbol_wgcna.edges.txt", sep=""),
#nodeFile = paste("0.25Symbol_wgcna.nodes.txt", sep=""),
#weighted = TRUE, threshold = 0,
#nodeNames = probes, nodeAttr = moduleColors)

