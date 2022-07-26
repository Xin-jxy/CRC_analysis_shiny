## CRC_analysis_shiny

A project of lung cancer patients who have the Covid, that uses mainly bulk RNA-seq analysis and some public data, which contains different analyse scripts, this is the instruction of the utilization of script, the data is confidential. 
#### Authors
[@Xin-jxy](https://github.com/Xin-jxy)

## BASH
### RNAseq UP-stream analysis
The alignment and annotation process is done at server, find the data by taping `cd /data/data_cremer_lab/LungT_Covid_XinAnalysis` in Terminal. See what we have in this repertory with the command `ls`. (detail about the alignment and annotation, please read my report, the command bash is available in the repertory as well)

## R
### A script individual for the convenience of shiny or other utilizations

`Universal_Function.R` is a script that contains diverse functions for avoiding the repetition of code, which could be modified in different situations, the detail of utilization is written in the script.
 
### RNAseq DOWN-stream analysis

`analyse_DEGs_Deseq2-EdgeR.R` is a script that contains mainly the process of gene differentially expressed (DEGs) analysis with two toolkits DESeq2 and edgeR. 

`Enrichessement-analysis.R` is a script of enrichment analysis, containing the GO, KEGG, GESA analysis.

`mcp_counter.R` and `Consensus_TME.R` are two scripts similary by using the toolkit MCPcounter and ConsensusTME. ConsensusTME can use design signature genes for calculating the score.

`WGCNA_analyse.R` is a script that we use for WGCNA analysis, the notices are written in the script.

### other analysis with public data
`side_analysis_immunology.R` contains the analysis of survie by using the database of TCGA, don't need to download the data, a package is availiable to query directly. The data used is `Lung Adenocarcinoma (TCGA, Nature 2014)- 230 samples`
There is a side analysis of heatmap, is used to find out the correaltion of selected gene with differents patients.

`scRNAseq_luad.R` uses the public data from [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) with `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz`, the code bash for extracting only the Tumor tissue is accomplished with Terminal:

`zcat GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz | awk 'BEGIN { FS=OFS="\t" }
    NR==1 {
        for (inFldNr=1; inFldNr<=NF; inFldNr++) {
            if ($inFldNr ~ /LUNG_T/) {
                out2inFldNr[++numOutFlds] = inFldNr
            }
        }
    }
    {
        for (outFldNr=1; outFldNr<=numOutFlds; outFldNr++) {
            inFldNr = out2inFldNr[outFldNr]
            printf "%s%s", $inFldNr, (outFldNr<numOutFlds ? OFS : ORS)
        }
    }' |gzip > GSE131907_LungTumor_log2TPM.txt.gz `

**Attention** This script can be used only in serveur and load the whole data sample of GSE131907 may full all the space with server. Nevertheless, this selection of data reduces the saturation of serveur, but for utilization of others, please clean the object after each utilization by using the command `rm(list=ls())` and then `.rs.restartR()` in console.

The data and figure are available at repertory`/data/data_cremer_lab/LungT_Covid_XinAnalysis` as well in our server.