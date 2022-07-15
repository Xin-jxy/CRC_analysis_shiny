server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  
  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = read.csv(input$upload$datapath, sep = ",",row.names = 1),
           tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
    
  })
  output$summary <- renderPrint({
    if (is.null(input$upload)){return()}
    else { return(summary(data()))
    }
  })
  output$table <- renderDataTable({
    if (is.null(input$upload)){return()}
    else { return(cbind(names = rownames(data()), data()))
      
    }
  })
  
  output$volcanoplot<-renderPlotly({
    ggplotly(EnhancedVolcano(data(),
                             lab = rownames(data()),
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
                             title = "Volcano plot without sample POS with Deseq2",
                             legendLabels = c("NS", expression(Log[2] ~ FC), "adjusted p-value", expression(p - adj ~ and
                                                                                                            ~ log[2] ~ FC)),
                             ylab = bquote(~-Log[10] ~ italic(Padj)),
                             axisLabSize = 18)+aes(x=log2FoldChange,y=-log10(padj)),
             height = 1000, width=600)
    
  })
  observeEvent(data(),{ 
    data_heatmap <- reactive({
      
      if (is.null(data())){return()}
      else {
        data=data()
        getname<-function(input){
          if (!is.null(input$datapath)) {
            # Extract file name (additionally remove file extension using sub)
            sub=sub(".csv$", "", basename(input$name))
            getName=strsplit(sub, split = "_")[[1]]
            return(getName[length(getName)])
          } else {
            return(NULL)
          } }
        Names=getname(input$upload)
        NormData=Fun_loadNormFile(Names)
        samplename=Fun_annot("withoutCOV7")
        Fun_heatmap(data,NormData,samplename) }
    })
    makeInteractiveComplexHeatmap(input, output, session,data_heatmap(),"heatmap")
    
  })
  
  ##################### the part of Go enrichssement #####################
  data_go <- reactive({
    dataDEGs=data()
    go=Fun_enrichmment(dataDEGs)
    dataDEGs_up=dataDEGs[dataDEGs$regulated=="up",]
    go_up=Fun_enrichmment(dataDEGs_up)
    dataDEGs_down=dataDEGs[dataDEGs$regulated=="down",]
    go_down=Fun_enrichmment(dataDEGs_down)
    return(list(go,go_up,go_down))
  })
  output$dotplot_all<-renderPlotly({
    go=data_go()[[1]]
    cat=dim(go)[1]
    ggplotly(dotplot(go,showCategory=cat,split="ONTOLOGY",title=paste("GO enrichment of pathway of DEGs",sep=" "))+facet_grid(ONTOLOGY~.,scale="free"),height = 1000, width = 600)%>% layout(height = 800, width = 800)
  })
  
  output$dotplot_up<-renderPlotly({
    go=data_go()[[2]]
    cat=dim(go)[1]
    ggplotly(dotplot(go,showCategory=cat,split="ONTOLOGY",title=paste("GO enrichment of pathway of DEGs up regulated",sep=" "))+facet_grid(ONTOLOGY~.,scale="free"),height = 1000, width = 600)%>% layout(height = 800, width = 800)
  })
  
  output$dotplot_down<-renderPlotly({
    go=data_go()[[3]]
    cat=dim(go)[1]
    ggplotly(dotplot(go,showCategory=cat,split="ONTOLOGY",title=paste("GO enrichment of pathway of DEGs down regulated",sep=" "))+facet_grid(ONTOLOGY~.,scale="free"),height = 1000, width = 600)%>% layout(height = 800, width = 800)
  })
  
  output$barplot_all<-renderPlotly({
    go=data_go()[[1]]
    cat=dim(go)[1]
    ggplotly(barplot(go,showCategory=cat, split="ONTOLOGY",title="GO enrichment of DEGs split with ontology")+ facet_grid(ONTOLOGY~.,scale="free"),height = 1000, width = 600)%>% layout(height = 800, width = 800)
    
  })
  
  output$barplot_up<-renderPlotly({
    go=data_go()[[2]]
    cat=dim(go)[1]
    ggplotly(barplot(go,showCategory=cat, split="ONTOLOGY",title="GO enrichment of DEGs up regulated split with ontology")+ facet_grid(ONTOLOGY~.,scale="free"),height = 1000, width = 600)%>% layout(height = 800, width = 800)
    
  })
  
  output$barplot_down<-renderPlotly({
    go=data_go()[[3]]
    cat=dim(go)[1]
    ggplotly(barplot(go,showCategory=cat, split="ONTOLOGY",title="GO enrichment of DEGs down regulated split with ontology")+ facet_grid(ONTOLOGY~.,scale="free"),height = 1000, width = 600)%>% layout(height = 800, width = 800)
    
  })
  output$cnetplot_all<-renderPlot({
    dataDEGs=data()
    go=data_go()[[1]]
    Fun_cnetplot(dataDEGs,go)
  })
  output$cnetplot_up<-renderPlot({
    dataDEGs=data()
    go=data_go()[[2]]
    Fun_cnetplot(dataDEGs,go)
  })
  output$cnetplot_down<-renderPlot({
    dataDEGs=data()
    go=data_go()[[3]]
    Fun_cnetplot(dataDEGs,go)
  })
  output$treeplot_all<-renderPlot({
    go=data_go()[[1]]
    cat=dim(go)[1]
    edox=pairwise_termsim(go, showCategory = cat)
    treeplot(edox, hclust_method = "average")
  })
  output$treeplot_up<-renderPlot({
    go=data_go()[[2]]
    cat=dim(go)[1]
    edox=pairwise_termsim(go, showCategory = cat)
    treeplot(edox, hclust_method = "average")
  })
  output$treeplot_down<-renderPlot({
    go=data_go()[[3]]
    cat=dim(go)[1]
    edox=pairwise_termsim(go, showCategory = cat)
    treeplot(edox, hclust_method = "average")
  })
  
  output$heatplot_all<-renderPlotly({
    go=data_go()[[1]]
    cat=dim(go)[1]
    genelist=data()$log2FoldChange
    names(genelist)=rownames(data())
    edox=pairwise_termsim(go, showCategory = cat)
    heatplot(edox, foldChange=genelist, showCategory=cat)
  })
  output$heatplot_up<-renderPlotly({
    go=data_go()[[2]]
    cat=dim(go)[1]
    genelist=data()$log2FoldChange
    names(genelist)=rownames(data())
    edox=pairwise_termsim(go, showCategory = cat)
    heatplot(edox, foldChange=genelist, showCategory=cat)
  })
  output$heatplot_down<-renderPlotly({
    go=data_go()[[3]]
    cat=dim(go)[1]
    genelist=data()$log2FoldChange
    names(genelist)=rownames(data())
    edox=pairwise_termsim(go, showCategory = cat)
    heatplot(edox, foldChange=genelist, showCategory=cat)
  })
  
  ##################### the part of gseGo enrichssement #####################
  data_gsego <- reactive({
    getname<-function(input){
      if (!is.null(input$datapath)) {
        # Extract file name (additionally remove file extension using sub)
        sub=sub(".csv$", "", basename(input$name))
        getName=strsplit(sub, split = "_")[[1]]
        return(getName[length(getName)])
      } else {
        return(NULL)
      } }
    Names=getname(input$upload)
    ele_genelist=Fun_genelist(Names)
    genelist=ele_genelist[[1]]
    data_deseq=ele_genelist[[2]]
    gsego=Fun_gseGO(genelist)
    return (list(gsego,data_deseq))
  })
  output$dotplot1<-renderPlotly({
    go=data_gsego()[[1]]
    cat=dim(go)[1]
    ggplotly(dotplot(go,showCategory=cat, title=paste("gseGO enrichment of pathway of DEGs")),height = 800, width = 600)%>% layout(height = 800, width = 800)
  })
  
  output$cnetplot1<-renderPlot({
    data=data_gsego()[[2]]
    go=data_gsego()[[1]]
    Fun_cnetplot(data,go)
  })
  
  output$treeplot1<-renderPlot({
    go=data_gsego()[[1]]
    cat=dim(go)[1]
    edox=pairwise_termsim(go, showCategory = cat)
    treeplot(edox, hclust_method = "average")
  })
  
  output$heatplot1<-renderPlot({
    go=data_gsego()[[1]]
    cat=dim(go)[1]
    genelist=data_gsego()[[2]]$log2FoldChange
    names(genelist)=rownames(data())
    edox=pairwise_termsim(go, showCategory = cat)
    heatplot(edox, foldChange=genelist, showCategory=cat)
  })
  ##################### the part of GESA enrichssement #####################
  data_gesa <- reactive({
    gesa=data_gsego()[[2]]
    df_symbol=normTransSymbol(gesa)
    geneList_gesa <- df_symbol$log2FoldChange
    names(geneList_gesa) <- rownames(df_symbol)
    geneList_gesa_sort=sort(geneList_gesa,decreasing = T)
    gesa=Fun_gesa(geneList_gesa_sort)
    return(list(gesa,df_symbol))
  })
  output$dotplot2<-renderPlot({
    go=data_gesa()[[1]]
    cat=dim(go@result)[1]
    dotplot(go,showCategory=cat, title=paste("GESA enrichment of pathway of DEGs"))
  })
  
  output$cnetplot2<-renderPlot({
    data=data_gesa()[[2]]
    go=data_gesa()[[1]]
    Fun_cnetplot(data,go)
  })
  
  output$treeplot2<-renderPlot({
    go=data_gesa()[[1]]
    cat=dim(go@result)[1]
    edox=pairwise_termsim(go, showCategory = cat)
    treeplot(edox, hclust_method = "average")
  })
  
  output$heatplot2<-renderPlot({
    go=data_gesa()[[1]]
    cat=dim(go@result)[1]
    genelist=data_gesa()[[2]]$log2FoldChange
    names(genelist)=rownames(data())
    edox=pairwise_termsim(go, showCategory = cat)
    heatplot(edox, foldChange=genelist, showCategory=cat)
  })
  
  ##################### the part of MCPcounter #####################
  
  ########## violin plot ########
  output$violin<-renderPlot({
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withoutCov7.csv",header=T,row.names=1)
    data_mcp=Fun_mcpcounter(df_normalized)
    df_cov7=data_mcp[[1]]
    Fun_violin(df_cov7)
    
  })
  
  output$violin1<-renderPlot({
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withinCov7.csv",header=T,row.names=1)
    data_mcp=Fun_mcpcounter(df_normalized)
    df_cov7=data_mcp[[1]]
    Fun_violin(df_cov7)
    
  })
  
  ###### heatmap for immunology analysis
  output$HeatMcp<-renderPlot({ 
    #preparation of gene
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withoutCov7.csv",header=T,row.names=1)
    sample_info_commun=Fun_annot("withoutCOV7")
    data_mcp=Fun_mcpcounter(df_normalized)
    df_norm=data_mcp[[2]]
    genes=data_mcp[[3]]
    Fun_heatmcp(df_norm,genes,sample_info_commun)
  })
  
  output$HeatMcp1<-renderPlot({ 
    #preparation of gene
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withinCov7.csv",header=T,row.names=1)
    sample_info_commun=Fun_annot("withCOV7")
    data_mcp=Fun_mcpcounter(df_normalized)
    df_norm=data_mcp[[2]]
    genes=data_mcp[[3]]
    Fun_heatmcp(df_norm,genes,sample_info_commun)
  })
  
  ##################### the part of TMEConsensus #####################
  ########## violin plot ########
  output$violin2<-renderPlot({
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withoutCov7.csv",header=T,row.names=1)
    data_mcp=Fun_TME(df_normalized)
    Fun_violin(data_mcp)
    
  })
  
  output$violin3<-renderPlot({
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withinCov7.csv",header=T,row.names=1)
    data_mcp=Fun_TME(df_normalized)
    Fun_violin(data_mcp)
    
  })
  
  ###### heatmap for immunology analysis
  output$HeatTME<-renderPlot({ 
    #preparation of gene
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withoutCov7.csv",header=T,row.names=1)
    logTPM <- function(x) {return(log2(x+1))}
    df_norm=df_normalized %>% mutate_if(is.numeric, logTPM)
    sample_info_commun=Fun_annot("withoutCOV7")
    genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
    Fun_heatmcp(df_norm,genes,sample_info_commun)
  })
  
  output$HeatTME1<-renderPlot({ 
    #preparation of gene
    df_normalized=read.csv("~/Documents/analysis_script_RNAseq/data/Normalized_TPM_withinCov7.csv",header=T,row.names=1)
    logTPM <- function(x) {return(log2(x+1))}
    df_norm=df_normalized %>% mutate_if(is.numeric, logTPM)
    sample_info_commun=Fun_annot("withCOV7")
    genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
    Fun_heatmcp(df_norm,genes,sample_info_commun)
  })
}