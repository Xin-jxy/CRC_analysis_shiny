X <- c("plyr", "dplyr", "tm", "readxl", "wordcloud", "SnowballC", "stringdist", "tidytext",
       "rmarkdown", "knitr", "quanteda", "reshape", "stringr", "RecordLinkage", "plotly",
       "data.table", "rvest",  "shiny", "shinydashboard", "shinyWidgets", "DT","shinythemes","InteractiveComplexHeatmap",
       "ggplot2" ,"clusterProfiler","pheatmap" ,"edgeR" , "statmod","DESeq2","NMF","ggbeeswarm",
       "genefilter","pheatmap","ade4","viridis","tidyverse","dplyr","tximeta","tximport","openxlsx",
       'EnhancedVolcano',"ggrepel","biomaRt","reshape","RColorBrewer","VennDiagram","ComplexHeatmap",
       "crosstalk","ggvenn","DT","curl","ConsensusTME","org.Hs.eg.db","GSVA","limma","tidyr",
       "enrichplot","fgsea","msigdbr","ggnewscale","devtools","MCPcounter","gridExtra","ggsignif") 


lapply(X, FUN = function(X){
  do.call("library", list(X))
})

ui <- dashboardPage(
  dashboardHeader(title = "Data visualization of DEGs, impact of Covid in lung cancer patients"),
  dashboardSidebar(
    sidebarMenu(
      ## Tab 1 -- Specify Task and View Raw Data Files
      menuItem("Select Task And Upload Files", tabName = "task", icon = icon("file-text-o")),
      ## Tab 2 -- Heatmap and volcanoplot of DEGs
      menuItem("DEGs analysis", tabName = "raw", icon = icon("file-text-o")),
      ## Tab 3 -- Enrichssement analysis
      menuItem("Enrichissement analysis", tabName = "processed", icon = icon("file-text-o")),
      ## Tab 4 -- Immunology analysis
      menuItem("Immunology analysis", tabName = "more", icon = icon("file-text-o"))
      
    )),
  dashboardBody(
    tabItems(
      ### Specify Task & Upload Files Tab
      tabItem(tabName = "task",
              fileInput("upload", label="Please load data for visualization", accept = c(".csv", ".tsv")),
              helpText(paste("Please upload a file.  Supported file types are:  .csv and .tsv")),
              mainPanel(
                
                # Output: Tabset w/ plot, summary, and table ----
                tabsetPanel(type = "tabs",
                            tabPanel("Table", dataTableOutput("table")),
                            tabPanel("DataSummary", verbatimTextOutput("summary"))
                )
              )
              
      ), # close first tabItem
      
      tabItem(tabName = "raw",
              helpText(paste("This tab displays the analysis of DEGs.")),
              # Main panel for displaying outputs ----
              mainPanel(
                
                # Output: Tabset w/ plot, summary, and table ----
                tabsetPanel(type = "tabs",
                            tabPanel("Heatmap", InteractiveComplexHeatmapOutput("heatmap")),
                            tabPanel("VolcanoPlot", plotlyOutput("volcanoplot"))
                )
              )
      ), # close tabItem
      
      tabItem(tabName = "processed",
              helpText(paste("This tab displays the enrichissement analysis.")),
              # Main panel for displaying outputs ----
              navbarPage("Enrichissement analysis", collapsible = TRUE, inverse = TRUE, theme = shinytheme("readable"),
                         tabPanel("GO enrichment analysis",
                                  fluidPage(
                                    tabsetPanel(
                                      tabPanel("Network of pathway",fluidRow(
                                        h4("Network for all DEGs"),
                                        plotOutput("cnetplot_all",width = "800px", height = "1000px"),
                                        h4("Network for DEGs up regulated"),
                                        plotOutput('cnetplot_up',width = "800px", height = "1000px"),
                                        h4("Network for DEGs down regulated"),
                                        plotOutput('cnetplot_down',width = "800px", height = "1000px"))),
                                      tabPanel("dotplot",fluidRow(
                                        plotlyOutput("dotplot_all",width = "800px", height = "1000px"),
                                        plotlyOutput('dotplot_up',width = "800px", height = "1000px"),
                                        plotlyOutput('dotplot_down',width = "800px", height = "1000px"))),
                                      tabPanel("barplot",fluidRow(
                                        plotlyOutput("barplot_all",width = "800px", height = "1000px"),
                                        plotlyOutput('barplot_up',width = "800px", height = "1000px"),
                                        plotlyOutput('barplot_down',width = "800px", height = "1000px"))),
                                      tabPanel("relation of pathway plot",fluidRow(
                                        h4("Treeplot for all DEGs"),
                                        plotOutput("treeplot_all",width = "800px", height = "1000px"),
                                        h4("Treeplot for DEGs up regulated"),
                                        plotOutput('treeplot_up',width = "800px", height = "1000px"),
                                        h4("Treeplot for DEGs down regulated"),
                                        plotOutput('treeplot_down',width = "800px", height = "1000px"))),
                                      tabPanel("Heatplot of pathway",fluidRow(
                                        h4("Heatplot for all DEGs"),
                                        plotlyOutput("heatplot_all",width = "800px", height = "1000px"),
                                        h4("Heatplot for DEGs up regulated"),
                                        plotlyOutput('heatplot_up',width = "800px", height = "1000px"),
                                        h4("Heatplot for DEGs down regulated"),
                                        plotlyOutput('heatplot_down',width = "800px", height = "1000px")))
                                    ))),
                              tabPanel("gseGO enrichment analysis",
                                  fluidPage(
                                    tabsetPanel(
                                      tabPanel("Network of pathway",
                                        plotOutput("cnetplot1",width = "800px", height = "1000px")),
                                      tabPanel("dotplot",
                                        plotlyOutput("dotplot1",width = "600px", height = "1000px")),
                                      tabPanel("relation of pathway plot",
                                        plotOutput("treeplot1",width = "800px", height = "1000px")),
                                      tabPanel("Heatplot of pathway",
                                        plotOutput("heatplot1",width = "600px", height = "1000px"))
                                    ))), 
                         tabPanel("GESA enrichment analysis",
                                  fluidPage(
                                    tabsetPanel(
                                      tabPanel("Network of pathway",
                                               plotOutput("cnetplot2",width = "800px", height = "1000px")),
                                      tabPanel("dotplot",
                                               plotOutput("dotplot2",width = "1000px", height = "1000px")),
                                
                                      tabPanel("relation of pathway plot",
                                               plotOutput("treeplot2",width = "800px", height = "1000px")),
                                      tabPanel("Heatplot of pathway",
                                               plotOutput("heatplot2",width = "600px", height = "1000px"))
                                    )))
              )
              
      ), # close tabItem
      tabItem(tabName = "more",
              helpText(paste("This tab displays the analysis of Immunology.")),
              # Main panel for displaying outputs ----
              mainPanel(
                
                # Output: Tabset w/ plot, summary, and table ----
                tabsetPanel(type = "tabs",
                            tabPanel("MCPcounter analysis",
                                     fluidPage(
                                       tabsetPanel(
                                         tabPanel("violinplot",fluidRow(
                                                  h4("Violin plots of Cell Types without COV7"),
                                                  plotOutput("violin",width = "1000px", height = "1000px"),
                                                  h4("Violin plots of Cell Types within COV7"),
                                                  plotOutput("violin1",width = "1000px", height = "1000px"))),
                                         tabPanel("Heatmap",fluidRow(
                                                  h4("Heatmap of Cell type for samples without COV7"),
                                                  plotOutput("HeatMcp",width = "600px", height = "1000px"),
                                                  h4("Heatmap of Cell type for samples within COV7"),
                                                  plotOutput("HeatMcp1",width = "600px", height = "1000px")))
                                       ))),
                            tabPanel("TMEconcensus", 
                                     fluidPage(
                                       tabsetPanel(
                                         tabPanel("violinplot",fluidRow(
                                           h4("Violin plots of Cell Types without COV7"),
                                           plotOutput("violin2",width = "1000px", height = "1000px"),
                                           h4("Violin plots of Cell Types within COV7"),
                                           plotOutput("violin3",width = "1000px", height = "1000px"))),
                                         tabPanel("Heatmap",fluidRow(
                                           h4("Heatmap of Cell type for samples without COV7"),
                                           plotOutput("HeatTME",width = "600px", height = "1000px"),
                                           h4("Heatmap of Cell type for samples within COV7"),
                                           plotOutput("HeatTME1",width = "600px", height = "1000px")))
                                       )))
                )
              )
      )# close tabItem
    ) # close tabItems
  ) 
)



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

shinyApp(ui, server)