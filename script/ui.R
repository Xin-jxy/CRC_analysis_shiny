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
