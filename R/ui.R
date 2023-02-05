library(shiny)
library(RColorBrewer)



# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  titlePanel(title = "DEasy", windowTitle = "DEasy!"),
  sidebarLayout(
    sidebarPanel(
      titlePanel("Upload files"),
      fileInput(inputId = "sampleinfo",label = "Sample imformation",),
      fileInput(inputId = "raw_counts",label = "Raw counts"),
      titlePanel("Parameters"),
      selectInput(inputId = "pgk",
                  label = "Differential expression package to use:",
                  choices = c("DESeq2","edgeR")),
      selectInput(inputId = "treatment", label = "Treatment",
                  choices = NULL),
      selectInput(inputId = "interaction", label = "Add interaction (Optional)",
                  choices = NULL),
      numericInput(inputId = "pvalue",
                   label = "p-value",
                   step = 0.01,
                   value = 0.05,
                   min = 0.01,
                   max = 0.99),
      numericInput(inputId = "treshold",
                   label = "Threshold",
                   step = 0.5,
                   value = 2,
                   min = 0.5,
                   max = 4),
      selectInput(inputId = "palette",
                  label = "Color palette",
                  choices = rownames(brewer.pal.info),
                  selected = "Set1"),
      # main action button
      actionButton(inputId = "run", label = "Run",
                   icon =  icon("dna", lib = "font-awesome"),
                   class = "btn-success btn-lg btn-block")
      
    ),
    ### tab panels
    mainPanel(
      tabsetPanel(
        type = "pills",
        navbarMenu(title = "Preview",
                             tabPanel(title = "Sample information", 
                                      DT::DTOutput("sampleinfo")),
                             tabPanel(title = "Raw counts",
                                      DT::DTOutput("raw_counts"))),
        tabPanel("Results",
                column(12,
                       dataTableOutput("DE_results")),
                       actionButton("downloadR","Download",
                                    icon = icon("download"),
                                    class = "btn-info")),
        tabPanel("Plots",
                column(6,
                  plotOutput("plot1"),
                  plotOutput("plot2")
                  ),
                column(6,
                  plotOutput("plot3"),
                  plotOutput("plot4"),
                  textInput(inputId = "geneID",
                            label =  "gene ID",width = "400px",
                            placeholder = "Write gene ID from raw counts",
                            value = NULL
                              )
                  
                ),
                actionButton("downloadP","Download",
                             icon = icon("download"),
                             class = "btn-info")),
        tabPanel("Help",
                 column(11,
                        htmlOutput(outputId = "readme")
                        )
                 )
      )
    ),     
  )

))
