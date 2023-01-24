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
      numericInput(inputId = "pvalue",
                   label = "p-value",
                   step = 0.01,
                   value = 0.05,
                   min = 0.01,
                   max = 0.99),
      numericInput(inputId = "Treshold",
                   label = "p-value",
                   step = 0.5,
                   value = 2,
                   min = 0.5,
                   max = 4),
      selectInput(inputId = "palette",
                  label = "Color palette",
                  choices = rownames(brewer.pal.info))
    ),
    mainPanel(
      tabsetPanel(type = "pills",
        tabPanel(title = "Input data",
                           DT::DTOutput("raw_counts")),
        tabPanel("Results",
                dataTableOutput("DE_results")),
        tabPanel("Plots",
                column(6,
                  plotOutput("plot1"),
                  plotOutput("plot2")
                  ),
                column(6,
                  plotOutput("plot3"),
                  plotOutput("plot4")
                  
                ))
      )
    ),
  )

))
