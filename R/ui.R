library(shiny)
library(RColorBrewer)
library(waiter)



# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  titlePanel(title = "DEasy", windowTitle = "DEasy!"),
  sidebarLayout(
    sidebarPanel(
      titlePanel("Upload files"),
      fileInput(inputId = "sampleinfo", label = "Sample information", ),
      fileInput(inputId = "raw_counts", label = "Raw counts"),
      titlePanel("Parameters"),
      selectInput(
        inputId = "pgk",
        label = "Differential expression package to use:",
        choices = c("DESeq2", "edgeR"),
        selected = "DESeq2"
      ),
      selectInput(
        inputId = "treatment",
        label = "Treatment",
        choices = NULL
      ),
      selectInput(
        inputId = "interaction",
        label = "Add interaction (Optional)",
        choices = NULL
      ),
      numericInput(
        inputId = "pvalue",
        label = "p-value",
        step = 0.01,
        value = 0.05,
        min = 0.01,
        max = 0.99
      ),
      numericInput(
        inputId = "treshold",
        label = "Threshold",
        step = 0.5,
        value = 0,
        min = 0.5,
        max = 4
      ),
      # main action button
      actionButton(
        inputId = "run",
        label = "Run",
        icon = icon("dna", lib = "font-awesome"),
        class = "btn-success btn-lg btn-block"
      )
    ),
    ### tab panels
    mainPanel(
      tabsetPanel(
        type = "pills",
        navbarMenu(
          title = "Preview",
          tabPanel(
            title = "Sample information",
            DT::DTOutput("sampleinfo")
          ),
          tabPanel(
            title = "Raw counts",
            DT::DTOutput("raw_counts")
          )
        ),
        tabPanel(
          "Results",
          column(
            12,
            DT::DTOutput("DE_results")
          ),
          downloadButton(
            "downloadR",
            "Download",
            icon = icon("download"),
            class = "btn-info"
          )
        ),
        tabPanel(
          "Plots",
          column(
            6,
            plotOutput("plot1"),
            downloadButton(
              "downloadP1",
              "Download",
              icon = icon("download"),
              class = "btn-info"
            ),
            plotOutput("plot2"),
            downloadButton(
              "downloadP2",
              "Download",
              icon = icon("download"),
              class = "btn-info"
            )
          ),
          column(
            6,
            plotOutput("plot3"),
            downloadButton(
              "downloadP3",
              "Download",
              icon = icon("download"),
              class = "btn-info"
            ),
            plotOutput("plot4"),
            textInput(
              inputId = "geneID",
              label = "gene ID",
              width = "400px",
              placeholder = "Write a geneID from raw counts",
              value = NULL
            ),
            downloadButton(
              "downloadP4",
              "Download",
              icon = icon("download"),
              class = "btn-info"
            )
          )
        ),
        tabPanel(
          "Help",
          column(
            11,
            htmlOutput(outputId = "readme")
          )
        )
      )
    ),
  )
))

