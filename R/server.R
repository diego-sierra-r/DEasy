library(shiny)
library(DT)
library(magrittr)
options(shiny.maxRequestSize = 100 *
          1024^2)

# define colors needed to customize tables
set.seed(123)
colors <- c("#F08080","#E9967A","#DC143C",
            "#FFC0CB","#FF69B4","#DB7093",
            "#FFA07A","#FF6347","#FF8C00",
            "#FFFF00","#FFFACD","#FFDAB9",
            "#D8BFD8","#EE82EE","#BA55D3",
            "#ADFF2F","#98FB98","#3CB371",
            "#00FFFF","#E0FFFF","#7FFFD4",
            "#48D1CC","#5F9EA0","#87CEFA")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  #circular reference to treatment select input
  observeEvent(input$sampleinfo,
               if (!is.null(input$sampleinfo$name)) {
                 updateSelectInput(inputId = "treatment",
                                   choices = colnames(sampleinfo_data()))
                 updateSelectInput(inputId = "interaction",
                                   choices = c("None",colnames(sampleinfo_data())))
               }
               )
  
  # load main imput data
  raw_counts_data <- reactive({
    req(input$raw_counts)
    file <- vroom::vroom(input$raw_counts$datapath)
  })
  sampleinfo_data <- reactive({
    req(input$sampleinfo)
    file2 <- vroom::vroom(input$sampleinfo$datapath, delim = ",")
  })
  treatment_choose <- reactive({
    file2 <- vroom::vroom(input$sampleinfo$datapath, delim = ",")
    column <- file2[[input$treatment]] 
  })
  
  # Define outputs
  output$sampleinfo <- DT::renderDataTable({
    DT::datatable(sampleinfo_data(),
                  rownames = FALSE,
                  options = list(
                    columnDefs = list(list(targets = '_all', width = "2px")),
                    columns.type = "num",
                    scrollCollapse = TRUE,
                    searching = FALSE,
                    lengthChange = TRUE,
                    ordering = FALSE,
                    paging = TRUE,
                    pageLength = 20,
                    scrollX = "200px",
                    scrollY = "530px")) %>% 

        DT::formatStyle(as.character(input$treatment),
                        backgroundColor = DT::styleEqual(
                        levels = unique(treatment_choose()),
                        values = sample(colors,
                                        length(unique(treatment_choose())),
                                        replace = T)
                        )
      )

  })

  output$raw_counts <- DT::renderDataTable({
    DT::datatable(
      raw_counts_data(),
      rownames = FALSE,
      options = list(
        columnDefs = list(list(targets = '_all', width = "2px")),
        columns.type = "num",
        scrollCollapse = TRUE,
        searching = FALSE,
        lengthChange = TRUE,
        ordering = FALSE,
        paging = TRUE,
        pageLength = 20,
        scrollX = "200px",
        scrollY = "530px")
    )
  })
    
  output$DE_results<- renderDataTable({
    mtcars
  })
  
  output$plot1 <- renderPlot({
    plot(1:3)
  })
  output$plot2 <- renderPlot({
    plot(3:6)
  })
  output$plot3 <- renderPlot({
    plot(6:9)
  })
  output$plot4 <- renderPlot({
    plot(9:12)
  })

})
