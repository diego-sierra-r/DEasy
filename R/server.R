library(shiny)
options(shiny.maxRequestSize = 100 *
          1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  raw_counts_data <- reactive({
    req(input$raw_counts)
    file <- vroom::vroom(input$raw_counts$datapath,delim = ",")
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
