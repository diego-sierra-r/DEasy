library(shiny)
library(DT)
library(magrittr)
library(dplyr)
options(shiny.maxRequestSize = 100 *
          1024^2)

# define colors needed to customize tables
set.seed(123)
colors <- c('#00BFFF','#1E90FF','#6495ED','#7B68EE',
            '#4169E1','#0000FF','#0000CD','#00008B',
            '#000080','#191970','#FFF8DC','#FFEBCD',
            '#FFE4C4','#FFDEAD','#F5DEB3','#DEB887',
            '#D2B48C','#BC8F8F','#F4A460','#DAA520',
            '#B8860B','#CD853F','#D2691E','#8B4513',
            '#A0522D','#A52A2A','#800000','#FFFFFF',
            '#FFFAFA','#F0FFF0','#F5FFFA','#F0FFFF',
            '#F0F8FF','#F8F8FF','#F5F5F5','#FFF5EE',
            '#F5F5DC','#FDF5E6','#FFFAF0','#808080',
            '#FFFFF0','#FAEBD7','#FAF0E6','#FFF0F5',
            '#FFE4E1','#DCDCDC','#D3D3D3','#C0C0C0',
            '#A9A9A9')
            

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
  interection_choose <- reactive({
    file2 <- vroom::vroom(input$sampleinfo$datapath, delim = ",")
    column <- file2[[input$interaction]] 
  })

  # Define outputs
  treatment_DT <- reactive(DT::datatable(sampleinfo_data(), # sample information table with colors according treatment
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
    
                             DT::formatStyle(
                               as.character(input$treatment),
                               backgroundColor = DT::styleEqual(
                                 levels = unique(treatment_choose()),
                                 values = sample(colors,
                                                 length(unique(treatment_choose(
                                                 ))),
                                                 replace = T)
                               )
                             ) 
  )
  
  output$sampleinfo <- DT::renderDataTable({
    if (input$interaction == "None") {
      treatment_DT()
    } else {
      treatment_DT() %>% DT::formatStyle(
        as.character(input$interaction),
        backgroundColor = DT::styleEqual(
          levels = unique(interection_choose()),
          values = sample(colors,
                          length(unique(interection_choose(
                          ))),
                          replace = T)
        )
      )
    }
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
