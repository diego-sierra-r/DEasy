library(shiny)
library(DT)
library(magrittr)
library(dplyr)
options(shiny.maxRequestSize = 100 *
          1024^2)

# define colors needed to customize tables
set.seed(123)
colors <- c(
        '#CD5C5C','#F08080','#FA8072','#E9967A','#FFA07A',
        '#FFC0CB','#FFB6C1','#FF69B4','#FF1493','#FFA07A',
        '#FF7F50','#FF6347','#FF8C00','#FFA500','#FFD700',
        '#FFFF00','#FFFFE0','#FFFACD','#FFDAB9','#EEE8AA',
        '#F0E68C','#BDB76B','#E6E6FA','#D8BFD8','#DDA0DD',
        '#EE82EE','#DA70D6','#FF00FF','#FF00FF','#BA55D3',
        '#9370DB','#ADFF2F','#7FFF00','#7CFC00','#00FF00',
        '#32CD32','#98FB98','#90EE90','#00FA9A','#00FF7F',
        '#3CB371','#00FFFF','#00FFFF','#E0FFFF','#AFEEEE',
        '#7FFFD4','#40E0D0','#48D1CC','#00CED1','#FFF8DC',
        '#FFEBCD','#FFE4C4','#FFDEAD','#F5DEB3','#DEB887',
        '#D2B48C','#BC8F8F','#F4A460','#DAA520','#B8860B',
        '#CD853F','#D2691E')

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  #circular reference between sampleinfo and treatment select input
  observeEvent(input$sampleinfo,
               if (!is.null(input$sampleinfo$name)) {
                 updateSelectInput(inputId = "treatment",
                                   choices = colnames(sampleinfo_data()))
                 updateSelectInput(inputId = "interaction",
                                   choices = c("None",colnames(sampleinfo_data())))
               }
               )



  # load main input data
  raw_counts_data <- reactive({
    req(input$raw_counts)
    file <- vroom::vroom(input$raw_counts$datapath)
  })
  sampleinfo_data <- reactive({
    req(input$sampleinfo)
    file2 <- vroom::vroom(input$sampleinfo$datapath, delim = ",")
  })
  treatment_choose <- reactive({
    req(input$sampleinfo)
    file2 <- vroom::vroom(input$sampleinfo$datapath, delim = ",")
    column <- file2[[input$treatment]] %>%  as.factor()
  })
  interection_choose <- reactive({
    req(input$sampleinfo)
    file2 <- vroom::vroom(input$sampleinfo$datapath, delim = ",")
    column <- file2[[input$interaction]] %>%  as.factor()
  })

  # Define outputs
  treatment_DT <- reactive(
    DT::datatable(sampleinfo_data(), 
                                # sample information table with colors according treatment
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
    req(input$interaction)
    # valdidate treatment vr interaction (can't be the same)
    if (input$treatment == input$interaction) {
      validate("Treatment an interaction can't be de same")
    }
    if (length(levels(treatment_choose())) != 2 ) {
      validate("Treatment must contain only 2 levels")
    }
    if (ncol(raw_counts_data())-1 != nrow(sampleinfo_data())) {
      validate("The number of rows in sample information and columns in raw counts without genID must be de same  ")
    }
    
    ## add colors to preview/sampleinfo if the user selected an interaction
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
  # render raw counts
  output$raw_counts <- DT::renderDataTable({
    req(input$raw_counts)
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
  output$readme <- renderUI({
    tags$div(includeHTML("~/Documentos/R/my_pkgs/DEasy/README.html"))
  })    
 
}) 
