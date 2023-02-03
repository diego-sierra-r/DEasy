library(shiny)
library(tidyr)
library(DT)
library(DESeq2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

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

# define plot funtions


MDS_plot <- function(countData,
                    colData,
                    treatment,
                    interac = "None") {
  
  if (interac == "None") {
    design <- as.formula(paste0("~ ", treatment))
  } else if (is.null(interac)) {
    design <- as.formula(paste0("~ ", treatment))
  } else {
    design =  as.formula(paste0("~ ", treatment, "+ ", interac))
    
  }
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData),
                                colData = colData,
                                design = design)
  vsd_0 <- vst(dds, blind = F) # calcualte dispersion trend
  sampleDists <- dist(t(assay(vsd_0))) #Calculate distance matrix
  sampleDistMatrix <- as.matrix( sampleDists ) # Create distance matrix
  mdsData <- data.frame(cmdscale(sampleDistMatrix)) #perform MDS
  mds <- cbind(mdsData, as.data.frame(colData(vsd_0)))
  
  F_vr_M_DESeq2_MDS <-  ggplot(mds, aes(X1,X2,color=SEX)) +
    geom_label_repel(aes(label = rownames(mds)), size = 3) +
    geom_point(size=3) +
    scale_color_manual(values =  c("#B22222","#8B008B"),
                       labels = c("Female", "Male"),
                       name = "Sex") +
    labs(title = "Females vr Males DESeq2",
         x = "Dim 1",
         y = "Dim 2") +
    theme_classic2()
  return(F_vr_M_DESeq2_MDS)
}

# create funtions for DE with DESeq2


validate_row_cols  <- function(df_s, df_r) {
  if (identical(row.names(df_s), colnames(df_r)) == FALSE) {
    stop("ERROR: Sample IDs on sample information rows and Raw counts columns must be the same \n see example data. ")
  }
}


DE_DESeq2_design  <- function(countData,
                              colData,
                              treatment,
                              alpha,
                              threshold,
                              interac = NULL) {
  countData <- as.matrix(countData)
  treatment =  treatment
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = colData,
                                        design = treatment)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  res <- results(dds, alpha = alpha, lfcThreshold = threshold)
  return(as.data.frame(res))
}

DE_DESeq2_main <- function(
    countData,
    colData,
    treatment,
    interac = NULL,
    alpha,
    threshold) {
  
  validate_row_cols(df_s =  colData,
                    df_r = countData)
  
  countData <- as.matrix(countData)
  
  if (is.null(interac)) {
    
    treatment =  as.formula(paste0("~ ", treatment))
    
    res <- DE_DESeq2_design(countData =  countData,
                            colData =   colData,
                            treatment =  treatment,
                            alpha = alpha,
                            threshold = threshold)
    
  } else if (interac == "None") {
    
    treatment =  as.formula(paste0("~ ", treatment))
    res <- DE_DESeq2_design(countData =  countData,
                            colData =   colData,
                            treatment =  treatment,
                            alpha = alpha,
                            threshold = threshold)
  } else {
    treatment =  as.formula(paste0("~ ", treatment, "+ ", interac))
    res <- DE_DESeq2_design(countData =  countData,
                            colData =   colData,
                            treatment =  treatment,
                            interac =  interac,
                            alpha = alpha,
                            threshold = threshold)
    
  }
  res <- as.data.frame(res)
  
  res <- res %>% 
    mutate(., Difference =(case_when(.data$log2FoldChange >= threshold & .data$padj <= 0.05 ~ "UP",
                                     .data$log2FoldChange <= -threshold & .data$padj <= 0.05 ~ "DOWN",
                                     .data$log2FoldChange <= threshold | .data$log2FoldChange >= 2 & .data$padj >0.05 ~ "Not significant"))) %>% drop_na()
  return(res)
}




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
    file <- read.csv(input$raw_counts$datapath, row.names = 1)
  })
  sampleinfo_data <- reactive({
    req(input$sampleinfo)
    file2 <- read.csv(input$sampleinfo$datapath, row.names = 1)
  })
  treatment_choose <- reactive({
    req(input$sampleinfo)
    file2 <- read.csv(input$sampleinfo$datapath)
    column <- file2[[input$treatment]] %>%  as.factor()
  }) 
  interection_choose <- reactive({
    req(input$sampleinfo)
    file2 <- read.csv(input$sampleinfo$datapath)
    column <- file2[[input$interaction]] %>%  as.factor()
  })

  # Define outputs
  treatment_DT <- reactive(
    DT::datatable(sampleinfo_data(), 
                                # sample information table with colors according treatment
                                rownames = TRUE,
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
      rownames = TRUE,
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
  
  results_DE_DESeq2 <- reactive({
    DE_DESeq2_main(
      countData = raw_counts_data(),
      colData = sampleinfo_data(),
      treatment = input$treatment,
      interac = input$interaction,
      alpha = 0.05,
      threshold = 2
    )
    
     
  })
  
  output$DE_results <-  renderDataTable({results_DE_DESeq2()})
  
  output$plot1 <- renderPlot({
    if (input$treatment == input$interaction) {
      validate("Treatment an interaction can't be de same")
    }
    MDS_plot(
      countData = raw_counts_data(),
      colData = sampleinfo_data(),
      treatment =  input$treatment,
      interac = input$interaction
    )
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
