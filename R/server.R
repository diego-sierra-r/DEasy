library(DESeq2)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(shiny)
library(pasilla)
library(ggpubr)
library(ggrepel)
library(rlang)
library(reshape2)
library(viridis)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(magrittr)
library(glue)
library(waiter)

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
    labs(title = "DESeq2",
         x = "Dim 1",
         y = "Dim 2") +
    theme_classic2()
  return(F_vr_M_DESeq2_MDS)
}

## MA-plot function 


MA_plot  <- function(df) {
  if (length(colnames(df)) == 7) {
    plot <- ggplot(data = df, 
                   aes(x = log2(baseMean),
                       y = log2FoldChange, 
                       color = Difference))  +
      geom_point(alpha = 0.75) +
      #        geom_hex(bins = 30) +
      labs(color = "Differentially expressed",
           fill = "Number of transcripts")+ 
      xlab("Log2 Base-mean")+ ylab("log2 Fold-change")+
      theme_classic2() +
      coord_flip() +
      ggtitle("DESeq2")
    
    return(plot)
  } else if (length(colnames(df)) == 6) {
    plot <- ggplot(data = df, 
                   aes(x = logCPM,
                       y = logFC, 
                       color = Difference))  +
      geom_point(alpha = 0.75) +
      #        geom_hex(bins = 30) +
      labs(color = "Differentially expressed",
           fill = "Number of transcripts")+ 
      xlab("Log2 Base-mean")+ ylab("log2 Fold-change")+
      theme_classic2() +
      coord_flip() +
      ggtitle("edgeR")
    return(plot)
  }
}

# single_GenID_plot
single_gen_plot <- function(geneID, countData, ColData, treatment, interaction = "None") {
  if (interaction == "None") {
    design = as.formula(paste0("~ ", treatment))
  } else if (is.null(interaction)) {
    design = as.formula(paste0("~ ", treatment))
  } else {
    design = as.formula(paste0("~ ", treatment, "+ ", interaction))
  }
  
  deseq2Data <-DESeqDataSetFromMatrix(countData = countData,
                                      colData = ColData,
                                      design = design)
  intgroup = c(colnames(deseq2Data@colData))
  otop2Counts <-
    plotCounts(
      deseq2Data,
      gene = as.character(geneID),
      intgroup = intgroup,
      returnData = TRUE
    )
  
  
  plot <- ggplot(otop2Counts,
                 aes(
                   x =  !! sym(names(otop2Counts)[2]),
                   y = count,
                   color = !! sym(names(otop2Counts)[2]),
                 )) + geom_point()  + 
    theme_classic2() +
    labs(title = NULL,
         y = "Counts",
         x = "Treatment",
         color = "")
  return(plot)
}

## heatmap plot 

ggheatmap <- function(countData,
                      colData,
                      treatment,
                      interac,
                      threshold,
                      alpha
) {
  if (interac == "None") {
    design <- as.formula(paste0("~ ", treatment))
  } else if (is.null(interac)) {
    design <- as.formula(paste0("~ ", treatment))
  } else {
    design =  as.formula(paste0("~ ", treatment, "+ ", interac))
  }
  deseq2Data <- DESeqDataSetFromMatrix(countData = as.matrix(countData),
                                       colData = colData,
                                       design = design)
  deseq2Data <- DESeq(deseq2Data)
  deseq2Results <- results(deseq2Data)
  deseq2ResDF <- as.data.frame(deseq2Results)
  
  deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 10, ]
  deseq2VST <- vst(deseq2Data)
  deseq2VST <- assay(deseq2VST)
  deseq2VST <- as.data.frame(deseq2VST)
  deseq2VST$Gene <- rownames(deseq2VST)
  head(deseq2VST)
  
  sigGenes <-
    rownames(deseq2ResDF[deseq2ResDF$padj <= as.numeric(alpha) &
                           abs(deseq2ResDF$log2FoldChange) > as.numeric(threshold), ])
  deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
  deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))
  
  
  heatmap <- ggplot(deseq2VST, 
                    aes(x=variable, y=Gene, fill=value)) +
    geom_raster() + scale_fill_viridis(trans="sqrt") + 
    theme(axis.text.x=element_text(angle=65, hjust=1),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(title = NULL,
         y = "Gene",
         x = NULL)
  
  deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
  rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
  deseq2VSTMatrix$Gene <- NULL
  distanceGene <- dist(deseq2VSTMatrix)
  distanceSample <- dist(t(deseq2VSTMatrix))
  clusterGene <- hclust(distanceGene, method="average")
  clusterSample <- hclust(distanceSample, method="average")
  sampleModel <- as.dendrogram(clusterSample)
  sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
  sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()
  deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])
  
  heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  heatmap
  grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))
  
  sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
  heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
  sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
  heatmapGrob <- ggplotGrob(heatmap_1)
  sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
  sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)
  maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
  sampleDendrogramGrob$widths <- as.list(maxWidth)
  heatmapGrob$widths <- as.list(maxWidth)
  
  heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
  finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))
  grid.draw(finalGrob)
  return(grid.draw(finalGrob))

}

# create funtions for DE with edgeR

get_model <- function(df, var) {
  if (length(var) == 1) {
    factor1 <- as.factor(df[[1]])
    model <- model.matrix(~ factor1)
    return(model)
  } else {
    factor1 <- as.factor(df[[1]])
    factor2 <- as.factor(df[[2]])
    model <- model.matrix(~ factor1 + factor2)
    return(model)
  }
}




DE_edgeR_design <- function(countData,
                            colData,
                            treatment,
                            threshold,
                            alpha
) {
  dgList <- DGEList(counts = countData, genes = rownames(countData))
  y <- dgList
  design <- get_model(colData, treatment)
  keep <- filterByExpr(y, design = design)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y) #
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  qlf$table$FDR <- p.adjust(qlf$table$PValue, method = "BH")
  qlf$table <- qlf$table %>%
    mutate(., Difference = (
      case_when(
        .data$logFC >= threshold & .data$FDR <= alpha ~ "UP",
        .data$logFC <= -threshold &
          .data$FDR <= alpha ~ "DOWN",
        .data$logFC <= threshold |
          .data$FDR > alpha ~ "Not significant",
        is.na(.data$FDR) ~ "Not significant"
      )
    ))
  
  return(qlf$table)
}

DE_edgeR_main <- function(countData,
                          colData,
                          interac = NULL,
                          treatment,
                          threshold,
                          alpha) {
  validate_row_cols(df_s =  colData,
                    df_r = countData)
  if (is.null(interac)) {
    
    treatment = c(treatment,interac)
    qlf <- DE_edgeR_design(countData = countData,
                           colData = colData,
                           treatment = treatment,
                           threshold = threshold,
                           alpha = alpha )
    
  } else if (interac == "None") {
    treatment =  treatment
    qlf <- DE_edgeR_design(countData = countData,
                           colData = colData,
                           treatment = treatment,
                           threshold = threshold,
                           alpha = alpha)
    
  } else {
    treatment =  c(treatment,interac)
    qlf <- DE_edgeR_design(countData = countData,
                           colData = colData,
                           treatment = treatment,
                           threshold = threshold,
                           alpha = alpha)
    
    
  }
  return(qlf)
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
  res <- results(dds, alpha = alpha, lfcThreshold = 0)
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
    mutate(., Difference = (
      case_when(
        .data$log2FoldChange >= threshold & .data$padj <= alpha ~ "UP",
        .data$log2FoldChange <= -threshold &
          .data$padj <= alpha ~ "DOWN",
        .data$log2FoldChange <= threshold |
          .data$padj > alpha ~ "Not significant",
        is.na(.data$padj) ~ "Not significant"
      )
    ))
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
  
  observeEvent(input$geneID,
               if (is.null(input$geneID)) {
                 updateSelectInput(inputId = "geneID",
                                   value = rownames(sampleinfo_data())[1])
               })



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
      validate("Treatment an interaction can't be the same")
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
  

  ## I have to conect the selectInput and add a conditional to choose between DEseq2 and edgeR
  
  


  results <- eventReactive(input$run,{
    withProgress(message = "Processing...", {
      if (input$treatment == input$interaction) {
        validate("Treatment an interaction can't be the same")
      }
      incProgress(1/8)
      if (input$pgk == "DESeq2") {
        incProgress(1/7)
        DE_DESeq2_main(
          countData = raw_counts_data(),
          colData = sampleinfo_data(),
          treatment = input$treatment,
          interac = input$interaction,
          alpha = input$pvalue,
          threshold = input$treshold
        )
      }else if (input$pgk == "edgeR") {
        incProgress(1/7)
        DE_edgeR_main(
          countData = raw_counts_data(),
          colData = sampleinfo_data(),
          treatment = input$treatment,
          interac = input$interaction,
          threshold = input$treshold,
          alpha = input$pvalue
        )
      }
      
    })

  })
  
  output$DE_results <- renderDataTable({
    results()
  })

  plot1 <- eventReactive(input$run,{

    if (input$treatment == input$interaction) {
      validate("Treatment an interaction can't be the same")
    }
    MDS_plot(
      countData = raw_counts_data(),
      colData = sampleinfo_data(),
      treatment =  input$treatment,
      interac = input$interaction
    )
  })
  
  plot2 <- eventReactive(input$run,{
    if (input$treatment == input$interaction) {
      validate("Treatment an interaction can't be the same")
    }
    MA_plot(results()) 
  })
  
  plot3 <- eventReactive(input$run,{
    if (input$treatment == input$interaction) {
      validate("Treatment an interaction can't be the same")
    }
    ggheatmap(countData = raw_counts_data(),
              colData = sampleinfo_data(),
              treatment = input$treatment,
              interac = input$interaction,
              threshold = input$treshold,
              alpha = input$pvalue
    )
  })
  
  plot4 <- eventReactive(input$geneID,{
    if (input$treatment == input$interaction) {
      validate("Treatment an interaction can't be the same")
    }
    if (is.null(input$geneID)) {
      validate("You have to provide a genID from raw counts data.")
      
    }
    if (input$geneID == "") {
      validate("You have to provide a genID from raw counts data.")
      
    }
    single_gen_plot(geneID = input$geneID,
                    countData = raw_counts_data(),
                    ColData = sampleinfo_data(),
                    treatment = input$treatment,
                    interaction = input$interaction)
  })
  
  output$plot1 <- renderPlot({
    plot1()
  })

  output$plot2 <- renderPlot({
    plot2()
  })
  
  output$plot3 <- renderPlot({
    plot3()
  }) 
  output$plot4 <- renderPlot({
    plot4()
  })
  
  output$downloadR <- downloadHandler(
    filename  = function() {
      paste0("Results_", input$pgk, ".csv")
    },
    content =  function(file) {
      write.csv(results(), file)
    }
  )


  output$downloadP1 <- downloadHandler(
    filename = function(){paste("plot_01",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plot1())
    }
  )
  output$downloadP2 <- downloadHandler(
    filename = function(){paste("plot_02",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plot2())
    }
  )
  output$downloadP3 <- downloadHandler(
    filename = function(){paste("plot_03",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plot3())
    }
  )
  output$downloadP4 <- downloadHandler(
    filename = function(){paste("plot_04",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plot4())
    }
  )
  
  
  output$readme <- renderUI({
    tags$div(includeHTML("~/Documentos/R/my_pkgs/DEasy/README.html"))
  })    
 
}) 
