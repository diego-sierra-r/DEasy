library(DESeq2)
library(shiny)
library("pasilla")


data <- DT::datatable(iris)
data |>
  DT::formatStyle("Species",
                  backgroundColor = DT::styleEqual(unique(iris$Species),
                                                   values = c("gray","yellow","#76D7C4")
                                                   ))
set.seed(123)
colors <- c("#F08080","#E9967A","#DC143C",
                 "#FFC0CB","#FF69B4","#DB7093",
                 "#FFA07A","#FF6347","#FF8C00",
                 "#FFFF00","#FFFACD","#FFDAB9",
                 "#D8BFD8","#EE82EE","#BA55D3",
                 "#ADFF2F","#98FB98","#3CB371",
                 "#00FFFF","#E0FFFF","#7FFFD4",
                 "#48D1CC","#5F9EA0","#87CEFA	")
sample(list_colors,(2+2),replace = TRUE,)

rmarkdown::render("~/Documentos/R/my_pkgs/DEasy/README.Rmd",output_format = "html_document")
rmarkdown::render("~/Documentos/R/my_pkgs/DEasy/README.Rmd",output_format = "github_document")

# create skeleton of main DESeq2 function 

## import data form package
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)

rename_rows  <- function(df_s) {
  rownames(df_s) <- df_s[,1]
  df_s <- df_s[,2:ncol(df_s)]
  return(df_s)
}

validate_row_col <- function(df_s, df_r) {
  if (ncol(df_r) != nrow(df_s)) {
    checker <- colnames(df_r) %in% df_s[1]
    if (checker[1] == FALSE & checker[2] == TRUE) {
      df_s <- rename_rows(df_s = df_s)
      if (identical(rownames(df_s), colnames(df_r)) == TRUE) {
        return(df_s)
      } else {
        stop("Raw_counts columns and sample information colums must be excatly the same.")
      } 
    } else {
      stop("Raw_counts columns and sample information colums must be excatly the same...")
    } 

  } else {
    return(df_s)
  }
}

DE_DESeq2_design  <- function(countData,
                              colData,
                              treatment,
                              interac = NULL,
                              alpha,
                              threshold) {
  colData <- validate_row_col(df_s = colData, df_r = countData)
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
  colData <- validate_row_col(df_s = colData, df_r = countData)
  countData <- as.matrix(countData)
  if (all(rownames(colData) == colnames(countData)) == FALSE) {
    stop("Row names on Sample information and col names on raw count are not the same")
  }
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
  return(as.data.frame(res))
}


cts <- read.csv(pasCts,sep="\t",row.names="gene_id") |> as.matrix()
#cts <- cts[1:3000,] %>% as.matrix()
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
rownames(coldata) <- sub("fb", "", rownames(coldata)) ## coldata y count deben tner los minos nombre
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds,alpha = 0.05,lfcThreshold = 0)
res


ix <- DE_DESeq2_main(
  countData = cts,
  colData = coldata,
  treatment = "condition",
  alpha = 0.05,
  threshold = 0,
  interac = "type"
)

sample_information <- read.csv("tests/Samples_information.csv")
raw_coutns <- read.csv("tests/ara_counts.csv")

iix <- DE_DESeq2_main(countData = raw_coutns,
                      colData = sample_information,
                      treatment = "SEX",
                      interac = "None",
                      alpha = 0.05,
                      threshold = 0)

## MA-plot function
## MDS plot function
## Individual count plot funtion
## heatmapplot function