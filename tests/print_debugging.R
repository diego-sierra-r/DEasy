 ### Renam2 HTSEQ file
# ara <-
#   read.csv(
#     "~/Escritorio/Anexos objetivo 04/R_Transcriptome/Resultados DEG/DEGs_preliminareszip/HTSeqcount_EDITED_mapeo_clean_ALL_AraBAN.txt",
#     sep = "\t"
#   )
# colnames(ara) <- colnames(ara) |>
#   stringr::str_replace_all("HTSeqcount_mapeo_cleanAraBAN", replacement = "Sample_") |> 
#   stringr::str_replace_all(".bam.txt",replacement = "")
# write.csv(
#   x = ara,
#   "~/Documentos/R/my_pkgs/DEasy/tests/ara_counts.csv",
#   row.names = F
# )

# custom datatables

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
