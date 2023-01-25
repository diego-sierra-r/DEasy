# Renam2 HTSEQ file
ara <-
  read.csv(
    "~/Escritorio/Anexos objetivo 04/R_Transcriptome/Resultados DEG/DEGs_preliminareszip/HTSeqcount_EDITED_mapeo_clean_ALL_AraBAN.txt",
    sep = "\t"
  )
colnames(ara) <- colnames(ara) |>
  stringr::str_replace_all("HTSeqcount_mapeo_cleanAraBAN", replacement = "Sample_") |> 
  stringr::str_replace_all(".bam.txt",replacement = "")
write.csv(
  x = ara,
  "~/Documentos/R/my_pkgs/DEasy/tests/ara_counts.csv",
  row.names = F
)
