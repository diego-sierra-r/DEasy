# Renam2 HTSEQ file
ara <- read.csv("~/Documentos/R/my_pkgs/DEasy/HTSeqcount_EDITED_mapeo_clean_ALL_AraBAN.txt",sep = "\t")
colnames(ara) <- colnames(ara) |> 
  stringr::str_replace_all("HTSeqcount_mapeo_cleanAraBAN",replacement = "Sample_")
write.csv(x = ara, "~/Documentos/R/my_pkgs/DEasy/ara_counts.csv",row.names = F)


