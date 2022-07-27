#!/bin/sh
#$ -S /bin/sh

# cluster isoforms by identical coding sequence
vsearch \
--cluster_fast /path/to/GRCh38_gencode38_classification.filtered_lite.cds.fasta \
--id 1 \
--centroids /path/to/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.centroids.fasta \
--clusterout_id \
--consout /path/to/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.consensus.fasta \
--msaout /path/to/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.msa.fasta \
--otutabout /path/to/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.table

arg = paste0("seqkit seq -n /path/to/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.msa.fasta > tmp.txt") 
system(arg , ignore.stdout=F,ignore.stderr=F)
cluster <- fread("tmp.txt", header = FALSE, fill = TRUE, sep = "\t") %>%
  dplyr::mutate(isoform = stringr::str_split(V1,pattern="_",simplify=TRUE) %>% .[,1],
                group = 1) %>%
  dplyr::select(-V1)
group_n = 1
for(i in 1:nrow(cluster)){
  print(i)
  if(cluster$isoform[i] == "consensus"){
    cluster$group[i] <- group_n
    group_n = group_n + 1
  } else {
    cluster$group[i] <- group_n
  }
}
centroids = cluster %>%
  dplyr::filter(grepl("^\\*",isoform))
cluster %>%
  dplyr::filter(isoform != "consensus") %>%
  merge(.,centroids,by="group") %>%
  dplyr::select(-group) %>%
  magrittr::set_colnames(c("isoform","centroid")) %>%
  dplyr::mutate(isoform = gsub("^\\*","",isoform),
                centroid = gsub("^\\*","",centroid)) %>%
  write.table(., paste0("/path/to/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.cluster_table.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
  
