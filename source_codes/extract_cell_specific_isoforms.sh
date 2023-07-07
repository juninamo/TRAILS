#!/usr/bin/Rscript

# ROKU funtion 
library(TCC)                           
source("http://www.iu.a.u-tokyo.ac.jp/~kadota/R/R_functions.R")#Pre-load a file containing the kadota_2006_bmc_bioinformatics function to run #ROKU.
library(affy)                          #Loading package that contains tukey.biweight function to calculate Tukey's Biweight
library(som)                           #Loading package containing normalize function to normalize (mean=0, standard deviation=1) each gene expression vector
library(genefilter)
library(data.table)
library(tidyverse)

# custom code fo accelerating ROKU function
## https://sites.google.com/site/scriptofbioinformatics/maikuroarei-guan-xi/rokuno-gao-su-hua-r
AIC_0.25 <- function(x){
  if(length(x) == sum(is.na(x))){
    x <- c(rep(0, length(x)))
  }else if(length(x) == sum(is.nan(x))){
    x <- c(rep(0, length(x)))
  }else{}
  x_org <- x            
  x <- x[(!is.na(x))]   
  x <- x[(!is.nan(x))]  

  n_plus_s <- length(x)  
  x.sort <- sort(x)
  x.order <- order(x)

  limit <- 0.25*n_plus_s
  target <- rev(order(abs(x)))[1:limit]
  
  s <- 1:limit
  n <- n_plus_s - s
  set_sd <- c()
  for(i in 1:limit){set_sd[i] <- sd(x[setdiff(1:n_plus_s,target[1:i])])*sqrt((n[i]-1)/n[i])}
  Ut <- n*log(set_sd) + sqrt(2)*s*lgamma(n+1)/n
  maice_Ut <- 0
  maice_i <- 0
  maice_j <- 0
  if(maice_Ut > min(Ut)){
    maice_Ut <- which(Ut==min(Ut))
    maice_i <- length(which(x[target[1:maice_Ut]]<0))
    maice_j <- length(which(x[target[1:maice_Ut]]>0))
  }
  flag <- c(rep(0, length=n_plus_s))
  if(sd(x) != 0){
    if(maice_i > 0){ flag[x.order[1:(maice_i)]] <- -1 }
    if(maice_j > 0){ flag[x.order[n_plus_s - maice_j + 1 :n_plus_s]] <- 1 }
    tmp <- replace(x_org, ((!is.nan(x_org)) & (!is.na(x_org))), flag)
    return(tmp)
  }else{
    tmp <- replace(x_org, ((!is.nan(x_org)) & (!is.na(x_org))), flag)
    return(tmp)
  }
}

tra2symbol <- fread("~/reference/ENSEMBL/GRCh38/GTF/Homo_sapiens.GRCh38.101_gene_annotation_table.txt")
read2gene <- fread("/path/to/GRCh38_all-sr-correct_counts_matrix.tsv") %>%
  .$ids %>%
  stringr::str_split(., pattern = "_", simplify = TRUE) %>%
  as.data.frame() %>%
  dplyr::filter(grepl("ENSG", .$V2)) %>%
  dplyr::mutate(gene_id = stringr::str_split(V2, pattern = "\\.", simplify = TRUE)[,1]) %>%
  merge(., tra2symbol, by="gene_id") %>%
  dplyr::select(V1, GeneSymbol) %>%
  dplyr::rename(GENEID = GeneSymbol,query_name = V1)

## isoform count (RPM)
cts = fread("/path/to/GRCh38_all-sr-correct_counts_matrix.renamed.RPM.tsv") %>%
  tibble::column_to_rownames("V1") %>%
  as.matrix() %>%
  .[matrixStats::rowVars(.) != 0, ]
head(cts)

hoge <- ROKU(cts)                     
outlier <- hoge$outlier               
modH <- hoge$modH                      
ranking <- hoge$rank                  

data.z <- som::normalize(data.frame(cts), byrow=TRUE)  #Perform normalization and store results in data.z
out <- t(apply(data.z, 1, AIC_0.25))
colnames(out) <- colnames(data.frame(cts))        
entropy_score <- apply(data.frame(cts)+1, 1, kadota_2006_bmc_bioinformatics) #Transform the gene expression vector x one row (one gene) at a time (x' = |x - Tbw|>), then calculate entropy H(x') and store in entropy_score 

#Save to file
tmp <- cbind(rownames(data.frame(cts)), out, entropy_score)#The leftmost column is the gene ID, the next is the "outlier matrix information" consisting of columns for the number of samples, and the last column is the entropy value H(x'), which is stored in tmp as a matrix. 
write.table(tmp, "/path/to/specific_isoforms/entropy_score_isoform_count.txt", sep="\t", append=F, quote=F, row.names=F)

for (cell in colnames(cts)){
  print(cell)
  tmp <- rbind(out, ifelse(colnames(cts) == cell, 1, 0) %>% as.data.frame() %>%t())           #  Add the template pattern #template to the last row of the obtained outlier matrix out and store the result in tmp
  closeg <- genefinder(tmp, nrow(out)+1, nrow(out)) # Ranked by degree of specific expression and stored in closeg
  obj <- closeg[[1]]$indices[closeg[[1]]$dists == 0] # Store line number information in obj for genes that are exactly the same as the template
  print(length(obj))                            # Indicate how many specifically expressed genes of interest were found
  tmp <- cbind(rownames(data.frame(cts)), data.frame(cts), entropy_score)#Combine the entropy value H(x') on the right side of the input data and store the result in tmp
  tmp2 <- tmp[obj,]                      # Extract only the rows corresponding to the row number indicated by obj from the matrix tmp and store them in tmp2
  tmp3 <- tmp2[order(tmp2$entropy_score),] %>%
    magrittr::set_colnames(gsub("rownames\\(data.frame\\(cts)\\)","transcript_id",colnames(.))) #Sorting groups of specifically expressed genes of interest in order of decreasing entropy and storing them in tmp3
  write.table(tmp3, paste0("/path/to/specific_isoforms/specific_",cell,"_isoform_count.txt"), sep="\t", append=F, quote=F, row.names=F) 
}


## isoform usage
cts <- fread("/path/to/transcript_ratio.txt") %>%
  tibble::column_to_rownames("query_name") %>%
  .[,-1] %>%
  as.matrix() %>%
  .[matrixStats::rowVars(.) != 0, ]
head(cts)

hoge <- ROKU(cts)                     
outlier <- hoge$outlier               
modH <- hoge$modH                      
ranking <- hoge$rank                  

data.z <- som::normalize(data.frame(cts), byrow=TRUE)  #Perform normalization and store results in data.z
out <- t(apply(data.z, 1, AIC_0.25))
colnames(out) <- colnames(data.frame(cts))        
entropy_score <- apply(data.frame(cts)+1, 1, kadota_2006_bmc_bioinformatics) #Transform the gene expression vector x one row (one gene) at a time (x' = |x - Tbw|>), then calculate entropy H(x') and store in entropy_score 

#Save to file
tmp <- cbind(rownames(data.frame(cts)), out, entropy_score)#The leftmost column is the gene ID, the next is the "outlier matrix information" consisting of columns for the number of samples, and the last column is the entropy value H(x'), which is stored in tmp as a matrix. 
write.table(tmp, "/path/to/specific_isoforms/entropy_score_isoform_usage.txt", sep="\t", append=F, quote=F, row.names=F)

for (cell in colnames(cts)){
  print(cell)
  tmp <- rbind(out, ifelse(colnames(cts) == cell, 1, 0) %>% as.data.frame() %>%t())           #  Add the template pattern #template to the last row of the obtained outlier matrix out and store the result in tmp
  closeg <- genefinder(tmp, nrow(out)+1, nrow(out)) # Ranked by degree of specific expression and stored in closeg
  obj <- closeg[[1]]$indices[closeg[[1]]$dists == 0] # Store line number information in obj for genes that are exactly the same as the template
  print(length(obj))                            # Indicate how many specifically expressed genes of interest were found
  tmp <- cbind(rownames(data.frame(cts)), data.frame(cts), entropy_score)#Combine the entropy value H(x') on the right side of the input data and store the result in tmp
  tmp2 <- tmp[obj,]                      # Extract only the rows corresponding to the row number indicated by obj from the matrix tmp and store them in tmp2
  tmp3 <- tmp2[order(tmp2$entropy_score),] %>%
    magrittr::set_colnames(gsub("rownames\\(data.frame\\(cts)\\)","transcript_id",colnames(.))) #Sorting groups of specifically expressed genes of interest in order of decreasing entropy and storing them in tmp3
  write.table(tmp3, paste0("/path/to/specific_isoforms/specific_",cell,"_isoform_usage.txt"), sep="\t", append=F, quote=F, row.names=F) #tmp3の中身を指定したファイル名で保存
}