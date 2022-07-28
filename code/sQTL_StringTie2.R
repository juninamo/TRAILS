#!/usr/bin/R

suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(BEDMatrix))
suppressMessages(library(dplyr))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(gassocplot2))
suppressMessages(library(coloc))

setDTthreads(1)
getDTthreads()


read2gene <- fread("/path/to/GRCh38_all-sr-correct_counts_matrix.tsv") %>%
  .$ids %>%
  stringr::str_split(., pattern = "_", simplify = TRUE) %>%
  as.data.frame() %>%
  dplyr::filter(grepl("ENSG", .$V2)) %>%
  dplyr::mutate(gene_id = stringr::str_split(V2, pattern = "\\.", simplify = TRUE)[,1]) %>%
  dplyr::inner_join(., fread("/path/to/Homo_sapiens.GRCh38.101_gene_annotation_table.txt"), by="gene_id") %>%
  dplyr::select(V1, GeneSymbol) %>%
  dplyr::rename(GENEID = GeneSymbol,query_name = V1)



## EGA
IID_tr <- read.table("/path/to/EGA/sample_Info.txt",header = TRUE)

# covariate file
SUBSET=1:5
POP = "EUR"
for (i in 1:length(SUBSET)){
  for (p in 1:length(POP)){
    tryCatch({
      
      cov = as.data.frame(fread(paste("/path/to/EGA/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",SUBSET[i],"_",POP[p],"_PC10.txt",sep=""),header = TRUE))
      cov$SampleID <- str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
      colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
      assign(paste0("cov_",SUBSET[i],"_",POP[p]), cov)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}



# ImmVar
IID_tr_i <- read.table("/path/to/ImmVar/sample_Info.txt",header = TRUE)

# covariate file
SUBSET <- c("CD4T_Activated","MoDC_unstim","MoDC_FLU","MoDC_IFNb")
POP = "EUR"
for (i in 1:length(SUBSET)){
  for (p in 1:length(POP)){
    tryCatch({
      
      cov = as.data.frame(fread(paste("/path/to/ImmVar/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",SUBSET[i],"_",POP[p],".hg38_PC10.txt",sep=""),header = TRUE))
      cov$SampleID <- str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
      colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
      assign(paste0("cov_",SUBSET[i],"_",POP[p]), cov)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}



# DICE
corres = read.table("/path/to/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE)

# covariate file
cov = as.data.frame(fread(paste("/path/to/DICE/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup.hg38_PC10.txt",sep=""),header = TRUE)) %>%
  dplyr::mutate(SampleID = str_split(SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(SampleID, pattern = "_", simplify = TRUE))]) %>%
  magrittr::set_colnames(gsub("SampleID","id",colnames(.)%>%gsub("^0_","",.))) %>%
  tibble::column_to_rownames("id") %>%
  t() %>%
  merge(.,read.table("/path/to/DICE/DICE_DNA_RNAseq_correspondance.txt", header = TRUE),by.x="row.names",by.y="IID") %>%
  tibble::column_to_rownames("Run") %>%
  dplyr::select(starts_with("PC")) %>%
  t() %>% as.data.frame() %>%
  dplyr::mutate(id = rownames(.)) %>%
  .[,c("id",read.table("/path/to/DICE/StringTie/sample_list.txt", header = TRUE)[,1])]
assign(paste0("cov_DICE"), cov)



# GEUV
IID_tr_g <- 
  read.table("/path/to/GEUV/all_sample.txt", sep="\t") %>%
  .[,c("V1","V3")] %>%
  magrittr::set_colnames(c("IID","stim")) %>%
  magrittr::set_rownames(.$IID) %>%
  .[read.table("/path/to/GEUV/StringTie/sample_list.txt", header = TRUE)[,1],] %>%
  dplyr::mutate(stim = ifelse(grepl("Yoruba", stim), "YRI", "EUR")) %>%
  .[,c("stim","IID")]


# covariate file
POP = "EUR"
for (p in 1:length(POP)){
  tryCatch({
    
    cov = as.data.frame(fread(paste("/path/to/GEUV/VCF/GEUVADIS_GRCh38_",POP[p],"_PC10.txt",sep=""),header = TRUE))
    cov$SampleID <- str_split(cov$SampleID, pattern = "_", simplify = TRUE)[,ncol(str_split(cov$SampleID, pattern = "_", simplify = TRUE))]
    colnames(cov) <- gsub("SampleID","id",colnames(cov)%>%gsub("^0_","",.))
    assign(paste0("cov_GEUV_",POP[p]), cov)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# normalized expression
isopos_c = fread("/path/to/eQTL/StringTie/geneInfo.txt") %>% 
  merge(read2gene,.,by.x="query_name",by.y="transcript_id") %>%
  dplyr::mutate(chr = paste0("chr",chr)) %>%
  as.data.frame()
psi_c <- fread(paste0("/path/to/sQTL/StringTie/quantile-rank-normalized_transcript_ratio.txt")) %>%
  dplyr::rename(gid = V1) %>%
  as.data.frame()
peer_factors <- fread(paste0("/path/to/sQTL/StringTie/PEERfactors.txt")) %>%
  as.data.frame()


chr = commandArgs(trailingOnly=TRUE)[1] %>% as.integer()
# chr=1

pop = "EUR"

# DICE

dir.create(paste0("/path/to/DICE/QTL/result/splicing_peer"), showWarnings = F, recursive = T)

for (stim in c("TFH","TH1","TH17","TH2","TH1-17","TREG_MEMORY","TREG_NAIVE","B_NAIVE","CD4_NAIVE","CD4_N_STIM","CD8_NAIVE","CD8_N_STIM","NONCLASSICAL_MONOCYTES","CLASSICAL_MONOCYTES","NK_CD16POS")){
  tryCatch({
    
    print(stim)
    variant = fread(paste("/path/to/DICE/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bim",sep="")) %>%
      dplyr::filter(V1 == chr) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/path/to/DICE/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.)))
    
    # Position files
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    
    cov <- dplyr::bind_rows(eval(parse(text=paste0("cov_DICE"))),
                            peer_factors %>%
                              dplyr::rename(id = PEERfactors) %>%
                              .[,grepl(paste0("^id$|^SRR"),colnames(.))])
    cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
      magrittr::set_colnames(colnames(cov))
    
    shared_samples = 
      intersect(
        corres[ grepl("_1_RNA_",corres$biospecimen_repository_sample_id) & corres$histological_type == stim , "Run"],
        colnames(psi_c) %>% .[grepl("^SRR",.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1]
      )  %>%
      intersect(
        .,
        colnames(cov)
      )
    
    shared_IID = 
      corres %>%
      dplyr::filter(Run %in% shared_samples) %>%
      .$IID %>% as.character() %>%
      intersect(
        .,
        colnames(snps_i)
      )
    shared = 
      corres %>%
      dplyr::filter(Run %in% shared_samples) %>%
      dplyr::filter(IID %in% shared_IID) %>%
      magrittr::set_rownames(.$IID) %>%
      .[shared_IID,]
    
    snps = snps_i[,c("id",shared$IID %>% as.character())] %>%
      magrittr::set_colnames(c("id",shared$Run %>% as.character()))
    expr2 <- psi_c[, c("gid", shared$Run %>% as.character()) ] %>%
      .[ which ( .$gid %in% isopos_c[isopos_c$chr == paste0("chr",chr) & 
                                       isopos_c$start > min(snpspos$pos-(5e+05)) & 
                                       isopos_c$end < max(snpspos$pos+(5e+05)), "query_name"] ), ]
    cov <- cov[,c("id",shared$Run %>% as.character())]
    
    
    if ( nrow(expr2) > 0 &&
         all(colnames(snps)[-1] == colnames(expr2)[-1] &&
             colnames(snps)[-1] == colnames(cov)[-1]) ) {
      
      # print(head(expr2))
      # print(dim(expr2))
      # print(head(snps))
      # print(dim(snps))
      # print(head(cov))
      # print(dim(cov))
      
      # paste("Mean value of expression for gene ",expr2$gid," is ", rowMeans(expr2[, -1]))
      # paste("Standard deviation of expression for gene ", expr2$gid," is ", t(apply(expr2[, -1], 1, sd)))
      # paste("Mean value for SNP ",snps$id," is ", rowMeans(snps[, -1]))
      # paste("Standard deviation for SNP ", snps$id," is ", t(apply(snps[, -1], 1, sd)))
      
      # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
      write.table(snps, paste0("/path/to/snps_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(expr2, paste0("/path/to/exp_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(cov, paste0("/path/to/cov_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      
      # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      
      # Genotype file name
      SNP_file_name = paste("/path/to/snps_chr",chr,"_str.txt", sep="");
      
      # Gene expression file name
      expression_file_name = paste("/path/to/exp_chr",chr,"_str.txt", sep="");
      
      # Covariates file name
      # Set to character() for no covariates
      covariates_file_name = paste("/path/to/cov_chr",chr,"_str.txt", sep="");
      
      
      # Error covariance matrix
      # Set to numeric() for identity.
      errorCovariance = numeric();
      
      pvOutputThreshold_cis = 1;
      pvOutputThreshold_tra = 0;
      
      # Distance for local gene-SNP pairs
      cisDist = 5e5;
      
      
      # Output file name
      output_file_name_cis = tempfile();
      output_file_name_tra = tempfile();
      
      ## Load genotype data
      SNPs = SlicedData$new();
      SNPs$fileDelimiter = "\t";      # the TAB character
      SNPs$fileOmitCharacters = "NA"; # denote missing values;
      SNPs$fileSkipRows = 1;          # one row of column labels
      SNPs$fileSkipColumns = 1;       # one column of row labels
      SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      SNPs$LoadFile(SNP_file_name);
      
      ## Load gene expression data
      gene = SlicedData$new();
      gene$fileDelimiter = "\t";      # the TAB character
      gene$fileOmitCharacters = "NA"; # denote missing values;
      gene$fileSkipRows = 1;          # one row of column labels
      gene$fileSkipColumns = 1;       # one column of row labels
      gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      gene$LoadFile(expression_file_name);
      
      ## Load covariates
      cvrt = SlicedData$new();
      cvrt$fileDelimiter = "\t";      # the TAB character
      cvrt$fileOmitCharacters = "NA"; # denote missing values;
      cvrt$fileSkipRows = 1;          # one row of column labels
      cvrt$fileSkipColumns = 1;       # one column of row labels
      if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
      }
      
      
      ## Run the analysis
      
      me = Matrix_eQTL_main(
        snps = SNPs,
        gene = gene,
        cvrt = cvrt,
        output_file_name     = output_file_name_tra,
        pvOutputThreshold     = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = isopos_c[,c(1,3,4,5)],
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
      
      unlink(output_file_name_cis);
      
      
      ## Results:
      cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
      # cat('Detected local eQTLs:', '\n');
      # show(me$cis$eqtls)
      # cat('Detected distant eQTLs:', '\n');
      # show(me$trans$eqtls)
      
      ## Results:
      
      me$cis$eqtls %>%
        dplyr::inner_join(.,
                          fread(paste("/path/to/DICE/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_chr",chr,".hg38.bim",sep="")) %>%
                            dplyr::filter(V1 == chr) %>%
                            dplyr::select(V2,V5,V6) %>%
                            magrittr::set_colnames(c("snps","A1","A0")),
                          by="snps") %>%
        merge(.,
              isopos_c,
              by.x="gene",by.y="query_name") %>%
        dplyr::arrange(pvalue) %>%
        write.table(., 
                    paste0("/path/to/DICE/QTL/result/splicing_peer/chr",chr,"_in_",stim,"_StringTie.txt"),
                    sep = "\t",quote=FALSE, row.names = FALSE)
      
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# EGA

dir.create(paste0("/path/to/EGA/QTL/result/splicing_peer"), showWarnings = F, recursive = T)

for (stim in 1:5){
  tryCatch({
    
    print(stim)
    variant = fread(paste("/path/to/EGA/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stim,"_",pop,"_chr",chr,".hg38.bim",sep="")) %>%
      dplyr::filter(V1 == chr) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/path/to/EGA/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stim,"_",pop,"_chr",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","@",colnames(.)))
    
    # Position files
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    
    stimulus = c("Mono_NS","Mono_LPS","Mono_Pam3CSK4","Mono_R848","Mono_IAV")[stim]
    
    cov <- dplyr::bind_rows(eval(parse(text=paste0("cov_",stim,"_",pop))),
                            peer_factors %>%
                              dplyr::rename(id = PEERfactors) %>%
                              .[,grepl(paste0("^id$|",stimulus),colnames(.))] %>%
                              magrittr::set_colnames(gsub(paste0(";",stimulus),"",colnames(.))))
    cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
      magrittr::set_colnames(colnames(cov))
    
    shared_samples = 
      intersect(
        colnames(psi_c) %>% .[grepl(stimulus,.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
        colnames(snps_i)
      ) %>%
      intersect(
        .,
        colnames(cov)
      )
    
    snps = snps_i[,c("id",shared_samples)]
    expr2 <- psi_c[, c("gid", paste0(shared_samples,";",stimulus)) ] %>%
      magrittr::set_colnames(gsub(paste0(";",stimulus),"",colnames(.))) %>%
      .[ which ( .$gid %in% isopos_c[isopos_c$chr == paste0("chr",chr) & 
                                       isopos_c$start > min(snpspos$pos-(5e+05)) & 
                                       isopos_c$end < max(snpspos$pos+(5e+05)), "query_name"] ), ]
    cov <- cov[,c("id",shared_samples)]
    
    
    if ( nrow(expr2) > 0 &&
         all(colnames(snps)[-1] == colnames(expr2)[-1] &&
             colnames(snps)[-1] == colnames(cov)[-1]) ) {
      
      # print(head(expr2))
      # print(dim(expr2))
      # print(head(snps))
      # print(dim(snps))
      # print(head(cov))
      # print(dim(cov))
      
      # paste("Mean value of expression for gene ",expr2$gid," is ", rowMeans(expr2[, -1]))
      # paste("Standard deviation of expression for gene ", expr2$gid," is ", t(apply(expr2[, -1], 1, sd)))
      # paste("Mean value for SNP ",snps$id," is ", rowMeans(snps[, -1]))
      # paste("Standard deviation for SNP ", snps$id," is ", t(apply(snps[, -1], 1, sd)))
      
      # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
      write.table(snps, paste0("/path/to/snps_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(expr2, paste0("/path/to/exp_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(cov, paste0("/path/to/cov_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      
      # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      
      # Genotype file name
      SNP_file_name = paste("/path/to/snps_chr",chr,"_str.txt", sep="");
      
      # Gene expression file name
      expression_file_name = paste("/path/to/exp_chr",chr,"_str.txt", sep="");
      
      # Covariates file name
      # Set to character() for no covariates
      covariates_file_name = paste("/path/to/cov_chr",chr,"_str.txt", sep="");
      
      
      # Error covariance matrix
      # Set to numeric() for identity.
      errorCovariance = numeric();
      
      pvOutputThreshold_cis = 1;
      pvOutputThreshold_tra = 0;
      
      # Distance for local gene-SNP pairs
      cisDist = 5e5;
      
      
      # Output file name
      output_file_name_cis = tempfile();
      output_file_name_tra = tempfile();
      
      ## Load genotype data
      SNPs = SlicedData$new();
      SNPs$fileDelimiter = "\t";      # the TAB character
      SNPs$fileOmitCharacters = "NA"; # denote missing values;
      SNPs$fileSkipRows = 1;          # one row of column labels
      SNPs$fileSkipColumns = 1;       # one column of row labels
      SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      SNPs$LoadFile(SNP_file_name);
      
      ## Load gene expression data
      gene = SlicedData$new();
      gene$fileDelimiter = "\t";      # the TAB character
      gene$fileOmitCharacters = "NA"; # denote missing values;
      gene$fileSkipRows = 1;          # one row of column labels
      gene$fileSkipColumns = 1;       # one column of row labels
      gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      gene$LoadFile(expression_file_name);
      
      ## Load covariates
      cvrt = SlicedData$new();
      cvrt$fileDelimiter = "\t";      # the TAB character
      cvrt$fileOmitCharacters = "NA"; # denote missing values;
      cvrt$fileSkipRows = 1;          # one row of column labels
      cvrt$fileSkipColumns = 1;       # one column of row labels
      if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
      }
      
      
      ## Run the analysis
      
      me = Matrix_eQTL_main(
        snps = SNPs,
        gene = gene,
        cvrt = cvrt,
        output_file_name     = output_file_name_tra,
        pvOutputThreshold     = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = isopos_c[,c(1,3,4,5)],
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
      
      unlink(output_file_name_cis);
      
      
      ## Results:
      cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
      # cat('Detected local eQTLs:', '\n');
      # show(me$cis$eqtls)
      # cat('Detected distant eQTLs:', '\n');
      # show(me$trans$eqtls)
      
      ## Results:
      me$cis$eqtls %>%
        dplyr::inner_join(.,
                          fread(paste("/path/to/EGA/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-",stim,"_",pop,"_chr",chr,".hg38.bim",sep="")) %>%
                            dplyr::filter(V1 == chr) %>%
                            dplyr::select(V2,V5,V6) %>%
                            magrittr::set_colnames(c("snps","A1","A0")),
                          by="snps") %>%
        merge(.,
              isopos_c,
              by.x="gene",by.y="query_name") %>%
        dplyr::arrange(pvalue) %>%
        write.table(., 
                    paste0("/path/to/EGA/QTL/result/splicing_peer/chr",chr,"_in_",stimulus,"_StringTie.txt"),
                    sep = "\t",quote=FALSE, row.names = FALSE)
      
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# ImmVar
dir.create(paste0("/path/to/ImmVar/QTL/result/splicing_peer"), showWarnings = F, recursive = T)

for (stim in c("CD4T_Activated","MoDC_unstim","MoDC_FLU","MoDC_IFNb")){
  tryCatch({
    
    print(stim)
    variant = fread(paste("/path/to/ImmVar/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",stim,"_",pop,"_chr",chr,".hg38.bim",sep="")) %>%
      dplyr::filter(V1 == chr) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/path/to/ImmVar/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",stim,"_",pop,"_chr",chr,".hg38.bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.)))
    
    # Position files
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    
    
    cov <- dplyr::bind_rows(eval(parse(text=paste0("cov_",stim,"_",pop))),
                            peer_factors %>%
                              dplyr::rename(id = PEERfactors) %>%
                              .[,grepl(paste0("^id$|",stim),colnames(.))] %>%
                              magrittr::set_colnames(gsub(paste0(";",stim),"",colnames(.))))
    cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
      magrittr::set_colnames(colnames(cov))
    
    shared_samples = 
      intersect(
        colnames(psi_c) %>% .[grepl(stim,.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1],
        colnames(snps_i)
      ) %>%
      intersect(
        .,
        colnames(cov)
      )
    
    snps = snps_i[,c("id",shared_samples)]
    expr2 <- psi_c[, c("gid", paste0(shared_samples,";",stim)) ] %>%
      magrittr::set_colnames(gsub(paste0(";",stim),"",colnames(.))) %>%
      .[ which ( .$gid %in% isopos_c[isopos_c$chr == paste0("chr",chr) & 
                                       isopos_c$start > min(snpspos$pos-(5e+05)) & 
                                       isopos_c$end < max(snpspos$pos+(5e+05)), "query_name"] ), ]
    cov <- cov[,c("id",shared_samples)]
    
    
    if ( nrow(expr2) > 0 &&
         all(colnames(snps)[-1] == colnames(expr2)[-1] &&
             colnames(snps)[-1] == colnames(cov)[-1]) ) {
      
      # print(head(expr2))
      # print(dim(expr2))
      # print(head(snps))
      # print(dim(snps))
      # print(head(cov))
      # print(dim(cov))
      
      # paste("Mean value of expression for gene ",expr2$gid," is ", rowMeans(expr2[, -1]))
      # paste("Standard deviation of expression for gene ", expr2$gid," is ", t(apply(expr2[, -1], 1, sd)))
      # paste("Mean value for SNP ",snps$id," is ", rowMeans(snps[, -1]))
      # paste("Standard deviation for SNP ", snps$id," is ", t(apply(snps[, -1], 1, sd)))
      
      # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
      write.table(snps, paste0("/path/to/snps_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(expr2, paste0("/path/to/exp_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(cov, paste0("/path/to/cov_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      
      # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      
      # Genotype file name
      SNP_file_name = paste("/path/to/snps_chr",chr,"_str.txt", sep="");
      
      # Gene expression file name
      expression_file_name = paste("/path/to/exp_chr",chr,"_str.txt", sep="");
      
      # Covariates file name
      # Set to character() for no covariates
      covariates_file_name = paste("/path/to/cov_chr",chr,"_str.txt", sep="");
      
      
      # Error covariance matrix
      # Set to numeric() for identity.
      errorCovariance = numeric();
      
      pvOutputThreshold_cis = 1;
      pvOutputThreshold_tra = 0;
      
      # Distance for local gene-SNP pairs
      cisDist = 5e5;
      
      
      # Output file name
      output_file_name_cis = tempfile();
      output_file_name_tra = tempfile();
      
      ## Load genotype data
      SNPs = SlicedData$new();
      SNPs$fileDelimiter = "\t";      # the TAB character
      SNPs$fileOmitCharacters = "NA"; # denote missing values;
      SNPs$fileSkipRows = 1;          # one row of column labels
      SNPs$fileSkipColumns = 1;       # one column of row labels
      SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      SNPs$LoadFile(SNP_file_name);
      
      ## Load gene expression data
      gene = SlicedData$new();
      gene$fileDelimiter = "\t";      # the TAB character
      gene$fileOmitCharacters = "NA"; # denote missing values;
      gene$fileSkipRows = 1;          # one row of column labels
      gene$fileSkipColumns = 1;       # one column of row labels
      gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      gene$LoadFile(expression_file_name);
      
      ## Load covariates
      cvrt = SlicedData$new();
      cvrt$fileDelimiter = "\t";      # the TAB character
      cvrt$fileOmitCharacters = "NA"; # denote missing values;
      cvrt$fileSkipRows = 1;          # one row of column labels
      cvrt$fileSkipColumns = 1;       # one column of row labels
      if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
      }
      
      
      ## Run the analysis
      
      me = Matrix_eQTL_main(
        snps = SNPs,
        gene = gene,
        cvrt = cvrt,
        output_file_name     = output_file_name_tra,
        pvOutputThreshold     = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = isopos_c[,c(1,3,4,5)],
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
      
      unlink(output_file_name_cis);
      
      
      ## Results:
      cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
      # cat('Detected local eQTLs:', '\n');
      # show(me$cis$eqtls)
      # cat('Detected distant eQTLs:', '\n');
      # show(me$trans$eqtls)
      
      ## Results:
      me$cis$eqtls %>%
        dplyr::inner_join(.,
                          fread(paste("/path/to/ImmVar/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_",stim,"_",pop,"_chr",chr,".hg38.bim",sep="")) %>%
                            dplyr::filter(V1 == chr) %>%
                            dplyr::select(V2,V5,V6) %>%
                            magrittr::set_colnames(c("snps","A1","A0")),
                          by="snps") %>%
        merge(.,
              isopos_c,
              by.x="gene",by.y="query_name") %>%
        dplyr::arrange(pvalue) %>%
        write.table(., 
                    paste0("/path/to/ImmVar/QTL/result/splicing_peer/chr",chr,"_in_",stim,"_StringTie.txt"),
                    sep = "\t",quote=FALSE, row.names = FALSE)
      
      
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# GEUV

dir.create(paste0("/path/to/GEUV/QTL/result/splicing_peer"), showWarnings = F, recursive = T)

for (stim in c("EUR")){
  tryCatch({
    
    print("LCL")
    variant = fread(paste("/path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_",pop,"_chr",chr,".bim",sep="")) %>%
      dplyr::filter(V1 == chr) %>%
      as.data.frame() %>%
      .[,"V2"]
    snps_i = BEDMatrix(paste("/path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_",pop,"_chr",chr,".bed",sep="")) %>%
      as.matrix() %>%
      magrittr::set_colnames(paste(stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,1],
                                   stringr::str_split(colnames(.), pattern = "_", simplify = TRUE) %>% .[,2],
                                   sep="_")) %>%
      .[, which (colnames(.) %in% variant)] %>%
      as.matrix() %>%
      t()%>%
      data.frame(id = variant, .) %>%
      magrittr::set_colnames(gsub("^X0.","",colnames(.))) %>%
      magrittr::set_colnames(gsub("\\.","-",colnames(.)))
    
    # Position files
    snpspos = data.frame(snp = snps_i$id,
                         chr = paste0("chr",str_split(snps_i$id, pattern = "_", simplify = TRUE)[,1]),
                         pos = as.integer(str_split(snps_i$id, pattern = "_", simplify = TRUE)[,2]))
    
    
    cov <- dplyr::bind_rows(eval(parse(text=paste0("cov_GEUV_EUR"))),
                            peer_factors %>%
                              dplyr::rename(id = PEERfactors) %>%
                              .[,grepl(paste0("^id$|^ERR"),colnames(.))])
    cov <-  data.frame(id = cov$id, impute::impute.knn(as.matrix(cov[,-1]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data ) %>%
      magrittr::set_colnames(colnames(cov))
    
    shared_samples = 
      intersect(
        IID_tr_g[IID_tr_g$stim == "EUR", "IID"],
        colnames(psi_c) %>% .[grepl("^ERR",.)] %>% stringr::str_split(., pattern = ";", simplify = TRUE) %>% .[,1]
      )  %>%
      intersect(
        .,
        colnames(cov)
      )
    
    snps = snps_i[, c("id",shared_samples) ]
    expr2 <- psi_c[, c("gid", shared_samples) ] %>%
      .[ which ( .$gid %in% isopos_c[isopos_c$chr == paste0("chr",chr) & 
                                       isopos_c$start > min(snpspos$pos-(5e+05)) & 
                                       isopos_c$end < max(snpspos$pos+(5e+05)), "query_name"] ), ]
    cov <- cov[,c( "id",shared_samples )]
    
    
    if ( nrow(expr2) > 0 &&
         all(colnames(snps)[-1] == colnames(expr2)[-1] &&
             colnames(snps)[-1] == colnames(cov)[-1]) ) {
      
      # print(head(expr2))
      # print(dim(expr2))
      # print(head(snps))
      # print(dim(snps))
      # print(head(cov))
      # print(dim(cov))
      
      # paste("Mean value of expression for gene ",expr2$gid," is ", rowMeans(expr2[, -1]))
      # paste("Standard deviation of expression for gene ", expr2$gid," is ", t(apply(expr2[, -1], 1, sd)))
      # paste("Mean value for SNP ",snps$id," is ", rowMeans(snps[, -1]))
      # paste("Standard deviation for SNP ", snps$id," is ", t(apply(snps[, -1], 1, sd)))
      
      # MatrixQTLでmQTL-SNPsと発現量が相関する遺伝子を抽出
      write.table(snps, paste0("/path/to/snps_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(expr2, paste0("/path/to/exp_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      write.table(cov, paste0("/path/to/cov_chr",chr,"_str.txt"), sep = "\t",quote=FALSE, row.names = FALSE)
      
      # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
      
      # Genotype file name
      SNP_file_name = paste("/path/to/snps_chr",chr,"_str.txt", sep="");
      
      # Gene expression file name
      expression_file_name = paste("/path/to/exp_chr",chr,"_str.txt", sep="");
      
      # Covariates file name
      # Set to character() for no covariates
      covariates_file_name = paste("/path/to/cov_chr",chr,"_str.txt", sep="");
      
      
      # Error covariance matrix
      # Set to numeric() for identity.
      errorCovariance = numeric();
      
      pvOutputThreshold_cis = 1;
      pvOutputThreshold_tra = 0;
      
      # Distance for local gene-SNP pairs
      cisDist = 5e5;
      
      
      # Output file name
      output_file_name_cis = tempfile();
      output_file_name_tra = tempfile();
      
      ## Load genotype data
      SNPs = SlicedData$new();
      SNPs$fileDelimiter = "\t";      # the TAB character
      SNPs$fileOmitCharacters = "NA"; # denote missing values;
      SNPs$fileSkipRows = 1;          # one row of column labels
      SNPs$fileSkipColumns = 1;       # one column of row labels
      SNPs$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      SNPs$LoadFile(SNP_file_name);
      
      ## Load gene expression data
      gene = SlicedData$new();
      gene$fileDelimiter = "\t";      # the TAB character
      gene$fileOmitCharacters = "NA"; # denote missing values;
      gene$fileSkipRows = 1;          # one row of column labels
      gene$fileSkipColumns = 1;       # one column of row labels
      gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
      gene$LoadFile(expression_file_name);
      
      ## Load covariates
      cvrt = SlicedData$new();
      cvrt$fileDelimiter = "\t";      # the TAB character
      cvrt$fileOmitCharacters = "NA"; # denote missing values;
      cvrt$fileSkipRows = 1;          # one row of column labels
      cvrt$fileSkipColumns = 1;       # one column of row labels
      if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
      }
      
      
      ## Run the analysis
      
      me = Matrix_eQTL_main(
        snps = SNPs,
        gene = gene,
        cvrt = cvrt,
        output_file_name     = output_file_name_tra,
        pvOutputThreshold     = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = isopos_c[,c(1,3,4,5)],
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
      
      unlink(output_file_name_cis);
      
      
      ## Results:
      cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
      # cat('Detected local eQTLs:', '\n');
      # show(me$cis$eqtls)
      # cat('Detected distant eQTLs:', '\n');
      # show(me$trans$eqtls)
      
      ## Results:
      me$cis$eqtls %>%
        dplyr::inner_join(.,
                          fread(paste("/path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_",pop,"_chr",chr,".bim",sep="")) %>%
                            dplyr::filter(V1 == chr) %>%
                            dplyr::select(V2,V5,V6) %>%
                            magrittr::set_colnames(c("snps","A1","A0")),
                          by="snps") %>%
        merge(.,
              isopos_c,
              by.x="gene",by.y="query_name") %>%
        dplyr::arrange(pvalue) %>%
        write.table(., 
                    paste0("/path/to/GEUV/QTL/result/splicing_peer/chr",chr,"_in_",stim,"_StringTie.txt"),
                    sep = "\t",quote=FALSE, row.names = FALSE)
      
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
