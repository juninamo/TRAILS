fl = list.files(path = paste0("../data/ATAC/"), pattern = paste0(".bdg.gz$"), full.names = T)
cell_name = stringr::str_split(fl,pattern="/",simplify = TRUE) %>% .[,5] %>% gsub("_merged_treat_pileup.bdg.gz","",.) %>% gsub("_no_treament"," [unstim]", .) %>% gsub("_treament"," [stim]", .) %>% gsub("_"," ",.)
import_ATAC <- function(path, region, verbose = FALSE){
  
  if(verbose){
    print(path)
  }
  
  #Fetch summary statistics with seqminer
  summary_stats = seqminer::tabix.read.table(tabixFile = path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(isoform = stringr::str_split(path,pattern="/",simplify = TRUE) %>% .[,5] %>% gsub("_merged_treat_pileup.bdg.gz","",.) %>% gsub("_no_treament"," [unstim]", .) %>% gsub("_treament"," [stim]", .) %>% gsub("_"," ",.),
                  type = "ATAC-seq")
  colnames(summary_stats) <- c("chr","start","end","height","isoform","type")
  return(summary_stats)
  
}
safe_import_ATAC = purrr::safely(import_ATAC)
# sum read for normalizing ATAC-seq peak to compare among cell types
atac_fl = list.files(path = paste0("../data/ATAC/"), pattern = paste0(".bdg.gz.sum$"), full.names = T)
for (i in 1:length(atac_fl)){
  tryCatch({
    if (i==1){
      tryCatch({
        #print(atac_fl[i])
        sum_atac = fread(atac_fl[i]) %>%
          dplyr::mutate(isoform = stringr::str_split(atac_fl[i],pattern="/",simplify = TRUE) %>% .[,5] %>% gsub("_merged_treat_pileup.bdg.gz.sum","",.) %>% gsub("_no_treament"," [unstim]", .) %>% gsub("_treament"," [stim]", .) %>% gsub("_"," ",.)) %>%
          dplyr::rename(sum_count = V1)
      }, error=function(e){cat("")})
    } else {
      tryCatch({
        #print(atac_fl[i])
        sum_atac = fread(atac_fl[i]) %>%
          dplyr::mutate(isoform = stringr::str_split(atac_fl[i],pattern="/",simplify = TRUE) %>% .[,5] %>% gsub("_merged_treat_pileup.bdg.gz.sum","",.) %>% gsub("_no_treament"," [unstim]", .) %>% gsub("_treament"," [stim]", .) %>% gsub("_"," ",.)) %>%
          dplyr::rename(sum_count = V1) %>%
          rbind(.,sum_atac)
      }, error=function(e){cat("")})
    }
  }, error=function(e){cat("")})
}

plot_iso_structure <- function(gene = NULL, extra = 5e03, repeat_draw = TRUE, atac_draw = TRUE, rip_draw = TRUE){
  
  region_chrom = isoform_Info %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    .$chrom %>% unique()
  region_start = isoform_Info %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    .$start %>% min() %>% -extra
  region_end = isoform_Info %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    .$end %>% max() %>% +extra
  
  domains = rbind(sigp[,c("isoform","associated_gene","predicted_NMD","chrom","start","end","strand","coding")] %>%
                    dplyr::mutate(hmm_name = "signal peptide"),
                  pfam[,c("isoform","associated_gene","predicted_NMD","chrom","start","end","hmm_name","strand","coding")]
  ) %>%
    rbind(.,
          iupred[,c("isoform","associated_gene","predicted_NMD","chrom","start","end","group","strand","coding")] %>%
            dplyr::rename(hmm_name = group)) %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    dplyr::rename(domain = hmm_name) %>%
    dplyr::arrange(chrom,start) %>%
    dplyr::mutate(region = "domain",
                  anno = "Isoform Atlas",
                  isoform = ifelse(strand == "+", paste0("->",isoform), paste0("<-",isoform)),
                  height = NA) %>%
    as.data.frame()
  
  # GENCODE38
  coord_gencode_orf = coord_gencode %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    dplyr::filter(region == "CDS") %>%
    dplyr::group_by(associated_gene) %>%
    valr::bed_merge() %>%
    dplyr::mutate(region = "CDS")
  coord_gencode_exon = coord_gencode %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    dplyr::filter(region == "exon") %>%
    dplyr::group_by(associated_gene) %>%
    valr::bed_merge() %>%
    dplyr::mutate(region = "exon")
  coord_gencode_intron = coord_gencode %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    dplyr::filter(region == "intron") %>%
    dplyr::group_by(associated_gene) %>%
    valr::bed_merge() %>%
    dplyr::mutate(region = "intron")
  coord_gencode_merged = rbind(coord_gencode_orf,coord_gencode_exon) %>%
    rbind(.,coord_gencode_intron) %>%
    dplyr::mutate(anno = "GENCODE",
                  domain = NA,
                  isoform = associated_gene,
                  height = NA) %>%
    as.data.frame()
  coord_ia = coord %>%
    dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
    dplyr::mutate(isoform = ifelse(strand == "+", paste0("->",isoform), paste0("<-",isoform)),
                  anno = "Isoform Atlas",
                  domain = NA,
                  height = NA) %>%
    as.data.frame()
  rm_tmp = rm %>%
    dplyr::filter(start > region_start & end < region_end)
  
  if (rip_draw){
    # RIP-seq of ELAVL1 [GSM944520]
    print("extracting RIP-seq peak...")
    range = paste0(region_chrom,":",region_start,"-",region_end)
    rip = seqminer::tabix.read.table(tabixFile = "../data/genome_uniq_for_UCSC.bg.gz", tabixRange = range, stringsAsFactors = FALSE) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(isoform = "ELAVL1 RIP-seq",
                    type = "ELAVL1 RIP-seq",
                    anno = "RIP",
                    domain = NA,
                    associated_gene = NA,
                    region = "RIP") %>%
      as.data.frame()
    colnames(rip) <- c("chr","start","end","height","isoform","type","anno","domain","associated_gene","region")
  }

  if (atac_draw){
    # ATAC-seq [Nat Genet 51, 1494–1505 (2019).]
    print("extracting ATAC-seq peak...")
    range = paste0(region_chrom,":",region_start,"-",region_end)
    ATAC_summary_list = purrr::map(fl, 
                                   ~safe_import_ATAC(., region = range))
    names(ATAC_summary_list) = cell_name
    result_list = purrr::map(ATAC_summary_list, ~.$result)
    result_list = result_list[!unlist(purrr::map(result_list, is.null))]
    dummy_atac = purrr::map_df(result_list, ~as.data.frame(.x)) %>%
      merge(.,sum_atac,by="isoform") %>%
      dplyr::mutate(height = 1e10*height/sum_count,
                    height = max(height),
                    anno = "ATAC",
                    domain = NA,
                    associated_gene = NA,
                    region = "ATAC_dummy",
                    chrom = region_chrom)
    atac = purrr::map_df(result_list, ~as.data.frame(.x)) %>%
      merge(.,sum_atac,by="isoform") %>%
      dplyr::mutate(height = 1e10*height/sum_count,
                    anno = "ATAC",
                    domain = NA,
                    associated_gene = NA,
                    region = "ATAC",
                    chrom = region_chrom)
  }
  if (rip_draw){
    if (atac_draw){
      shared = intersect(colnames(coord_ia),colnames(coord_gencode_merged)) %>% intersect(.,colnames(domains)) %>% intersect(.,colnames(atac)) %>% intersect(.,colnames(rip))
    } else {
      shared = intersect(colnames(coord_ia),colnames(coord_gencode_merged)) %>% intersect(.,colnames(domains)) %>% intersect(.,colnames(rip))
    }
  } else {
    if (atac_draw){
      shared = intersect(colnames(coord_ia),colnames(coord_gencode_merged)) %>% intersect(.,colnames(domains)) %>% intersect(.,colnames(atac))
    } else {
      shared = intersect(colnames(coord_ia),colnames(coord_gencode_merged)) %>% intersect(.,colnames(domains))
    }  
  }
  
  if (rip_draw){
    if (atac_draw){
      if (repeat_draw){
        isoform_list = c("ELAVL1 RIP-seq",
                         coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         rm_tmp[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         paste0(rep(c("Bulk B","Naive B","Mem B","Plasmablasts","CD8pos T","Naive CD8 T","Central memory CD8pos T","Effector memory CD8pos T","Effector CD4pos T",  
                                      "Naive Teffs","Memory Teffs","Th1 precursors","Th2 precursors","Th17 precursors","Follicular T Helper","Regulatory T","Naive Tregs","Memory Tregs","Gamma delta T",
                                      "Immature NK", "Mature NK","Memory NK",                
                                      "Monocytes","Myeloid DCs","pDCs"),each = 2),rep(c(" [unstim]", " [stim]"),each=1))
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      } else {
        isoform_list = c("ELAVL1 RIP-seq",
                         coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         # rm_tmp[,shared] %>%
                         #   .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         paste0(rep(c("Bulk B","Naive B","Mem B","Plasmablasts","CD8pos T","Naive CD8 T","Central memory CD8pos T","Effector memory CD8pos T","Effector CD4pos T",  
                                      "Naive Teffs","Memory Teffs","Th1 precursors","Th2 precursors","Th17 precursors","Follicular T Helper","Regulatory T","Naive Tregs","Memory Tregs","Gamma delta T",
                                      "Immature NK", "Mature NK","Memory NK",                
                                      "Monocytes","Myeloid DCs","pDCs"),each = 2),rep(c(" [unstim]", " [stim]"),each=1))
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      }
    } else {
      if (repeat_draw){
        isoform_list = c("ELAVL1 RIP-seq",
                         coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         rm_tmp[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort()
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      } else {
        isoform_list = c("ELAVL1 RIP-seq",
                         coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort()
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,rip[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "RIP"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            ggh4x::facet_nested(factor(anno,levels=c("RIP","Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      }
    }
    
  } else {
    if (atac_draw){
      if (repeat_draw){
        isoform_list = c(coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         rm_tmp[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         paste0(rep(c("Bulk B","Naive B","Mem B","Plasmablasts","CD8pos T","Naive CD8 T","Central memory CD8pos T","Effector memory CD8pos T","Effector CD4pos T",  
                                      "Naive Teffs","Memory Teffs","Th1 precursors","Th2 precursors","Th17 precursors","Follicular T Helper","Regulatory T","Naive Tregs","Memory Tregs","Gamma delta T",
                                      "Immature NK", "Mature NK","Memory NK",                
                                      "Monocytes","Myeloid DCs","pDCs"),each = 2),rep(c(" [unstim]", " [stim]"),each=1))
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      } else {
        isoform_list = c(coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         # rm_tmp[,shared] %>%
                         #   .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         paste0(rep(c("Bulk B","Naive B","Mem B","Plasmablasts","CD8pos T","Naive CD8 T","Central memory CD8pos T","Effector memory CD8pos T","Effector CD4pos T",  
                                      "Naive Teffs","Memory Teffs","Th1 precursors","Th2 precursors","Th17 precursors","Follicular T Helper","Regulatory T","Naive Tregs","Memory Tregs","Gamma delta T",
                                      "Immature NK", "Mature NK","Memory NK",                
                                      "Monocytes","Myeloid DCs","pDCs"),each = 2),rep(c(" [unstim]", " [stim]"),each=1))
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,atac[,shared]) %>%
            rbind(.,dummy_atac[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_ribbon(data = . %>%
                          dplyr::filter(region == "ATAC"),
                        aes(x = start, ymin = 0, ymax = height, 
                            # fill = isoform, 
                            color = isoform),
                        alpha = 0.5) +
            geom_blank(data = . %>%
                         dplyr::filter(region == "ATAC_dummy"),
                       aes(x = start, ymin = 0, ymax = height)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE","ATAC")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      }
    } else {
      if (repeat_draw){
        isoform_list = c(coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         rm_tmp[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort()
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,rm_tmp[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "repetitive elements"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "forestgreen") +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      } else {
        isoform_list = c(coord_ia[,shared] %>%
                           .$isoform %>% unique() %>% sort(),
                         coord_gencode_merged[,shared] %>%
                           .$isoform %>% unique() %>% sort()
        )
        if (nrow(domains)>0 &
            nrow(coord %>%
                 dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                 dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)>0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            rbind(.,domains[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "domain"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6, fill = domain)) +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))>0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "CDS"),
                      aes(xmin = start, xmax = end, ymin = 4, ymax = 6),
                      fill = "gray0") +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        } else if (nrow(domains)==0 &
                   nrow(coord %>%
                        dplyr::filter((associated_gene == gene)|grepl(paste0("^",gene,"_"),associated_gene)|grepl(paste0("_",gene,"$"),associated_gene)) %>%
                        dplyr::filter(region == "CDS"))==0){
          print("plotting...")
          g = rbind(coord_ia[,shared],coord_gencode_merged[,shared]) %>%
            ggplot(.) +
            geom_segment(aes(x = start, xend = end, y = 5, yend = 5), size = 0.5, linetype="dashed", color = "gray70") +
            geom_rect(data = . %>%
                        dplyr::filter(region == "exon"),
                      aes(xmin = start, xmax = end, ymin = 4.5, ymax = 5.5),
                      fill = "grey45") +
            ggh4x::facet_nested(factor(anno,levels=c("Isoform Atlas","repetitive elements","GENCODE")) + factor(isoform,levels = isoform_list) ~ ., scales = "free", nest_line = TRUE) +
            theme(strip.text.y = element_text(angle = 0)) +
            xlab(paste0("genomic position [", region_chrom, "]")) +
            ylab("") +
            labs(title = paste0("isoform structure [", gene, "]")) +
            scale_y_continuous(breaks=NULL) +
            guides(color = FALSE) +
            theme_minimal() +
            theme(strip.text.x=element_text(size=9, color="black"),
                  strip.text.y=element_text(angle = 0, size=10, color="black"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  plot.title = element_text(size=10),
                  axis.title.x = element_text(size=10),
                  axis.title.y = element_text(size =10),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text =  element_text(size = 10), 
                  legend.key.size = grid::unit(0.8, "lines"),
                  legend.title = element_text(size = 10, hjust = 0),
                  legend.direction = "horizontal", legend.box = "horizontal")
        }
      }
    }
  }
  return(g)
}



plot_iso_usage <- function(gene = NULL,
                           method = c("ratio","sum"),
                           legend = c("top","bottom","right","left","none")){
  library(ggplot2)
  library(ggsci)
  library(matrixStats)
  
  if(length(read2gene[read2gene$GeneSymbol == gene, ]$transcript_id) >1 ){
    tmp = cts[rownames(cts) %in% read2gene[read2gene$GeneSymbol == gene, ]$transcript_id, ] %>%
      as.data.frame()
  } else if (length(read2gene[read2gene$GeneSymbol == gene, ]$transcript_id) ==1){
    tmp = cts[rownames(cts) %in% read2gene[read2gene$GeneSymbol == gene, ]$transcript_id, ] %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
  }
  
  if (nrow(tmp) > 12){
    label = tmp %>%
      dplyr::mutate_all(funs(./sum(.))) %>%
      dplyr::mutate_all(funs(ifelse(is.na(.),0,.))) %>%
      dplyr::mutate_all(funs(ifelse(is.infinite(.),0,.))) %>% 
      dplyr::mutate_all(funs(ifelse(is.nan(.),0,.))) %>% 
      as.matrix() %>%
      matrixStats::rowMaxs(.) %>%
      rank() %>%
      ifelse(.>nrow(tmp)-12,"major")
    major = tmp[label,]
    other = tmp[-label,] %>%
      colSums(.) %>%
      t() %>%
      magrittr::set_rownames("others")
    tmp = rbind(major,other)
  }
  
  if (method == "sum"){
    label = "read count(RPM)"
    cts_long = 
      tmp %>%
      dplyr::mutate(isoform = rownames(.)) %>%
      reshape2::melt(.,
                     id.vars=c("isoform"),
                     variable.name="cell",
                     value.name="read_count",
                     na.rm=TRUE) %>%
      dplyr::mutate(read_count = read_count,
                    isoform = factor(isoform,levels=c(unique(rownames(tmp))[!grepl("others",unique(rownames(tmp)))],"others")))
  } else if (method == "ratio"){
    label = "isoform ratio (%)"
    cts_long = 
      tmp %>%
      dplyr::mutate_all(funs(./sum(.))) %>%
      dplyr::mutate_all(funs(ifelse(is.na(.),0,.))) %>%
      dplyr::mutate_all(funs(ifelse(is.infinite(.),0,.))) %>% 
      dplyr::mutate_all(funs(ifelse(is.nan(.),0,.))) %>% 
      dplyr::mutate(isoform = rownames(.)) %>%
      reshape2::melt(.,
                     id.vars=c("isoform"),
                     variable.name="cell",
                     value.name="read_count",
                     na.rm=TRUE) %>%
      dplyr::mutate(read_count = read_count*100 %>% as.numeric(),
                    isoform = factor(isoform,levels=c(unique(rownames(tmp))[!grepl("others",unique(rownames(tmp)))],"others")))
  }
  
  g <- ggplot(data=cts_long, aes(x=cell, y = read_count, fill=isoform)) + 
    geom_bar(stat = "identity") + 
    xlab("") +
    ylab(label) +
    coord_flip() +
    theme_minimal(base_size = 20) +
    scale_fill_manual(values = c(unique(c(ggsci::pal_nejm("default", alpha = .7)(8)[2:1],
                                          ggsci::pal_npg("nrc", alpha = .7)(10)))[1:(nrow(tmp)-1)],"grey45")) +
    theme(strip.text.x=element_text(size=9, color="black"),
          strip.text.y=element_text(angle = 0, size=9, color="black"),
          panel.grid = element_blank(),
          legend.position = legend,
          plot.title = element_text(size=15),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size =20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.text =  element_text(size = 20), 
          legend.key.size = grid::unit(0.8, "lines"),
          legend.title = element_text(size = 0, hjust = 0),
          legend.direction = "vertical", legend.box = "horizontal")
  return(g)
}
