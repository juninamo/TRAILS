#!/usr/bin/R

chr = commandArgs(trailingOnly=TRUE)[1]
# chr=21

txdb = GenomicFeatures::makeTxDbFromGFF(paste0("/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.utr.start.stop.sorted.chr",chr,".gtf"),
                                        format="gtf",
                                        organism="Homo sapiens")
RPFs <- dir(paste0("/path/to/GEUV/TE/BAM/ribo/chr",chr), ".bam$", full.names=TRUE)
RNAs <- dir(paste0("/path/to/GEUV/TE/BAM/rna/chr",chr), ".bam$", full.names=TRUE)
print(RPFs)
print(RNAs)
print("Calculating coverageDepth...")
cvgs <- ribosomeProfilingQC::coverageDepth(RPFs, RNAs, txdb)
print("Calculating Translation Efficiency...")
TE90 <- ribosomeProfilingQC::translationalEfficiency(cvgs, window = 90, normByLibSize=TRUE)

write.table(TE90$TE,paste0("/path/to/GEUV/TE/translation_efficiency_chr",chr,".txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
write.table(TE90$mRNA,paste0("/path/to/GEUV/TE/coverage_mRNA_chr",chr,".txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
write.table(TE90$RPFs,paste0("/path/to/GEUV/TE/coverage_RPFs_chr",chr,".txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)

print("ALL done!")
