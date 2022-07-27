#!/usr/bin/Rscript

# BiocManager::install("IsoformSwitchAnalyzeR")
# BiocManager::install("rhdf5")
# BiocManager::install("rhdf5filters")

library(data.table, lib.loc = "/data02/home/juninamo/miniconda3/envs/IsoformSwitchAnalyzeR/lib/R/library")
library(tidyverse, lib.loc = "/data02/home/juninamo/miniconda3/envs/IsoformSwitchAnalyzeR/lib/R/library")
library(BSgenome.Hsapiens.NCBI.GRCh38, lib.loc = "/data02/home/juninamo/miniconda3/envs/IsoformSwitchAnalyzeR/lib/R/library")
library(BiocParallel, lib.loc = "/data02/home/juninamo/miniconda3/envs/IsoformSwitchAnalyzeR/lib/R/library")
register(MulticoreParam(30))
library(IsoformSwitchAnalyzeR, lib.loc = "/data02/home/juninamo/miniconda3/envs/IsoformSwitchAnalyzeR/lib/R/library")
library(rhdf5, lib.loc = "/data02/home/juninamo/miniconda3/envs/IsoformSwitchAnalyzeR/lib/R/library")
packageVersion('IsoformSwitchAnalyzeR')

kallistoQuant <- importIsoformExpression(
    parentDir = "/path/to/kallisto"
)

print(head(kallistoQuant$abundance, 2))
print(head(kallistoQuant$counts, 2))

myDesign <- fread("/path/to/SraRunTable.txt",header = TRUE) %>%
  dplyr::mutate(condition = stringr::str_split(source_name, pattern=", ", simplify=TRUE) %>% .[,2] %>% stringr::str_split(., pattern=" ", simplify=TRUE) %>% .[,1]) %>%
  dplyr::select(c(Run,condition)) %>%
  magrittr::set_colnames(c("sampleID","condition"))

aSwitchList <- importRdata(
    isoformCountMatrix   = kallistoQuant$counts,
    isoformRepExpression = kallistoQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.gtf",
    isoformNtFasta       = "/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fasta",
    fixStringTieAnnotationProblem = TRUE,
    showProgress = FALSE
)

print(summary(aSwitchList))

exampleSwitchListPart1 <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist   = aSwitchList,
    pathToOutput         = '/path/to/IsoformSwitchAnalyzeR/',
    outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
    prepareForWebServers = FALSE
)

print(extractSwitchSummary( exampleSwitchListPart1 ))

save( exampleSwitchListPart1, file = paste0("/path/to/IsoformSwitchAnalyzeR/exampleSwitchListPart1.RData") )

load(file = paste0("/path/to/IsoformSwitchAnalyzeR/exampleSwitchListPart1.RData") )

exampleSwitchListPart2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = exampleSwitchListPart1, 
  codingCutoff              = 0.364,
  removeNoncodinORFs        = FALSE,
  pathToCPATresultFile      = "/path/to/IsoformSwitchAnalyzeR/cpat_all.result",
  pathToPFAMresultFile      = "/path/to/IsoformSwitchAnalyzeR/pfam_all.result",
  pathToIUPred2AresultFile  = '/path/to/IsoformSwitchAnalyzeR/iupred_all.result',
  pathToSignalPresultFile   = '/path/to/IsoformSwitchAnalyzeR/signalp_all.result',

  Analysis and output arguments
  n = Inf,
  consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'IDR_identified',
        'IDR_type',
        'signal_peptide_identified'
  ),
  pathToOutput = '/path/to/IsoformSwitchAnalyzeR/',
  fileType = 'pdf',
  outputPlots = TRUE,

  Other arguments
  quiet = FALSE
)

save( exampleSwitchListPart2, file = paste0("/path/to/IsoformSwitchAnalyzeR/exampleSwitchListPart2.RData") )
