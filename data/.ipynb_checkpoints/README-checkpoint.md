## colum of "isoform_info.txt"
### isoform: isoform ID
### associated_gene: gene symbol
### chrom: chrosome
### start: start position
### end: end position
### strand: strand
### five_utr_length: length of 5’UTR region
### CDS_length: length of coding region
### three_utr_length: length of 3’UTR region
### polyA_motif: motif of poly A signal (“no-PAS” means no canonical motif)
### kozak_score: kozak score [This is G c c A/G c c atg G. The most important nts are +4, -3 and -6.  Scoring these as +3 and the others as +1. Max score = 13]
### avg_codon_freq: codon frequency averaged across CDS (codon table: http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N retrieved on 11/20/2014)
### au_element_count: number of AU-stretches
### au_element_frac: percentage of UTR covered by AREs 
### max_au_length: longest A/U stretch
### fiveUTRcap_MFE: minimum folding energy at 5' end (for 5' UTR, specifically affects 43S loading). This is calculated using the sequence of the 50nt after the 5' end, or if the 5' UTR is less than 50nt just calculate using the whole 5’UTR sequence, using viennaRNA
### unique_TSS: transcription start site is specific to the isoform only
### unique_CDS coding sequence is specific to the isoform only
### unique_FE: no overlap with first exon of other isoforms
### unique_LE: no overlap with last exon of other isoforms
### te_rank: ranking according to translation efficiency (top10, others, and bottom 10: e.g., top10 means top 10% of translational efficiency). Translation efficiency is calculated using samples from 52 Yoruba (ribo-seq [GSE61742] and RNA-seq [GEUVADIS cohort, Nature 2013;501:506–511.]) 
### InnateDB: immune genes annotated by InnateDB (https://www.innatedb.com/annotatedGenes.do?type=innatedb)
### TF: transcription factors
### transmembrane: transmembrane proteins
### signal_peptide: isoforms containing signal peptide sequence
### idr: isoforms containing intrinsically disordered protein region
### anchor2: isoforms containing intrinsically disordered binding region
### uORF: isoforms containing predicted upstream open reading frame using ribo-TISH (ribo-seq datasets were downloaded from GSE39561, GSE56887, GSE61742, GSE74279, GSE75290, GSE81802, and GSE97140)
### predicted_NMD: isoforms predicted to cause nonsense-mediated decay 
### specificity_LR: specifically expressed isoforms in any of the long-read sequenced 29 cell-subsets based on both expression and transcript ratio using ROKU function in TCC package
### specific_cell_LR: specifically expressed cell in any of the long-read sequenced 29 cell-subsets
### specificity_LRgroup: specifically expressed isoforms in any of the long-read sequenced 8 cell-groups based on expression and transcript ratio using ROKU function in TCC package
### specific_cell_LRgroup: specifically expressed group in any of the long-read sequenced 8 cell-groups
### •	CD4T: NaiveCD4,Th1,Th2,Th17,Tfh,Fra1.Treg,Fra2.aTreg,Fra3.Treg,LAG3.Treg,MemoryCD4,Thx
### •	CD8T: NaiveCD8,CD8effector,CD8centralmem,CD8effectormem
### •	B: NaiveB,unswmemoryB,swmemoryB,DNB,plasmablast"
### •	DC: myeloidDC,plasmacytoidDC
### •	NK: NK
### •	monocyte: monocyteCD16,monocyteCD16minus,nonclassicalMonocyte,intermediateMonocyte
### •	PBMC: PBMC
### •	Neutrophil: Neutrophil
### specificity_SR_celltype: specifically expressed isoforms in any of the short-read sequenced nonstimulated cell-conditions based on both expression and transcript ratio using ROKU function in TCC package
### specificity_SR_stim: specifically expressed isoforms in any of the short-read sequenced stimulated cell-conditions based on both expression and transcript ratio using ROKU function in TCC package
### specific_cell_SR: specifically expressed cell in any of the short-read sequenced cell-subsets
### matching_repeat: repetitive elements contained in the isoform
### coloc_eQTL: colocalization between eQTL signal of associated gene and any of GWAS signal
### cell_disease_eQTL: cell condition and phenotype of colocalization
### coloc_sQTL: colocalization between sQTL signal of the isoform and any of GWAS signal
### cell_disease_sQTL: cell condition and phenotype of colocalization
### coloc_5UTR_QTL: colocalization between alternative aTSS-QTL signal of the isoform and any of GWAS signal
### cell_disease_5UTR_QTL: cell condition and phenotype of colocalization
### coloc_3UTR_QTL: colocalization between alternative polyadenylation-QTL signal of the isoform and any of GWAS signal
### cell_disease_3UTR_QTL: cell condition and phenotype of colocalization
