# Immune Isoform Atlas
[![DOI](https://zenodo.org/badge/518259643.svg)](https://zenodo.org/badge/latestdoi/518259643)
![](https://komarev.com/ghpvc/?username=juninamo&style=flat-square&color=green&label=REPOSITORY+VIEWS)

## Full-length transcript annotation focusing on human immune cells

Alternative splicing events are a major causal mechanism for complex traits, but they have been understudied due to the limitation of short-read sequencing. Here, we generated a comprehensive full-length isoform annotation of human immune cells, Immune Isoform Atlas, by long-read sequencing for 29 cell subsets. Our atlas contained a number of unannotated transcripts and functional characteristics of transcripts including encoded domains, inserted repetitive elements, cell-type specific expression, and translational efficiency. Further, we identified a number of disease-associated isoforms by isoform-switch analysis and by integration of several quantitative trait loci analyses with genome-wide association study data. These results are open on [the web](http://gfdweb.tmd.ac.jp:3838/?) and the [genome browser](https://genome.ucsc.edu/s/juninamo/Immune%20Isoform%20Atlas).

## Citation 
Jun Inamo, Akari Suzuki, Mahoko Ueda, Kensuke Yamaguchi, Hiroshi Nishida, Katsuya Suzuki, Yuko Kaneko, Tsutomu Takeuchi, Yasushi Ishihama, Kazuhiko Yamamoto, Yuta Kochi. Immune Isoform Atlas: Landscape of alternative splicing in human immune cells. [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2022.09.13.507708v1), doi:[10.1101/2023.07.03.547507](https://www.biorxiv.org/content/10.1101/2022.09.13.507708v1)

<kbd>
<img src="https://github.com/juninamo/isoform_atlas/blob/main/images/project_image.png" width="800" align="center">
</kbd>

&nbsp;&nbsp;

**Figure 1. Study design**

&nbsp;&nbsp;

## Sequenced cell-subsets
|  Subset name  |  Abbreviation  |
| :---: | :---: |
|  Naïve CD4 T cells |  Naïve CD4  |
|  Memory CD4 cells |  Mem CD4  |
|  Fraction I  naive regulatory T cells |  FraI nTreg  |
|  Fraction II  effector regulatory T cells |  FraII aTreg  |
|  Fraction III  non-regulatory T cells |  FraIII non-Treg  |
|  [Low-Density Granulocytes regulatory T cells](https://www.nature.com/articles/ncomms7329) |  LAG3 Treg  |
|  T helper 1 cells |  Th1  |
|  T helper 2 cells |  Th2  |
|  T helper 17 cells |  Th17  |
|  [CXCR3 +/−CCR6− T cells](https://www.jimmunol.org/content/200/6/2090.long) |  X3lowR6negT  |
|  T follicular helper cells |  Tfh  |
|  Naïve CD8 T cells |  Naïve CD8  |
|  Effector CD8 cells |  Eff CD8  |
|  Central Memory CD8 T cells |  CM CD8  |
|  Effector Memory CD8 T cells |  EM CD8  |
|  Naïve B cells |  Naïve B  |
|  Switched memory B cells |  SM B  |
|  Unswitched memory B cells | USM B  |
|  Double Negative B cells | DN B  |
|  Plasmablasts | Plasmablast  |
|  Natural Killer cells | NK  |
|  Classical Monocytes | CL Mono  |
|  NonClassical Monocytes | NC Mono  |
|  Intermediate Monocytes | Int Mono  |
|  CD16 positive Monocytes | CD16p Mono  |
|  Myeloid Dendric cells | mDC  |
|  Plasmacytoid Dendric cells | pDC  |
|  Neutroohils | Neu  |
|  Peripheral blood mononuclear cells | PBMC  |


## Highlights
- An atlas of full-length isoforms for 29 immune cell types.
- Transposable elements comprise a major fraction of isoform diversity.
- Alternative 3’UTR usage contributes to cell-type specific expression of isoforms.
- Integrated analysis of genetic and transcriptomic data with the atlas reveals unknown pathogenesis of diseases.

## How can users utilize Isoform Atlas?
- [User-friendly web app is available](http://gfdweb.tmd.ac.jp:3838/?)
- Users can remapping own RNA-seq datasets to Immune Isoform Atlas (isoform_atlas.gtf.gz, GRCh38) and investigate expression of isoforms in interested phenotypes.

## Files
- isoform_atlas.gtf.gz: GTF file of Immune Isoform Atlas.
- isoform_info.txt.gz: Detailed information for each isoform in Immune Isoform Atlas.
- male_PBMC.gtf.gz: GTF file of independent dataset by long-read RNA sequencing using the ONT platform (PromethION R10.4.1, V14 chemistry). We used this dataset to validate isoforms in Immune Isoform Atlas and investigate sex-differences.
- IsoQuant.transcript_models.gtf.gz: GTF file using IsoQuant pipeline (Prjibelski, A.D., et al. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol (2023).)
- IsoQuant.transcript_models.extended_annotation.gtf.gz: GTF file using IsoQuant pipeline (Extended).
- Espresso.gtf.gz: GTF file using ESPRESSO pipeline (without junction correction) (Gao Y, et al. ESPRESSO: Robust discovery and quantification of transcript isoforms from error-prone long-read RNA-seq data. Sci Adv. 2023 Jan 20;9(3):eabq5072.).
- Espresso_SJ.gtf.gz: GTF file using ESPRESSO pipeline (with junction correction) 
- code/: Source codes for each analysis.
- Figures/: Source codes for generating figures in our paper.

## Colums of "isoform_info.txt"
- isoform: isoform ID
- associated_gene: gene symbol
- chrom: chrosome
- start: start position
- end: end position
- strand: strand
- 5'UTR_length: length of 5’UTR region
- ORF_length: length of coding region
- 3'UTR_length: length of 3’UTR region
- polyA_motif: motif of poly A signal (“no-PAS” means no canonical motif)
- kozak_score: kozak score [This is G c c A/G c c atg G. The most important nts are +4, -3 and -6.  Scoring these as +3 and the others as +1. Max score = 13]
- avg_codon_freq: codon frequency averaged across CDS ([codon table](http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N) retrieved on 11/20/2014)
- AU_element_count: number of AU-stretches
- AU_element_frac: percentage of UTR covered by AREs 
- max_AU_length: longest A/U stretch
- 5'UTRcap_MFE: minimum folding energy at 5' end (for 5' UTR, specifically affects 43S loading). This is calculated using the sequence of the 50nt after the 5' end, or if the 5' UTR is less than 50nt just calculate using the whole 5’UTR sequence, using viennaRNA
- unique_TSS: transcription start site is specific to the isoform only
- unique_ORF coding sequence is specific to the isoform only
- unique_FE: no overlap with first exon of other isoforms
- unique_LE: no overlap with last exon of other isoforms
- translational_efficiency_rank: ranking according to translation efficiency (top10, others, and bottom 10: e.g., top10 means top 10% of translational efficiency). Translation efficiency is calculated using samples from 52 Yoruba (ribo-seq [GSE61742] and RNA-seq [GEUVADIS cohort, Nature 2013;501:506–511.]) 
- [immune_genes: immune genes annotated by InnateDB](https://www.innatedb.com/annotatedGenes.do?type=innatedb)
- TF: [transcription factors](http://humantfs.ccbr.utoronto.ca/download.php)
- transmembrane: transmembrane proteins
- signal_peptide: isoforms containing signal peptide sequence
- IDR: isoforms containing intrinsically disordered protein region
- ANCHOR2: isoforms containing intrinsically disordered binding region
- uORF: isoforms containing predicted upstream open reading frame using ribo-TISH (ribo-seq datasets were downloaded from GSE39561, GSE56887, GSE61742, GSE74279, GSE75290, GSE81802, and GSE97140)
- predicted_NMD: isoforms predicted to cause nonsense-mediated decay 
- specificity_LR: specifically expressed isoforms in any of the long-read sequenced 29 cell-subsets based on both expression and transcript ratio using ROKU function in TCC package
- specific_cell_LR: specifically expressed cell in any of the long-read sequenced 29 cell-subsets
- specificity_LRgroup: specifically expressed isoforms in any of the long-read sequenced 8 cell-groups based on expression and transcript ratio using ROKU function in TCC package
- specific_cell_LRgroup: specifically expressed group in any of the long-read sequenced 8 cell-groups
    - CD4T: NaiveCD4, Th1, Th2, Th17, Tfh, FraI nTreg, FraII aTreg, FraIII non-Treg, LAG3 Treg, Mem CD4, X3lowR6negT
    - CD8T: NaiveCD8, Eff CD8, CM CD8, EM CD8 
    - B: NaiveB, USMB, SMB,DNB, plasmablast
    - DC: mDC, pDC
    - NK: NK
    - monocyte: CL Mono, NC Mono, Int Mono, CD16p Mono 
    - PBMC: PBMC
    - Neutrophil: Neu
- repeat_elements: repetitive elements contained in the isoform
- coloc_eQTL: colocalization between eQTL signal of associated gene and any of GWAS signal
- cell_disease_eQTL: cell condition and phenotype of colocalization
- coloc_sQTL: colocalization between sQTL signal of the isoform and any of GWAS signal
- cell_disease_sQTL: cell condition and phenotype of colocalization
- coloc_3'aQTL: colocalization between 3'aQTL signal of the isoform and any of GWAS signal
- cell_disease_3'aQTL: cell condition and phenotype of colocalization
- male_PBMC: validated isoforms by PBMC sample (PromethION R10.4.1, V14 chemistry) from 40 y/o male 


## Contact us
Please contact us (Jun Inamo: juninamo.gfd@mri.tmd.ac.jp) with any questions or comments.

The data presented here comes from the laboratory of [Yuta Kochi](https://www.tmd.ac.jp/gfd/english/) through collaborating with the laboratory of [Yasushi Ishihama](https://www.pharm.kyoto-u.ac.jp/seizai/index_e.html), [RIKEN](https://www.riken.jp/en/), and [Keio University](https://www.keio.ac.jp/en/).

## Acknowledgments
Computations were partially performed on the NIG supercomputer at the ROIS National Institute of Genetics. This study was supported by (the Japan Society for the Promotion of Science)[https://www.jsps.go.jp/english/], [the MEXT of Japan](https://www.mext.go.jp/en/), and grants from [Nanken-Kyoten, TMDU](https://www.tmd.ac.jp/english/) and Medical Research Center Initiative for High Depth Omics.

<kbd>
<img src="https://github.com/juninamo/isoform_atlas/blob/main/images/README_logo.png" width="800" align="center">
</kbd>

&nbsp;&nbsp;

