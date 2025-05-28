#!/bin/bash
##$ -S /bin/sh

# GEUV
cohort="GEUV"
mkdir /path/to/${cohort}/QTL/result/p1e05/splicing_peer
for cell in EUR; do
for tool in leafcutter; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i chrom\tstart_clu\tend_clu\tstrand_clu\tsnp\tclu\tstatistic\tpval\tFDR\tbeta\tstart_gene\tend_gene\tclass\tgene\toverlap\tdist\tA1\tA0" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
for tool in StringTie kallisto kallisto-cds-clusterd-vsearch; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i gene\tsnps\tstatistic\tpvalue\tFDR\tbeta\tA1\tA0\tGENEID\tchr\tstart\tend\tstrand\tgene_id" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
done

# EGA
cohort="EGA"
mkdir /path/to/${cohort}/QTL/result/p1e05/splicing_peer
for cell in Mono_NS Mono_LPS Mono_Pam3CSK4 Mono_R848 Mono_IAV; do
for tool in leafcutter; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i chrom\tstart_clu\tend_clu\tstrand_clu\tsnp\tclu\tstatistic\tpval\tFDR\tbeta\tstart_gene\tend_gene\tclass\tgene\toverlap\tdist\tA1\tA0" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
for tool in StringTie kallisto kallisto-cds-clusterd-vsearch; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i gene\tsnps\tstatistic\tpvalue\tFDR\tbeta\tA1\tA0\tGENEID\tchr\tstart\tend\tstrand\tgene_id" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
done

# ImmVar
cohort="ImmVar"
mkdir /path/to/${cohort}/QTL/result/p1e05/splicing_peer
for cell in CD4T_Activated MoDC_unstim MoDC_FLU MoDC_IFNb; do
for tool in leafcutter; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i chrom\tstart_clu\tend_clu\tstrand_clu\tsnp\tclu\tstatistic\tpval\tFDR\tbeta\tstart_gene\tend_gene\tclass\tgene\toverlap\tdist\tA1\tA0" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
for tool in StringTie kallisto kallisto-cds-clusterd-vsearch; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i gene\tsnps\tstatistic\tpvalue\tFDR\tbeta\tA1\tA0\tGENEID\tchr\tstart\tend\tstrand\tgene_id" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
done

# DICE
cohort="DICE"
mkdir /path/to/${cohort}/QTL/result/p1e05/splicing_peer
for cell in TFH TH1 TH17 TH1-17 TH2 TREG_MEMORY TREG_NAIVE B_NAIVE CD4_NAIVE CD4_N_STIM CD8_NAIVE CD8_N_STIM NONCLASSICAL_MONOCYTES CLASSICAL_MONOCYTES NK_CD16POS; do
for tool in leafcutter; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i chrom\tstart_clu\tend_clu\tstrand_clu\tsnp\tclu\tstatistic\tpval\tFDR\tbeta\tstart_gene\tend_gene\tclass\tgene\toverlap\tdist\tA1\tA0" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
for tool in StringTie kallisto kallisto-cds-clusterd-vsearch; do
echo -e "${cohort}, ${cell}, ${tool}" 
awk '$8<1e-05 {print $0}' /path/to/${cohort}/QTL/result/splicing_peer/chr*_in_${cell}_${tool}.txt | sed -e "1i gene\tsnps\tstatistic\tpvalue\tFDR\tbeta\tA1\tA0\tGENEID\tchr\tstart\tend\tstrand\tgene_id" > /path/to/${cohort}/QTL/result/p1e05/splicing_peer/${cell}_${tool}.txt
done
done



# make annotation


## remove variants with duplicate position
for chr in {1..22}; do
echo $chr
awk '++a[$4] == 2' tools/ldsc/GRCh38/plink_files/1000G.EUR.hg38.${chr}.bim | awk '{print "chr" $1 "\t" $4-1 "\t" $4+1 }'> pos_${chr}.bed
done
cat pos_*.bed | LANG=C sort -V -k 1,1 -k 2,2n > pos.bed

## variant ± 500 bp
awk '{print "chr" $1 "\t" $2-500  "\t" $2+500 }' ${HOME}/tools/DaPars2/QTL/FIVE_UTR-PDUI_1E05_annotated.vcf | grep -v "#" | LANG=C sort -V -k 1,1 -k 2,2n | bedtools merge | bedtools subtract -a - -b pos.bed > ${HOME}/tools/DaPars2/QTL/FIVE_UTR-PDUI_1E05_annotated.bed
awk '{print "chr" $1 "\t" $2-500  "\t" $2+500 }' ${HOME}/tools/DaPars2/QTL/THREE_UTR-PDUI_1E05_annotated.vcf | grep -v "#" | LANG=C sort -V -k 1,1 -k 2,2n | bedtools merge | bedtools subtract -a - -b pos.bed > ${HOME}/tools/DaPars2/QTL/THREE_UTR-PDUI_1E05_annotated.bed
awk '{print "chr" $1 "\t" $2-500  "\t" $2+500 }' ${HOME}/tools/DaPars2/QTL/eQTL_annotated.vcf | grep -v "#" | LANG=C sort -V -k 1,1 -k 2,2n | bedtools merge | bedtools subtract -a - -b pos.bed | uniq > ${HOME}/tools/DaPars2/QTL/eQTL_annotated.bed
awk '{print "chr" $1 "\t" $2-500  "\t" $2+500 }' ${HOME}/tools/DaPars2/QTL/sQTL_annotated.vcf | grep -v "#" | LANG=C sort -V -k 1,1 -k 2,2n | bedtools merge | bedtools subtract -a - -b pos.bed | uniq > ${HOME}/tools/DaPars2/QTL/sQTL_annotated.bed


for chr in {1..22}; do
for feature in FIVE_UTR-PDUI_1E05 THREE_UTR-PDUI_1E05 eQTL sQTL; do
mkdir --parents $HOME/tools/ldsc/GRCh38/ldsc/${feature}
echo "chr${chr}, ${feature}"
$HOME/miniconda3/envs/ldsc/bin/python $HOME/tools/ldsc/make_annot.py \
--bed-file $HOME/tools/DaPars2/QTL/${feature}_annotated.bed \
--bimfile $HOME/tools/ldsc/GRCh38/plink_files/1000G.EUR.hg38.${chr}.bim \
--annot-file $HOME/tools/ldsc/GRCh38/ldsc/${feature}/${chr}.annot.gz
done
done


for chr in {1..22}; do
for feature in FIVE_UTR-PDUI_1E05 THREE_UTR-PDUI_1E05 eQTL sQTL; do
$HOME/miniconda3/envs/ldsc/bin/python $HOME/tools/ldsc/ldsc.py \
--l2 \
--bfile $HOME/tools/ldsc/GRCh38/plink_files/1000G.EUR.hg38.${chr} \
--ld-wind-cm 1 \
--annot $HOME/tools/ldsc/GRCh38/ldsc/${feature}/${chr}.annot.gz \
--thin-annot \
--out $HOME/tools/ldsc/GRCh38/ldsc/${feature}/${chr} \
--print-snps $HOME/tools/ldsc/list.txt # https://github.com/bulik/ldsc/issues/233: The baseline annotations were generated using SNPs from 'list.txt' (https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/list.txt) instead of the hapmap3_snps
done
done


# make sumstats
mkdir --parents tools/ldsc/sumstats
for disease in COVID19_HGI_A2; do
echo $disease
LANG=C join -1 2 -2 2 <(LANG=C sort -k 2 /data02/home/juninamo/reference/KGP/AF/all_rs_autosome_EUR_GRCh38_MAF0.01_short.txt) <(LANG=C sort -k 2 /path/to/nanopore/RNA/210931_PB29_fastq/shiny/data/gwas/${disease}_GRCh38.txt) | awk '{print $2 "\t" $6 "\t" $7 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' | tr " " "\t" | sed -e "1i RSID\tA2\tA1\tBETA\tP\tMAF\tN\tZ" > tools/ldsc/sumstats/${disease}.txt
$HOME/miniconda3/envs/ldsc/bin/python $HOME/tools/ldsc/munge_sumstats.py \
--sumstats tools/ldsc/sumstats/${disease}.txt \
--merge-alleles tools/ldsc/w_hm3.snplist \
--out reference/sumstats_formatted/all_sumstats_alkes/${disease}
done

disease="COVID19_HGI_A2"
$HOME/miniconda3/envs/ldsc/bin/python $HOME/tools/ldsc/munge_sumstats.py \
--sumstats COVID19_HGI_A2.txt \
--merge-alleles tools/ldsc/w_hm3.snplist \
--out reference/sumstats_formatted/all_sumstats_alkes/${disease}

COVID19_HGI_A2 Ishigaki_RA Kottgen_Gout Lopez_SSc Patsopoulos_MS Rheenen_ALS Taylor_SjS UKB_460K.disease_AOSD UKB_460K.disease_AS_SELF_REP UKB_460K.disease_SpA_SELF_REP Zeggini_OA 

# Partitioned Heritability
for feature in FIVE_UTR-PDUI_1E05 THREE_UTR-PDUI_1E05 eQTL sQTL; do
for disease in PASS_ADHD_Demontis2018 PASS_AdultOnsetAsthma_Ferreira2019 PASS_AgeFirstBirth PASS_AgeOfInitiation_Liu2019 PASS_Alzheimer PASS_AlzheimersMaternal_Marioni2018 PASS_AlzheimersPaternal_Marioni2018 PASS_AlzheimersProxy_Marioni2018 PASS_Alzheimers_Jansen2019 PASS_Anorexia PASS_AnorexiaNervosa_Duncan_2017 PASS_AtrialFibrillation_Nielsen2018 PASS_Autism PASS_Autism_Grove2019 PASS_BDSCZ_Ruderfer2018 PASS_BIP_Stahl2019 PASS_BMI1 PASS_BipolarDisorder_Ruderfer2018 PASS_Bipolar_Disorder PASS_Bone_Mineral_Density_Kemp_2017 PASS_BreastCancer PASS_CD_deLange2017 PASS_CardioembolicStroke_Malik2018 PASS_Celiac PASS_ChildOnsetAsthma_Ferreira2019 PASS_CigarettesPerDay_Liu2019 PASS_Coronary_Artery_Disease PASS_Coronary_Artery_Disease_Howson_2017 PASS_Crohns_Disease PASS_DS PASS_DepressedAffect_Nagel2018 PASS_Depression_Nagel2018 PASS_DrinksPerWeek_Liu2019 PASS_ENIGMA2_MeanPutamen PASS_Epilepsy_Anney_2014 PASS_Ever_Smoked PASS_FastingGlucose_Manning PASS_FetalBirthWeight_Warrington2019 PASS_GeneralRiskTolerance_KarlssonLinner2019 PASS_HDL PASS_HbA1C PASS_Height1 PASS_IBD PASS_IBD_deLange2017 PASS_Insomnia_Jansen2019 PASS_Intelligence_SavageJansen2018 PASS_IschemicStroke_Malik2018 PASS_LDL PASS_LargeArteryStroke_Malik2018 PASS_Lean_Body_Mass_Zillikens_2017 PASS_LongSleepDuration_Dashti2019 PASS_Lupus PASS_MDD_Howard2019 PASS_MDD_Wray2018 PASS_MaternalBirthWeight_Warrington2019 PASS_MedicationUse_Wu2019 PASS_Menarche2017 PASS_Multiple_sclerosis PASS_Neuroticism PASS_Neuroticism_Nagel2018 PASS_NumberChildrenEverBorn PASS_OvarianCancer PASS_Parkinsons_23andMe_Corces2020 PASS_Primary_biliary_cirrhosis PASS_ProstateCancer PASS_ReactionTime_Davies2018 PASS_Rheumatoid_Arthritis PASS_SCZvsBD_Ruderfer2018 PASS_SWB PASS_Schizophrenia PASS_Schizophrenia_Pardinas2018 PASS_Schizophrenia_Ruderfer2018 PASS_ShortSleepDuration_Dashti2019 PASS_SleepDuration_Dashti2019 PASS_SmallVesselStroke_Malik2018 PASS_SmokingCessation_Liu2019 PASS_SmokingInitiation_Liu2019 PASS_Stroke_Malik2018 PASS_Triglycerides PASS_Type_1_Diabetes PASS_Type_2_Diabetes PASS_UC_deLange2017 PASS_Ulcerative_Colitis PASS_VerbalNumericReasoning_Davies2018 PASS_Worry_Nagel2018 PASS_Years_of_Education1 PASS_Years_of_Education2 UKB_460K.biochemistry_AlanineAminotransferase UKB_460K.biochemistry_Albumin UKB_460K.biochemistry_AlkalinePhosphatase UKB_460K.biochemistry_ApolipoproteinA UKB_460K.biochemistry_ApolipoproteinB UKB_460K.biochemistry_AspartateAminotransferase UKB_460K.biochemistry_Calcium UKB_460K.biochemistry_Cholesterol UKB_460K.biochemistry_CreactiveProtein UKB_460K.biochemistry_Creatinine UKB_460K.biochemistry_CystatinC UKB_460K.biochemistry_DirectBilirubin UKB_460K.biochemistry_GammaGlutamyltransferase UKB_460K.biochemistry_Glucose UKB_460K.biochemistry_HDLcholesterol UKB_460K.biochemistry_HbA1c UKB_460K.biochemistry_IGF1 UKB_460K.biochemistry_LDLdirect UKB_460K.biochemistry_LipoproteinA UKB_460K.biochemistry_Phosphate UKB_460K.biochemistry_SHBG UKB_460K.biochemistry_Testosterone_Male UKB_460K.biochemistry_TotalBilirubin UKB_460K.biochemistry_TotalProtein UKB_460K.biochemistry_Triglycerides UKB_460K.biochemistry_Urate UKB_460K.biochemistry_Urea UKB_460K.biochemistry_VitaminD UKB_460K.blood_EOSINOPHIL_COUNT UKB_460K.blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT UKB_460K.blood_LYMPHOCYTE_COUNT UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN UKB_460K.blood_MEAN_PLATELET_VOL UKB_460K.blood_MEAN_SPHERED_CELL_VOL UKB_460K.blood_MONOCYTE_COUNT UKB_460K.blood_PLATELET_COUNT UKB_460K.blood_PLATELET_DISTRIB_WIDTH UKB_460K.blood_RBC_DISTRIB_WIDTH UKB_460K.blood_RED_COUNT UKB_460K.blood_WHITE_COUNT UKB_460K.bmd_HEEL_TSCOREz UKB_460K.body_BALDING1 UKB_460K.body_BALDING4 UKB_460K.body_BMIz UKB_460K.body_HEIGHTz UKB_460K.body_LEFT_HANDED UKB_460K.body_WHRadjBMIz UKB_460K.bp_DIASTOLICadjMEDz UKB_460K.bp_SYSTOLICadjMEDz UKB_460K.cancer_ALL UKB_460K.cancer_BREAST UKB_460K.cancer_MELANOMA UKB_460K.cancer_PROSTATE UKB_460K.cov_EDU_COLLEGE UKB_460K.cov_EDU_YEARS UKB_460K.cov_SMOKING_STATUS UKB_460K.disease_AID_ALL UKB_460K.disease_AID_SURE UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED UKB_460K.disease_ASTHMA_DIAGNOSED UKB_460K.disease_CARDIOVASCULAR UKB_460K.disease_DERMATOLOGY UKB_460K.disease_DIABETES_ANY_DIAGNOSED UKB_460K.disease_ENDOCRINE_DIABETES UKB_460K.disease_HI_CHOL_SELF_REP UKB_460K.disease_HYPERTENSION_DIAGNOSED UKB_460K.disease_HYPOTHYROIDISM_SELF_REP UKB_460K.disease_PSORIASIS UKB_460K.disease_RESPIRATORY_ENT UKB_460K.disease_T2D UKB_460K.disease_THYROID_ANY_SELF_REP UKB_460K.impedance_BASAL_METABOLIC_RATEz UKB_460K.lung_FEV1FVCzSMOKE UKB_460K.lung_FVCzSMOKE UKB_460K.mental_NEUROTICISM UKB_460K.other_MORNINGPERSON UKB_460K.pigment_HAIR UKB_460K.pigment_HAIR_blackmale UKB_460K.pigment_HAIR_blonde UKB_460K.pigment_HAIR_darkbrown UKB_460K.pigment_SKIN UKB_460K.pigment_SUNBURN UKB_460K.pigment_TANNING UKB_460K.repro_AgeFirstBirth_Female UKB_460K.repro_MENARCHE_AGE UKB_460K.repro_MENOPAUSE_AGE UKB_460K.repro_NumberChildrenEverBorn_Female UKB_460K.repro_NumberChildrenEverBorn_Male UKB_460K.repro_NumberChildrenEverBorn_Pooled; do
echo $disease
$HOME/miniconda3/envs/ldsc/bin/python $HOME/tools/ldsc/ldsc.py \
--h2 reference/sumstats_formatted/all_sumstats_alkes/${disease}.sumstats.gz \
--ref-ld-chr $HOME/tools/ldsc/GRCh38/ldsc/${feature}/,$HOME/tools/ldsc/GRCh38/baseline_v1.2/ \
--w-ld-chr $HOME/tools/ldsc/GRCh38/weights/ \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--frqfile-chr $HOME/tools/ldsc/EUR/1000G_Phase3_frq/1000G.EUR.QC. \
--out $HOME/tools/ldsc/results/Partitioned_Heritability/${disease}_${feature}
done
done

for chr in {1..22}; do
cp $HOME/tools/ldsc/GRCh38/baseline_v1.2/baseline.${chr}.l2.ldscore.gz $HOME/tools/ldsc/GRCh38/baseline_v1.2/${chr}.l2.ldscore.gz
cp $HOME/tools/ldsc/GRCh38/baseline_v1.2/baseline.${chr}.l2.M_5_50 $HOME/tools/ldsc/GRCh38/baseline_v1.2/${chr}.l2.M_5_50
cp $HOME/tools/ldsc/GRCh38/baseline_v1.2/baseline.${chr}.l2.M $HOME/tools/ldsc/GRCh38/baseline_v1.2/${chr}.l2.M
cp $HOME/tools/ldsc/GRCh38/baseline_v1.2/baseline.${chr}.annot.gz $HOME/tools/ldsc/GRCh38/baseline_v1.2/${chr}.annot.gz
done


for chr in {1..22}; do
echo "chr${chr}"
zcat $HOME/tools/ldsc/GRCh38/baseline_v1.2/${chr}.l2.ldscore.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/FIVE_UTR-PDUI_1E05/${chr}.l2.ldscore.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/THREE_UTR-PDUI_1E05/${chr}.l2.ldscore.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/eQTL/${chr}.l2.ldscore.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/sQTL/${chr}.l2.ldscore.gz | wc -l
echo -e "\n"
wc -l $HOME/tools/ldsc/GRCh38/plink_files/1000G.EUR.hg38.${chr}.bim
zcat $HOME/tools/ldsc/GRCh38/baseline_v1.2/${chr}.annot.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/FIVE_UTR-PDUI_1E05/${chr}.annot.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/THREE_UTR-PDUI_1E05/${chr}.annot.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/eQTL/${chr}.annot.gz | wc -l
zcat $HOME/tools/ldsc/GRCh38/ldsc/sQTL/${chr}.annot.gz | wc -l
echo -e "\n"
done
