#!/bin/bash
##$ -S /bin/sh

# hmmer
## https://kazumaxneo.hatenablog.com/entry/2017/07/31/114955
## https://ncrna.jp/blog/item/450-drawing-domain-structures

# generate index
cd /path/to/Pfam35.0/
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gzip -d Pfam-A.hmm.gz
hmmpress /path/to/Pfam35.0/Pfam-A.hmm


# SQANTI3 output
# generate amino acid sequence fasta from SQANTI3 output
sed -e "1d" /path/to/GRCh38_gencode38_classification.filtered_lite_classification_modified.txt \
| awk '{print ">" $1 " " $2 " strand=" $3 " isoform-length=" $4 " structural_category=" $6 " associated_gene=" $7 " associated_transcript=" $8 " combination_of_known_junctions=" $15 " coding=" $30 " ORF_length=" $31 " CDS_length=" $32 " CDS_start=" $33 " CDS_end=" $34 " CDS_genomic_start=" $35 " CDS_genomic_end=" $36 " predicted_NMD=" $37 "\n" $45}' \
| grep -w -v "NA" \
> /path/to/GRCh38_gencode38_classification.filtered_lite_aa.fasta

cd /path/to/Pfam35.0/
hmmscan \
--domtblout /path/to/hmmscan_PfamA_output.txt \
--cpu 35 \
-E 1e-10 \
/path/to/Pfam35.0/Pfam-A.hmm \
/path/to/GRCh38_gencode38_classification.filtered_lite_aa.fasta


# http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
## (1) target name: The name of the target sequence or profile.
## (2) target accession: Accession of the target sequence or profile, or ’-’ if none is available.
## (3) tlen: Length of the target sequence or profile, in residues. This (together with the query length) is useful for interpreting where the domain coordinates (in subsequent columns) lie in the sequence.
## (4) query name: Name of the query sequence or profile.
## (5) accession: Accession of the target sequence or profile, or ’-’ if none is available.
## (6) qlen: Length of the query sequence or profile, in residues.
## (7) E-value: E-value of the overall sequence/profile comparison (including all domains).
## (8) score: Bit score of the overall sequence/profile comparison (including all domains), inclusive of a null2 bias composition correction to the score.
## (9) bias: The biased composition score correction that was applied to the bit score.
## (10) #: This domain’s number (1..ndom).
## (11) of: The total number of domains reported in the sequence, ndom.
## (12) c-Evalue: The “conditional E-value”, a permissive measure of how reliable this particular domain
## may be. The conditional E-value is calculated on a smaller search space than the independent Evalue. The conditional E-value uses the number of targets that pass the reporting thresholds. The
## null hypothesis test posed by the conditional E-value is as follows. Suppose that we believe that
## there is already sufficient evidence (from other domains) to identify the set of reported sequences as
## homologs of our query; now, how many additional domains would we expect to find with at least this
## particular domain’s bit score, if the rest of those reported sequences were random nonhomologous
## sequence (i.e. outside the other domain(s) that were sufficient to identified them as homologs in the
##           first place)?
## (13) i-Evalue: The “independent E-value”, the E-value that the sequence/profile comparison would have
## received if this were the only domain envelope found in it, excluding any others. This is a stringent
## measure of how reliable this particular domain may be. The independent E-value uses the total
## number of targets in the target database.
## (14) score: The bit score for this domain.
## (15) bias: The biased composition (null2) score correction that was applied to the domain bit score.
## (16) from (hmm coord): The start of the MEA alignment of this domain with respect to the profile, numbered 1..N for a profile of N consensus positions.
## (17) to (hmm coord): The end of the MEA alignment of this domain with respect to the profile, numbered
## 1..N for a profile of N consensus positions.
## (18) from (ali coord): The start of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
## (19) to (ali coord): The end of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
## (20) from (env coord): The start of the domain envelope on the sequence, numbered 1..L for a sequence of L residues. The envelope defines a subsequence for which their is substantial probability
## mass supporting a homologous domain, whether or not a single discrete alignment can be identified.
## The envelope may extend beyond the endpoints of the MEA alignment, and in fact often does, for
## weakly scoring domains.
## (21) to (env coord): The end of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
## (22) acc: The mean posterior probability of aligned residues in the MEA alignment; a measure of how reliable the overall alignment is (from 0 to 1, with 1.00 indicating a completely reliable alignment according to the model).
## (23) description of target: The remainder of the line is the target’s description line, as free text.
