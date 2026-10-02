```bash

time ./ADVTIG-untargeted/split_fa.py
BASE="/home/james/SequencingData/"
find /home/james/SequencingData/NCBI_RVDB/virus_RVDBs/URVDB_split_fa/MH5904* -type f -iname "*Human_gammaherpesvirus_4*" | head

SPLIT_SAVE="${BASE}NCBI_RVDB/virus_RVDBs/URVDB_split_fa/"
INDEX="${BASE}NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix2-index.json"

# used the original centrifuge index
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.txt"
DIRECTORY="/mnt/usersData/ADVTIG_v2_untargeted/"
CENT_SAVE="~/SequencingData/Centrifuge_libraries/viral/"
CE_DICT_ID="uviral25"
time ~/SMART-CAMP/mp_metagenomic_assessment_v4.py -d "${DIRECTORY}" -t 5 -ci "${CENT_SAVE}" -c -bl "${CE_DICT_ID}" -m -rs 1 -fd "${SAMPLES}" -skip

# use the new db
cd
cd SMART-CAMP
CE_DICT_ID="uviral25-2"
DIRECTORY="/mnt/usersData/ADVTIG_v2_untargeted/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.txt"
CENT_SAVE="~/SequencingData/Centrifuge_libraries/viral/"
time ~/SMART-CAMP/mp_metagenomic_assessment_v4.py -d "${DIRECTORY}" -t 5 -ci "${CENT_SAVE}" -c -bl "${CE_DICT_ID}" -m -rs 1 -fd "${SAMPLES}" -skip

awk 'BEGIN{FS="\t"} $2 == 11768 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # felv
awk 'BEGIN{FS="\t"} $2 == 61673 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # perv
awk 'BEGIN{FS="\t"} $2 == 133704 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # pcv1
awk 'BEGIN{FS="\t"} $2 == 11856 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # smrv
awk 'BEGIN{FS="\t"} $2 == 28285 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # adv5
awk 'BEGIN{FS="\t"} $2 == 10376 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # ebv/hh4
awk 'BEGIN{FS="\t"} $2 == 12814 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # rsv
awk 'BEGIN{FS="\t"} $2 == 11250 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # Human orthopneumovirus
awk 'BEGIN{FS="\t"} $2 == 10891 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_mini-u-viral2_centrifuge_report.tsv # Reovirus

# DIRECTORY="/mnt/usersData/ADVTIG_v2_untargeted/"
# SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
# CE_DICT_ID="uviral25-2"
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 -e # asm5; extract fastq [HISTORICAL command, preserved as originally run: pre-1c782e3e run_coverage CLI where -e = extract FASTQ (current: -x); -c removed in aac5714b]
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 -r # asm5; run analysis
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 # asm5; stats and concatenate
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 -l # add labels to the mapped csv

time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/Viral_FDA/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all16-2.txt" -t lambda -p 6
time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/ADVTIG_v2_untargeted/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all9-2.txt" -t lambda -p 6

# first run
# 1X coverage from the pileup file
# percent identity from the consensus pileup file
# pileup is from raw alignment (pre-filtering)

# second run
# pileup is regenerated from the mapped reads.

DIRECTORY="/mnt/usersData/ADVTIG_v2_untargeted/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
CE_DICT_ID="uviral25-2"
# "First Run on unfiltered reads"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 -x # map-ont; model 2; extract fastq - 2400m [current CLI 2026-10-02: --edit-distance 0.15 -x; as committed in 1c782e3e this line read '-e 0.15 -e', which never parsed]
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 -r # map-ont; model 2; run analysis - 1200m [current CLI 2026-10-02: --edit-distance 0.15; as committed in 1c782e3e this line read '-e 0.15 -r', which never parsed]
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 # map-ont; model 2; stats and concatenate - 1000m [current CLI 2026-10-02: --edit-distance 0.15; as committed in 1c782e3e this line read '-e 0.15', which never parsed]
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 -l # add labels to the mapped csv [current CLI 2026-10-02: --edit-distance 0.15; as committed in 1c782e3e this line read '-e 0.15 -l', which never parsed]

# "Second Run on filtered reads"
DIRECTORY="/mnt/usersData/ADVTIG_v3_untargeted/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
CE_DICT_ID="uviral25-2"
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e 0.15 -dd # deduplicate fq - 1124m
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e 0.15 -dd -x # extract stats - 1024m (rerun) # rerun for intermediate json if changing ED
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e 0.15 -dd -a # extract stats - 34m

time ./ADVTIG-untargeted/filt_low_complexity.py -o "${DIRECTORY}"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"

time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES}" -p "15" -c -e 0.15 # 15% edit threshold
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES}" -p "15" -c --panel "new_combined" -e 0.15 # 15% edit threshold# 10m

time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES}" -p "15" -c --panel "refseq" -e 0.15 # 15% edit threshold# 1m # cat fq vs panel
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES}" -p "15" --panel "refseq" -r -e 0.15 # 15% edit threshold# m # subset fq vs panel
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES}" -p "15" -c --panel "bacterial_refseq" -e 0.15 # 15% edit threshold# 1m # cat fq vs panel
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES}" -p "15" --panel "bacterial_refseq" -r -e 0.15 # 15% edit threshold# m # subset fq vs panel

time ./ADVTIG-untargeted/count_the_screen.py  -o "${DIRECTORY}" -s "${SAMPLES}" -e 0.15
time ./ADVTIG-untargeted/count_the_screen.py  -o "${DIRECTORY}" -s "${SAMPLES}" -e 0.15 --skip

# run 3X
time ./ADVTIG-untargeted/collect_taxid_primary_counts.py
time ./ADVTIG-untargeted/collect_taxid_primary_counts.py
time ./ADVTIG-untargeted/collect_taxid_primary_counts.py

# time ./ADVTIG-untargeted/deep_cov.py # might need running twice to generate Third-Pass-FULL-untargeted-TaxID_taxID-uviral25-2_mapping_stats_labelled_filtered.csv

#########################

# Viral_FDA
DIRECTORY="/mnt/usersData/Viral_FDA/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all16-2.txt"
CE_DICT_ID="uviral25-2"
CENT_SAVE="~/SequencingData/Centrifuge_libraries/viral/"

# Untargeted metagenomics analysis
time ~/SMART-CAMP/mp_metagenomic_assessment_v4.py -d "${DIRECTORY}" -t 5 -ci "${CENT_SAVE}" -c -bl "${CE_DICT_ID}" -m -rs 1 -fd "${SAMPLES}" -skip

# "First Run on unfiltered reads"
SAMPLES_CSV="/home/james/SMART-CAMP/configs/viral_DNA_all16-2.csv"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 -x # map-ont; model 2; extract fastq - 2400m [current CLI 2026-10-02: --edit-distance 0.15 -x; as committed in 1c782e3e this line read '-e 0.15 -e', which never parsed]
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 -r # map-ont; model 2; run analysis - 1200m [current CLI 2026-10-02: --edit-distance 0.15; as committed in 1c782e3e this line read '-e 0.15 -r', which never parsed]
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 # map-ont; model 2; stats and concatenate - 1000m [current CLI 2026-10-02: --edit-distance 0.15; as committed in 1c782e3e this line read '-e 0.15', which never parsed]
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" --edit-distance 0.15 -l # add labels to the mapped csv [current CLI 2026-10-02: --edit-distance 0.15; as committed in 1c782e3e this line read '-e 0.15 -l', which never parsed]

# "Second Run on filtered reads"
DIRECTORY="/mnt/usersData/Viral_FDA/"
SAMPLES_CSV="/home/james/SMART-CAMP/configs/viral_DNA_all16-2.csv"
CE_DICT_ID="uviral25-2"
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e 0.15 -dd # deduplicate fq - 1124m
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e 0.15 -dd -x # extract stats - 1024m (rerun) # rerun for intermediate json if changing ED
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e 0.15 -dd -a # extract stats - 34m

# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e # map-ont; model 2; extract fastq - 1800m [HISTORICAL command, preserved as originally run: pre-1c782e3e run_coverage CLI where -e = extract FASTQ (current: -x)]
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -r # map-ont; model 2; run analysis - 900m
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" # map-ont; model 2; stats and concatenate - 700m
# time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -l # add labels to the mapped csv - 292m

# # run first twice
# time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -dd # deduplicate fq - 159m
# time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -dd -e # extract stats [HISTORICAL command, preserved as originally run: pre-1c782e3e deep_cov CLI where -e = extract (current: -x)]
# time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES_CSV}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -dd -a # extract stats

time ./ADVTIG-untargeted/filt_low_complexity.py -o "${DIRECTORY}"

time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -p "15" -c -e 0.15 # 15% edit threshold
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -p "15" -c --panel "new_combined" -e 0.15 # 15% edit threshold# 10m
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -p "15" -c --panel "refseq" -e 0.15 # 15% edit threshold# 1m # cat fq vs panel
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -p "15" --panel "refseq" -r -e 0.15 # 15% edit threshold# m # subset fq vs panel

time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -p "15" -c --panel "bacterial_refseq" -e 0.15 # 15% edit threshold# 1m # cat fq vs panel
time ./ADVTIG-untargeted/counter_screen.py -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -p "15" --panel "bacterial_refseq" -r -e 0.15 # 15% edit threshold# m # subset fq vs panel

time ./ADVTIG-untargeted/count_the_screen.py  -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -e 0.15
time ./ADVTIG-untargeted/count_the_screen.py  -o "${DIRECTORY}" -s "${SAMPLES_CSV}" -e 0.15 --skip



/mnt/usersData/Viral_FDA//map-ont/uviral25-2_mapping_stats_labelled.csv
/mnt/usersData/ADVTIG_v2_untargeted//map-ont/uviral25-2_mapping_stats_labelled.csv

rm /mnt/usersData/Viral_FDA//map-ont/2*/**/*sam 
rm /mnt/usersData/ADVTIG_v2_untargeted//map-ont/2*/**/*sam 

rm /mnt/usersData/Viral_FDA//map-ont/2*/**/*bam 
rm /mnt/usersData/ADVTIG_v2_untargeted//map-ont/2*/**/*bam 


cat \
Cat_GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fa \
Human_GCF_000001405.40_GRCh38.p14_genomic.fa \
Monkey_GCF_003339765.1_Mmul_10_genomic.fa \
Pig_GCF_000003025.6_Sscrofa11.1_genomic.fa \
GCA_021498475.1_ASM2149847v1/ncbi_dataset/data/GCA_021498475.1//GCA_021498475.1_ASM2149847v1_genomic.fna \
> all_mammalian.fa

cat \
Cat_GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fa \
Human_GCF_000001405.40_GRCh38.p14_genomic.fa \
Monkey_GCF_003339765.1_Mmul_10_genomic.fa \
Pig_GCF_000003025.6_Sscrofa11.1_genomic.fa \
GCA_021498475.1_ASM2149847v1/ncbi_dataset/data/GCA_021498475.1//GCA_021498475.1_ASM2149847v1_genomic.fna \
RSV-A.fa \
'REOvirus segments.fa' \
PCV1.fa \
NC_001510.1.fna \
EBV-B95-8.fa \
FeLV.fa \
ADV5.fa \
NC_001514.1.fna \
FeLV_Kawakami-Theilen_strain.fa \
NC_005831.2_HCoV-NL63.fa \
AF038600.1.fna \
> all_mammalian_targets.fa

cat all_mammalian_targets.fa AF038600.1.fna /home/james/SequencingData/Centrifuge_libraries/ecoli-lambda-phage.NC_001416.1.fasta > all_mammalian_targets2.fa

minimap2 -d all_mammalian_targets2.mmi all_mammalian_targets2.fa

grep ">" /mnt/usersData/ADVTIG_fa/Cat_GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fa \
  | awk '{print $1"\tCat"}' \
  > /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/Human_GCF_000001405.40_GRCh38.p14_genomic.fa \
  | awk '{print $1"\tHuman"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/Monkey_GCF_003339765.1_Mmul_10_genomic.fa \
  | awk '{print $1"\tMonkey"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/Pig_GCF_000003025.6_Sscrofa11.1_genomic.fa \
  | awk '{print $1"\tPig"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/GCA_021498475.1_ASM2149847v1/ncbi_dataset/data/GCA_021498475.1/GCA_021498475.1_ASM2149847v1_genomic.fna \
  | awk '{print $1"\tTamarin"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/RSV-A.fa \
  | awk '{print $1"\tRSV_A"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" "/mnt/usersData/ADVTIG_fa/REOvirus segments.fa" \
  | awk '{print $1"\tREOvirus"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/PCV1.fa \
  | awk '{print $1"\tPCV1"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/NC_001510.1.fna \
  | awk '{print $1"\tNC_001510_1"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/EBV-B95-8.fa \
  | awk '{print $1"\tEBV_B95_8"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/FeLV.fa \
  | awk '{print $1"\tFeLV"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/ADV5.fa \
  | awk '{print $1"\tADV5"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/NC_001514.1.fna \
  | awk '{print $1"\tNC_001514_1"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/FeLV_Kawakami-Theilen_strain.fa \
  | awk '{print $1"\tFeLV_Kawakami_Theilen"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv && \
grep ">" /mnt/usersData/ADVTIG_fa/NC_005831.2_HCoV-NL63.fa \
  | awk '{print $1"\tHCoV_NL63"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv

cp /mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv \
   /mnt/usersData/ADVTIG_fa/mammalian_ref_species2.tsv

grep '^>' /mnt/usersData/ADVTIG_fa/AF038600.1.fna \
  | awk '{print $1"\tAF038600_1"}' \
  >> /mnt/usersData/ADVTIG_fa/mammalian_ref_species2.tsv

########

cd /home/james/SequencingData/Centrifuge_libraries/
minimap2 -d ecoli-lambda-phage.NC_001416.1.mmi ecoli-lambda-phage.NC_001416.1.fasta

in=/home/james/SequencingData/Centrifuge_libraries/ecoli-lambda-phage.fasta
out=/home/james/SequencingData/Centrifuge_libraries/ecoli-lambda-phage.NC_001416.1.fasta

sed '1s/^>.*/>NC_001416.1 Enterobacteria_phage_lambda_complete_genome/' "$in" > "$out"

in="/mnt/usersData/ADVTIG_fa/mammalian_ref_species.tsv"
out="/mnt/usersData/ADVTIG_fa/mammalian_ref_species2.tsv"

cp "$in" "$out"
echo -e "NC_001416.1\tDNA_CS" >> "$out"

########

FASTA="/home/james/SequencingData/Centrifuge_libraries/viral/viral_input-sequences.fna"
MMI="/mnt/usersData/ADVTIG_fa/viral-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/viral-input-sequences.tsv"
minimap2 -d ${MMI} ${FASTA}

grep '^>' "$FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"

FASTA="/home/james/SequencingData/Centrifuge_libraries/v_f_b/fungal-seq-2022.fna"
MMI="/mnt/usersData/ADVTIG_fa/fungal-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/fungal-input-sequences.tsv"
minimap2 -d ${MMI} ${FASTA}

grep '^>' "$FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"

FASTA="/home/james/SequencingData/Centrifuge_libraries/v_f_b/ss_bacteria_sequences.fna"
MMI="/mnt/usersData/ADVTIG_fa/bacterial-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/bacterial-input-sequences.tsv"
minimap2 -d ${MMI} ${FASTA}

grep '^>' "$FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"

FASTA="/mnt/usersData/SILVA138/SILVA_138.1_LSURef_NR99_tax_silva.fasta"
MMI="/mnt/usersData/ADVTIG_fa/LSU-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/LSU-input-sequences.tsv"
minimap2 -d ${MMI} ${FASTA}

grep '^>' "$FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"

FASTA="/mnt/usersData/SILVA138/library/SILVA_138.1_SSURef_NR99_tax_silva.clean.fasta"
MMI="/mnt/usersData/ADVTIG_fa/SSU-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/SSU-input-sequences.tsv"
minimap2 -d ${MMI} ${FASTA}

grep '^>' "$FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"


FASTA="/mnt/usersData/v-f-b/NC_009656.1_Pseudomonas-aeruginosa-PA7--complete-gen.fa"
MMI="/mnt/usersData/ADVTIG_fa/NC_009656.1-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/NC_009656.1-input-sequences.tsv"
minimap2 -d ${MMI} ${FASTA}

grep '^>' "$FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"

###########

FASTA="/home/james/SequencingData/Centrifuge_libraries/v_f_b/ss_bacteria_sequences.fna"
PA_FASTA="/mnt/usersData/ADVTIG_fa/bacterial-Pa-input-sequences.fasta"
MMI="/mnt/usersData/ADVTIG_fa/bacterial-Pa-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/bacterial-Pa-input-sequences.tsv"

awk 'BEGIN{IGNORECASE=1}
     /^>/ {keep = ($0 ~ /Pseudomonas([;[:space:]]+|_)+aeruginosa/)}
     keep {print}' "$FASTA" > "$PA_FASTA"

minimap2 -d ${MMI} ${PA_FASTA}

grep '^>' "$PA_FASTA" \
  | sed 's/^>//' \
  | awk '{
      acc = $1;           # first token = accession
      $1 = "";            # drop accession from the line
      sub(/^ +/, "");     # trim leading spaces
      desc = $0;
      gsub(/[ ,]/, "_", desc);       # spaces/commas -> _
      gsub(/[^A-Za-z0-9_]/, "", desc); # strip quotes etc.
      print acc "\t" desc;
  }' > "$TSV"

############

FASTA="/mnt/usersData/SILVA138/SILVA_138.1_LSURef_NR99_tax_silva.fasta"
PA_FASTA="/mnt/usersData/ADVTIG_fa/LSU-Pa-input-sequences.fasta"
MMI="/mnt/usersData/ADVTIG_fa/LSU-Pa-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/LSU-Pa-input-sequences.tsv"

# 1) subset FASTA (keep full records)
awk 'BEGIN{IGNORECASE=1}
     /^>/ {keep = ($0 ~ /Pseudomonas([;[:space:]]+|_)+aeruginosa/)}
     keep {print}' "$FASTA" > "$PA_FASTA"

minimap2 -d "${MMI}" "${PA_FASTA}"

grep '^>' "$PA_FASTA" \
  | awk '{
      line=$0
      sub(/^>/,"",line)

      # accession = first whitespace-delimited token
      acc=line
      sub(/[[:space:]].*$/,"",acc)

      # desc = rest after accession + spaces
      desc=line
      sub(/^[^[:space:]]+[[:space:]]+/,"",desc)

      # separators -> underscore
      gsub(/[;[:space:],]+/, "_", desc)

      # strip anything else (keep only A-Za-z0-9 and _)
      gsub(/[^A-Za-z0-9_]/, "", desc)

      print acc "\t" desc
  }' > "$TSV"

############

FASTA="/mnt/usersData/SILVA138/library/SILVA_138.1_SSURef_NR99_tax_silva.clean.fasta"
PA_FASTA="/mnt/usersData/ADVTIG_fa/SSU-Pa-input-sequences.fasta"
MMI="/mnt/usersData/ADVTIG_fa/SSU-Pa-input-sequences.mmi"
TSV="/mnt/usersData/ADVTIG_fa/SSU-Pa-input-sequences.tsv"

# 1) subset FASTA (keep full records)
awk 'BEGIN{IGNORECASE=1}
     /^>/ {keep = ($0 ~ /Pseudomonas[[:space:]]+aeruginosa/)}
     keep {print}' "$FASTA" > "$PA_FASTA"

minimap2 -d "${MMI}" "${PA_FASTA}"

grep '^>' "$FASTA" \
  | awk -F'\t' 'BEGIN{IGNORECASE=1}
    /Pseudomonas[[:space:]]+aeruginosa/ {
      line=$0
      sub(/^>/,"",line)                         # drop >
      split(line, a, /[[:space:]]+/)            # split on whitespace
      acc=a[1]                                  # accession.token
      desc=substr(line, length(acc)+2)          # rest of header after space
      gsub(/[ ,]/, "_", desc)
      gsub(/[^A-Za-z0-9_;]/, "", desc)          # keep alnum/_/; (optional)
      gsub(/;/, "_", desc)                      # if you want ; -> _
      print acc "\t" desc
    }' > "$TSV"
############

fq4="/mnt/usersData/DNA/analysis/sample_data/20251024_aDNA_TC-Paeruginosa-16S-23S-LSK114_3CFU_9125_12/trimmed/trimmed_20251024_aDNA_TC-Paeruginosa-16S-23S-LSK114_3CFU_9125_12.Pseudomonas-aeruginosa.fastq"
fq5="/mnt/usersData/methylfilter_5mC/ReadLengths/20251024_aDNA_TC-Paeruginosa-16S-23S-LSK114_3CFU_9125_12.Pseudomonas-aeruginosa.0-5000bp.fastq"
cp ${fq4} ${fq5}

fq="/mnt/usersData/methylfilter_5mC/ReadLengths/20251028_aDNA_TC-Paeruginosa-16S-23S-5mC-M.SssI-LSK114_3CFU_9127_12.Pseudomonas-aeruginosa.2450-2650bp.fastq"
fq2="/mnt/usersData/methylfilter_5mC/ReadLengths/20251028_aDNA_TC-Paeruginosa-16S-23S-5mC-M.SssI-LSK114_3CFU_9127_12.Pseudomonas-aeruginosa.0-5000bp.fastq"
PAF="/mnt/usersData/methylfilter_5mC/ReadLengths/20251028_aDNA_TC-Paeruginosa-16S-23S-5mC-M.SssI-LSK114_3CFU_9127_12.Pseudomonas-aeruginosa.2450-2650bp.paf"
PAF2="/mnt/usersData/methylfilter_5mC/ReadLengths/20251028_aDNA_TC-Paeruginosa-16S-23S-5mC-M.SssI-LSK114_3CFU_9127_12.Pseudomonas-aeruginosa.0-5000bp.paf"

fq3="/mnt/usersData/DNA/analysis/sample_data/20251028_aDNA_TC-Paeruginosa-16S-23S-5mC-M.SssI-LSK114_3CFU_9127_12/trimmed/trimmed_20251028_aDNA_TC-Paeruginosa-16S-23S-5mC-M.SssI-LSK114_3CFU_9127_12.Pseudomonas-aeruginosa.fastq"
cp ${fq3} ${fq2}

minimap2 -cx map-ont -t 4 ${MMI} ${fq} > ${PAF}
minimap2 -cx map-ont -t 4 ${MMI} ${fq2} > ${PAF2}

# 1) header/plus sanity
awk 'NR%4==1 && $0 !~ /^@/ {print "BAD_HEADER",NR,$0; exit}
     NR%4==3 && $0 !~ /^\+/ {print "BAD_PLUS",NR,$0; exit}' "$fq"

# 2) seq/qual length match
awk 'NR%4==2{s=length($0)}
     NR%4==0{q=length($0); if(s!=q){print "LEN_MISMATCH read",NR/4,"seq",s,"qual",q; exit}}' "$fq"

# 3) non-ACGTN fraction (if huge → garbage/contaminated)
awk 'NR%4==2{c+=length($0); x=$0; gsub(/[ACGTN]/,"",x); b+=length(x)}
     END{print "bases",c,"non_ACGTN",b,"frac", (c?b/c:0)}' "$fq"

FASTA="/mnt/usersData/v-f-b/NC_009656.1_Pseudomonas-aeruginosa-PA7--complete-gen.fa"
MMI="/mnt/usersData/ADVTIG_fa/NC_009656.1-k19-w15-input-sequences.mmi"
minimap2 -k19 -w15 -d ${MMI} ${FASTA}
minimap2 -cx map-ont -t 4 -k19 ${MMI} "${fq}" 

minimap2 -cx map-pb -t 4 ${MMI} "${fq}" | head
minimap2 -cx map-pb -t 4 ${MMI} "${fq}" | wc -l

minimap2 -cx map-ont -t 4 -N 20 --secondary=yes ${MMI} "${fq}" | head
minimap2 -cx map-ont -t 4 -N 20 --secondary=yes ${MMI} "${fq}" | wc -l


```
