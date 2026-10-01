```bash
###############################################################################################

# BASE="/home/james/SequencingData/"

# # FA="${BASE}NCBI_RVDB/virus_RVDBs/abrdiged-U.fasta"
# FA="${BASE}NCBI_RVDB/virus_RVDBs/U-RVDB-clean.fasta"


# # MAP="${BASE}Centrifuge_libraries/viral/u-mini.map"
# # awk -F'|' '/^>/{print $3 "\t" 100000 + NR}' ${FA} > ${MAP} # random seqID

# # use the taxID map from the RVDB2
# # MAP=/mnt/usersData/RVDB2/nucl_gb.accession2taxid
# # ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

# # return only 2 columns
awk 'NR > 1 {print $2 "\t" $3}' /mnt/usersData/RVDB2/nucl_gb.accession2taxid > /mnt/usersData/RVDB2/RVDB2_seqID.map

# FULLFA="${BASE}NCBI_RVDB/virus_RVDBs/U-RVDBvCurrent.fasta"
# SAVESEQID="${BASE}NCBI_RVDB/virus_RVDBs/URVDB_seqID.map"
# awk -F'|' 'BEGIN {print "code,name"} /^>/ {code=$3; $1=$2=$3=$4=""; gsub(/^[\ \t]+|[\|]+[^|]*\|[^|]*$/, ""); print code "," $0}' "${FULLFA}" > "${SAVESEQID}"
# head "${SAVESEQID}"

MAP=/mnt/usersData/RVDB2/RVDB2_seqID.map
head ${MAP}

# clean up fasta file removing fluff and retaining accession number
NEWFA="${BASE}NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix.fasta"
sed '/^>/s/^>[^|]*|[^|]*|\([^|]*\)|.*/>\1/' "${FA}" > "${NEWFA}"
head ${NEWFA}

NODE="${BASE}Centrifuge_libraries/viral/taxonomy/nodes.dmp"
NAME="${BASE}Centrifuge_libraries/viral/taxonomy/names.dmp"
OUT="mini-u-viral"

# # build a functioning U-VIRAL DB
cd ~/SequencingData/Centrifuge_libraries/viral/
centrifuge-build --conversion-table ${MAP} \
                 --taxonomy-tree ${NODE} --name-table ${NAME} \
                 ${NEWFA} ${OUT} -p 8

###############################################################################################
###############################################################################################
###############################################################################################
###################################> CLEAN UP <################################################
###############################################################################################

BASE="/home/james/SequencingData/"
FA="${BASE}NCBI_RVDB/virus_RVDBs/U-RVDBvCurrent.fasta"
NEWFA="${BASE}NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix2.fasta"

awk '/^>/ {
    split($0, fields, "|");
    if (fields[3] != "") {
        print ">" fields[3];
    } else {
        print $0;
    }
    next;
}
!/^>/ { print $0 }' "${FA}" > "${NEWFA}"
head "${NEWFA}"
rg -i ">acc" "${NEWFA}" | head

# note, there are some seqIDs missing from /home/james/SequencingData/NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map
# they can be found here:
# rg -i "KX868466.2" /home/james/SequencingData/NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix2-index.json
# they might need to be tweaked if CRVDB v25 is used again

# remake the seqID.map
SAVESEQID="${BASE}NCBI_RVDB/virus_RVDBs/URVDB_seqID2.map"
awk -F'|' 'BEGIN {print "code,name"}
/^>/ {
    code=$3;
    $1=$2=$3="";
    gsub(/^[ \t]+|[ \t]+$/, "");  # Remove leading/trailing spaces
    gsub(/\|$/, "");              # Remove trailing pipes
    print code "," $0;
}' "${FA}" > "${SAVESEQID}"
head "${SAVESEQID}"
rg -i "NC_001798.2" "${SAVESEQID}" | head

# clean up seqID.map
awk -F',' '{
    if (NF > 2) {
        $2 = $2 "." substr($3, 2);  # Replace the first extra comma with a period
        $3 = "";  # Clear the extra field created
    }
    print $1 "," $2;
}' OFS=',' "${BASE}NCBI_RVDB/virus_RVDBs/URVDB_seqID2.map" > "${BASE}NCBI_RVDB/virus_RVDBs/fixed_URVDB_seqID2.map"
sed -n '37908p' "${BASE}NCBI_RVDB/virus_RVDBs/fixed_URVDB_seqID2.map"
rg -i "LC012610.1" "${BASE}NCBI_RVDB/virus_RVDBs/fixed_URVDB_seqID2.map"

NODE="${BASE}Centrifuge_libraries/viral/taxonomy/nodes.dmp"
NAME="${BASE}Centrifuge_libraries/viral/taxonomy/names.dmp"
OUT="mini-u-viral2"

# build a functioning U-VIRAL DB
BASE="/home/james/SequencingData/"
MAP="${BASE}NCBI_RVDB/virus_RVDBs/fixed_URVDB_seqID2.map" # this would require modification to be correct
MAP=/mnt/usersData/RVDB2/RVDB2_seqID.map # this is more complete and has the correct format
CENT_SAVE="~/SequencingData/Centrifuge_libraries/viral/"
cd "${CENT_SAVE}"
centrifuge-build --conversion-table ${MAP} \
                 --taxonomy-tree ${NODE} --name-table ${NAME} \
                 ${NEWFA} ${OUT} -p 8

###############################################################################################
############################## UNCLUSTERED RVDB ################################
###############################################################################################

wget https://rvdb.dbi.udel.edu/download/U-RVDBv29.0.fasta.gz
gunzip U-RVDBv29.0.fasta.gz

BASE1="/home/james/SequencingData/"
BASE2="/mnt/usersData/"

NEWFA="${BASE2}NCBI_RVDB/virus_RVDBs/U-RVDBv29.0.fasta"

NODE="${BASE2}Centrifuge_libraries/viral/taxonomy/nodes.dmp"
NAME="${BASE2}Centrifuge_libraries/viral/taxonomy/names.dmp"
# awk 'NR > 1 {print $2 "\t" $3}' /mnt/usersData/RVDB2/nucl_gb.accession2taxid > /mnt/usersData/RVDB2/RVDB2_seqID.map
MAP=/mnt/usersData/RVDB2/RVDB2_seqID.map
OUT="u-viral29"

sed '/^>/s/^>[^|]*|[^|]*|\([^|]*\)|.*/>\1/' "${NEWFA}" > "${BASE2}NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed.fasta"

awk '/^>/ {if (seqlen){print seq}; print; seq=""; seqlen=0; next} {seq = seq $0; seqlen += length($0)} END {print seq}' "${BASE2}NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed.fasta" > "${BASE2}NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed2.fasta"

# preprocess the FASTA file to ensure each sequence is on a single line.
head -n 3 "${NEWFA}"
head -n 3 "${BASE2}NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed2.fasta"
#######

CENT_SAVE="/${BASE2}Centrifuge_libraries/viral/"
CRNEWFA="${BASE2}NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed2.fasta"
cd "${CENT_SAVE}"

centrifuge-build --conversion-table ${MAP} \
                 --taxonomy-tree ${NODE} --name-table ${NAME} \
                 ${CRNEWFA} ${OUT} -p 2 --dcv 4096

# #######

cd /mnt/usersData/Centrifuge_libraries/viral/taxonomy/
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
rm taxdump.tar.gz

SAVESEQID="${BASE2}NCBI_RVDB/virus_RVDBs/URVDB_V29_seqID.map"
awk -F'|' 'BEGIN {print "code,name"}
/^>/ {
    code=$3;
    $1=$2=$3="";
    gsub(/^[ \t]+|[ \t]+$/, "");  # Remove leading/trailing spaces
    gsub(/\|$/, "");              # Remove trailing pipes
    print code "," $0;
}' "${NEWFA}" > "${SAVESEQID}"
head "${SAVESEQID}"
rg -i "NC_001798.2" "${SAVESEQID}" | head

awk -F',' '{
    if (NF > 2) {
        $2 = $2 "." substr($3, 2);  # Replace the first extra comma with a period
        $3 = "";  # Clear the extra field created
    }
    print $1 "," $2;
}' OFS=',' "${BASE2}NCBI_RVDB/virus_RVDBs/URVDB_V29_seqID.map" > "${BASE2}NCBI_RVDB/virus_RVDBs/fixed_URVDB_V29_seqID.map"
sed -n '37908p' "${BASE2}NCBI_RVDB/virus_RVDBs/fixed_URVDB_V29_seqID.map"
rg -i "LC012610.1" "${BASE2}NCBI_RVDB/virus_RVDBs/fixed_URVDB_V29_seqID.map"

time ./ADVTIG-untargeted/split_fa.py

###############################################################################################
############################## CLUSTERED RVDB ################################
###############################################################################################

BASE="/home/james/SequencingData/"

NEWFA="${BASE}NCBI_RVDB/virus_RVDBs/C-RVDBv29.0.fasta"

NODE="${BASE}Centrifuge_libraries/viral/taxonomy/nodes.dmp"
NAME="${BASE}Centrifuge_libraries/viral/taxonomy/names.dmp"
MAP=/mnt/usersData/RVDB2/RVDB2_seqID.map
OUT="c-viral29"

sed '/^>/s/^>[^|]*|[^|]*|\([^|]*\)|.*/>\1/' "${NEWFA}" > "${BASE}NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed3.fasta"

awk '/^>/ {if (seqlen){print seq}; print; seq=""; seqlen=0; next} {seq = seq $0; seqlen += length($0)} END {print seq}' "${BASE}NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed3.fasta" > "${BASE}NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed4.fasta"

# preprocess the FASTA file to ensure each sequence is on a single line.
head -n 3 "${NEWFA}"
head -n 3 "${BASE}NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed4.fasta"
#######

CENT_SAVE="/${BASE}Centrifuge_libraries/viral/"
CRNEWFA="${BASE}NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed4.fasta"
cd "${CENT_SAVE}"
centrifuge-build --conversion-table ${MAP} \
                 --taxonomy-tree ${NODE} --name-table ${NAME} \
                 ${CRNEWFA} ${OUT} -p 8

# #######

SAVESEQID="${BASE}NCBI_RVDB/virus_RVDBs/CRVDB_V29_seqID.map"
awk -F'|' 'BEGIN {print "code,name"}
/^>/ {
    code=$3;
    $1=$2=$3="";
    gsub(/^[ \t]+|[ \t]+$/, "");  # Remove leading/trailing spaces
    gsub(/\|$/, "");              # Remove trailing pipes
    print code "," $0;
}' "${NEWFA}" > "${SAVESEQID}"
head "${SAVESEQID}"
rg -i "NC_001798.2" "${SAVESEQID}" | head

awk -F',' '{
    if (NF > 2) {
        $2 = $2 "." substr($3, 2);  # Replace the first extra comma with a period
        $3 = "";  # Clear the extra field created
    }
    print $1 "," $2;
}' OFS=',' "${BASE}NCBI_RVDB/virus_RVDBs/CRVDB_V29_seqID.map" > "${BASE}NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map"
sed -n '37908p' "${BASE}NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map"
rg -i "LC012610.1" "${BASE}NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map"

time ./ADVTIG-untargeted/split_fa.py
BASE="/home/james/SequencingData/"
find /home/james/SequencingData/NCBI_RVDB/virus_RVDBs/CRVDB_V29_split_fa/MH5904* -type f -iname "*Human_gammaherpesvirus_4*" | head

SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.txt"
DIRECTORY="/mnt/usersData/ADVTIG/"
CENT_SAVE="~/SequencingData/Centrifuge_libraries/viral/"
CE_DICT_ID="c-viral29"
time ~/SMART-CAMP/mp_metagenomic_assessment_v4.py -d "${DIRECTORY}" -t 5 -ci "${CENT_SAVE}" -c -bl "${CE_DICT_ID}" -m -rs 1 -fd "${SAMPLES}" -skip

awk 'BEGIN{FS="\t"} $2 == 11768 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # felv
awk 'BEGIN{FS="\t"} $2 == 61673 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # perv
awk 'BEGIN{FS="\t"} $2 == 133704 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # pcv1
awk 'BEGIN{FS="\t"} $2 == 11856 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # smrv
awk 'BEGIN{FS="\t"} $2 == 28285 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # adv5
awk 'BEGIN{FS="\t"} $2 == 10376 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # ebv/hh4
awk 'BEGIN{FS="\t"} $2 == 12814 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # rsv
awk 'BEGIN{FS="\t"} $2 == 11250 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # Human orthopneumovirus
awk 'BEGIN{FS="\t"} $2 == 10891 {print FILENAME ": " $0}' /mnt/usersData/ADVTIG/analysis/sample_data/202*/centrifuge/202*_c-viral29_centrifuge_report.tsv # Reovirus

DIRECTORY="/mnt/usersData/ADVTIG/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
CE_DICT_ID="c-viral29"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -e # asm5; extract fastq
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -r # asm5; run analysis
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a # asm5; stats and concatenate
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -l # add labels to the mapped csv (TAXID)

DIRECTORY="/mnt/usersData/ADVTIG/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
CE_DICT_ID="c-viral29"
CE_DICT_ID="uviral25-2"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e # map-ont; model 2; extract fastq
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -r # map-ont; model 2; run analysis
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2-4.csv"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" # map-ont; model 2; stats and concatenate
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2-5.csv"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" # map-ont; model 2; stats and concatenate
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2-6.csv"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" # map-ont; model 2; stats and concatenate
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -l # add labels to the mapped csv


###############################################################################################
############################## CLUSTERED RVDB ################################
###############################################################################################

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

DIRECTORY="/mnt/usersData/ADVTIG_v2_untargeted/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
CE_DICT_ID="uviral25-2"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 -e # asm5; extract fastq
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 -r # asm5; run analysis
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 # asm5; stats and concatenate
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -a -c 2 -l # add labels to the mapped csv

DIRECTORY="/mnt/usersData/ADVTIG_v2_untargeted/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
# SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2-1.csv"
CE_DICT_ID="uviral25-2"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -e # map-ont; model 2; extract fastq
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -r # map-ont; model 2; run analysis
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" # map-ont; model 2; stats and concatenate
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -l # add labels to the mapped csv

time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -dd # deduplicate fq
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -dd -e # extract stats
time ./ADVTIG-untargeted/deep_cov.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" --database "${CE_DICT_ID}" -dd -a # extract stats

time ./ADVTIG-untargeted/filt_low_complexity.py -o "${DIRECTORY}"

# what about a filtered version using edit distance?

time ./ADVTIG-untargeted/summarise_urvdb.py


| **U-RVDB Version** | **Sequence Count (`>` headers)** |
| ------------------ | -------------------------------- |
| v18.0              | 3,076,419                        |
| v19.0              | 3,084,319                        |
| v20.0              | 3,180,577                        |
| v22.0              | 3,551,887                        |
| v25.0              | 8,079,113                        |
| v29.0              | 10,000,605                       |
| vCurrent (flitered)| 2,957,858                        |


gunzip U-RVDBv22.0.fasta.gz
rg -i ">" U-RVDBv22.0.fasta | wc -l

gunzip U-RVDBv18.0.fasta.gz
rg -i ">" U-RVDBv18.0.fasta | wc -l

gunzip U-RVDBv25.0.fasta.gz
rg -i ">" U-RVDBv25.0.fasta | wc -l

gunzip U-RVDBv29.0.fasta.gz
rg -i ">" U-RVDBv29.0.fasta | wc -l

DIRECTORY="/mnt/usersData/ADVTIG/"
SAMPLES="/home/james/SMART-CAMP/configs/viral_DNA_all9-2.csv"
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" -f -e # map-ont; use full fastq file; extract fastq
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" -f -r # map-ont; use full fastq file; run analysis
time ./ADVTIG-untargeted/run_coverage.py -s "${SAMPLES}" -d "${DIRECTORY}" -o "${DIRECTORY}" -f # map-ont; use full fastq file; stats and concatenate

find /mnt/usersData/ADVTIG/asm5/20221222_DNA_advtig-Rep5B-E50_10CFU_41_36 -type f -name "*map.bam" -size +0c | xargs ls -alh

find /mnt/usersData/ADVTIG/analysis/sample_data/20230119_DNA_advtig-Rep1B-E50_10CFU_41_36/fq_seqID -type f -name "*q" -size +0c
find /mnt/usersData/ADVTIG/analysis/sample_data/20230119_DNA_advtig-Rep1B-E50_10CFU_41_36/fq_seqID -type f -name "*q" -size +1M | xargs ls -alh

find /mnt/usersData/ADVTIG/asm5/202*Rep* -type f -name "*uviral2*map-c*t" -size +0c -print > /mnt/usersData/ADVTIG/asm5/names_combined_output_spikes.txt # get file names
find /mnt/usersData/ADVTIG/asm5/202*Rep* -type f -name "*uviral2*map-c*t" -size +0c -print0 | xargs -0 cat > /mnt/usersData/ADVTIG/asm5/combined_output_spikes.txt # concat mapping file summaries

rg -i "primary mapped" /mnt/usersData/ADVTIG/asm5/202* | awk -F '[ :+]' '
{
    count = $2;
    if (count > 0) {
        print $0;
    }
}'

rg -i "Number of reads" /mnt/usersData/ADVTIG/map-ont/202*1B*/**/*uviral2*primary-map*t | less
rg -i "Number of reads" /mnt/usersData/ADVTIG/map-ont/202*2B*/**/*uviral2*primary-map*t | less
rg -i "primary mapped" /mnt/usersData/ADVTIG/map-ont/202*1B*/**/*uviral2*flag*t | less
rg -A 7 -i "AB361382.1" /mnt/usersData/ADVTIG/map-ont/202*1B*/AB361382*/*t

find /mnt/usersData/ADVTIG/map-ont/202*Rep* -type f -name "*uviral2*map-c*t" -size +0c -print > /mnt/usersData/ADVTIG/map-ont/names_combined_output_spikes.txt # get file names
find /mnt/usersData/ADVTIG/map-ont/202*Rep* -type f -name "*uviral2*map-c*t" -size +0c -print0 | xargs -0 cat > /mnt/usersData/ADVTIG/map-ont/combined_output_spikes.txt # concat mapping file summaries
find /mnt/usersData/ADVTIG/map-ont/202*Rep* -type f -name "f*primary-map*t" -size +0c -print0 | xargs -0 cat > /mnt/usersData/ADVTIG/map-ont/prim_full_combined_output_spikes.txt # concat mapping file summaries

##################################################
# Define the output file
OUTPUT_FILE="/mnt/usersData/ADVTIG/map-ont/coverage_primary_read_counts.csv"

# Add the header to the output file
echo "File,Read_Count" > "$OUTPUT_FILE"

# Use `rg` and `awk` to process and format the output
rg -i "primary mapped" /mnt/usersData/ADVTIG/map-ont/202* | awk -F '[:+()]' '
{
    file_path = $1;  # Extract the file path before the ":"
    read_count = $2;  # Extract the primary mapped read count
    gsub(/^[ \t]+|[ \t]+$/, "", read_count);  # Trim leading/trailing whitespace from read_count

    if (read_count > 0) {
        print file_path "," read_count;  # Print file path and read count, separated by a comma
    }
}' >> "$OUTPUT_FILE"

echo "Formatted output saved to $OUTPUT_FILE"

# Verify the output
head "$OUTPUT_FILE"

# Define the input and output files
OUTPUT_FILE="/mnt/usersData/ADVTIG/map-ont/coverage_primary_read_counts.csv"
CLEAN_OUTPUT_FILE="/mnt/usersData/ADVTIG/map-ont/separated_read_counts.csv"

# Add a new header to the output file
echo "Data,Analysis,Sample,SeqID,SeqName,DB,Read_Count" > "$CLEAN_OUTPUT_FILE"

# Process the input file
awk -F ',' 'NR > 1 {
    split($1, path_parts, "/");   # Split the path by "/"
    read_count = $2;             # Read count is the second column

    # Extract parts of the file path
    root = path_parts[2];
    directory = path_parts[4];
    subdirectory = path_parts[5];
    sample = path_parts[6];
    seqID = path_parts[7];
    full_seqName = path_parts[8];

    # Split seqName to extract DB and clean name
    split(full_seqName, name_parts, "-");
    seqName = name_parts[1];
    db = name_parts[2];

    # Print formatted output to the new file
    print directory "," subdirectory "," sample "," seqID "," seqName "," db "," read_count >> "'$CLEAN_OUTPUT_FILE'"
}' "$OUTPUT_FILE"

echo "Processed data saved to $CLEAN_OUTPUT_FILE"

# Verify the new file structure
head "$CLEAN_OUTPUT_FILE"


###############################################################################################

#!/bin/bash

# Directory containing FASTQ files
OUTPUT_FILE="/mnt/usersData/ADVTIG/asm5/read_counts.csv";
echo "File,Read_Count" > "$OUTPUT_FILE";
# Enable recursive globbing with wildcards (**)
shopt -s globstar

# Loop through matching files
for fq in /mnt/usersData/ADVTIG/analysis/sample_data/**/fq_seqID/*.fastq; do
    if [[ -f "$fq" ]]; then
        total_lines=$(wc -l < "$fq");
        read_count=$((total_lines / 4));
        echo $fq,"$(basename "$fq"),$read_count" >> "$OUTPUT_FILE";
    fi;
done;
echo "Read counts saved to $OUTPUT_FILE"

cat "${OUTPUT_FILE}"
awk -F',' '$3 > 100' "${OUTPUT_FILE}"

###############################################################################################

```
