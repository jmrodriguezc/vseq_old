#!/bin/bash

# init vars
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
DATA_DIR="${SCRIPT_DIR}/../data"

# Extract the columns from the "sticker outputs"
# - Scan
# - FileName (Raw)
# - CorXcor
# - FinalSeq_Mass

cd ${DATA_DIR} &&
cut -f 1,16,17,36 heart_isotopWithQuant_stickerOUT.txt > aa &&
for A in $(ls -d */ | sed 's/\/\s*/\\|/g' | tr '\n' ' ' | sed 's/ //g' | sed s'/\\|$//g'); do \
    echo "grep -e '$A' aa > heart_isotopWithQuant_stickerOUT.impCols.txt" && \
    grep -e "$A" aa > heart_isotopWithQuant_stickerOUT.impCols.txt; done && \
rm aa

# cd ${DATA_DIR} &&
# cut -f 1,16,17,36 liver_isotopWithQuant_stickerOUT.data.txt > aa &&
# for A in $(ls -d */ | sed 's/\/\s*/\\|/g' | tr '\n' ' ' | sed 's/ //g' | sed s'/\\|$//g'); do \
#     echo "grep -e '$A' aa > liver_isotopWithQuant_stickerOUT.data.impCols.txt" && \
#     grep -e "$A" aa > liver_isotopWithQuant_stickerOUT.data.impCols.txt; done && \
# rm aa

cd ${DATA_DIR} &&
cut -f 1,16,17,36 liver_isotopWithQuant_stickerOUT.stickerOUT.txt > aa &&
for A in $(ls -d */ | sed 's/\/\s*/\\|/g' | tr '\n' ' ' | sed 's/ //g' | sed s'/\\|$//g'); do \
    echo "grep -e '$A' aa > liver_isotopWithQuant_stickerOUT.stickerOUT.impCols.txt" && \
    grep -e "$A" aa > liver_isotopWithQuant_stickerOUT.stickerOUT.impCols.txt; done && \
rm aa



# Extract the columns:
# - Sequence (FinalSeq_Mass)
# - FastaProteinDescription
# - cXcorr
# - Final-PTM-labels

cd ${DATA_DIR} &&
cut -f 1,2,23,32 Supplementary_Table_Heart_changes-PTMs.txt  | tail -n+4 > Supplementary_Table_Heart_changes-PTMs.impCols.txt

cd ${DATA_DIR} &&
cut -f 1,2,23,29 Supplementary_Table_Liver_changes-PTMs.txt  | tail -n+4 > Supplementary_Table_Liver_changes-PTMs.impCols.txt

# Create the Vseq data files into the input directory for the use in the website
cd ${DATA_DIR} && rm *.json minDeltaScans.txt

cd ${SCRIPT_DIR} &&
${SCRIPT_DIR}/create_vseq_data.py \
    -v \
    -i ${DATA_DIR} \
    -x \
       ${DATA_DIR}/heart_isotopWithQuant_stickerOUT.impCols.txt \
       ${DATA_DIR}/liver_isotopWithQuant_stickerOUT.stickerOUT.impCols.txt \
    -x2 \
       ${DATA_DIR}/Supplementary_Table_Heart_changes-PTMs.impCols.txt \
       ${DATA_DIR}/Supplementary_Table_Liver_changes-PTMs.impCols.txt
