#!/bin/bash

# init vars
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
DATA_DIR="${SCRIPT_DIR}/../data"
WS_DIR="${SCRIPT_DIR}/.."

echo "Extract the columns for BASAL files..."
# - Scan
# - FileName
# - CorXcor
# - Corr_Seq-mass
# - Assigned Modification
# - Tissue
cd ${DATA_DIR} &&
cut -f 1,16,17,37,50,51 Supplementary_Table_1_basalPeptidomes.txt > Supplementary_Table_1_basalPeptidomes.impCols.txt
echo "finished"

echo "Extract the columns for CHANGES-PTM files..."
# - FastaProteinDescription
# - Sequence (FinalSeq_Mass)
# - compl/control
# - p-value (Comp vs Cont)
# - heteropl/control
# - p-value (Heterop vs Cont)
# - Assigned Modification (Final-PTM-labels)
# - Tissue
cd ${DATA_DIR} &&
cut -f 1,2,17,18,19,20,22,23 Supplementary_Table_2_Changes.txt | sed 's/\"//g' > Supplementary_Table_2_Changes.impCols.txt
echo "finished"

echo "Create the relationship of Basal and Changes peptides..."
cd ${WS_DIR} && ./scripts/create_relpep_data.py \
    -ib data/Supplementary_Table_1_basalPeptidomes.impCols.txt \
    -ic data/Supplementary_Table_2_Changes.impCols.txt \
    -o  data/peptide_changes.tsv
echo "finished"

read -n 1 -s -r -p "Press any key to continue"
