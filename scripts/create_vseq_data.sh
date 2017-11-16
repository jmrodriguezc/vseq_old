#!/bin/bash

# init vars
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
DATA_DIR="${SCRIPT_DIR}/../data"
WS_DIR="${SCRIPT_DIR}/.."

echo "Extract the matched peptides for Vseq..."
cd ${DATA_DIR} && grep -v "^NA" peptide_changes.tsv > peptide_changes.matched.tsv
echo "finished"

echo "Create the Vseq data files..."
# Create the Vseq data files into the input directory for the use in the website
cd ${DATA_DIR} && rm *.json peptide_changes.not_matched.txt

cd ${WS_DIR} && ./scripts/create_vseq_data.py \
    -i data \
    -f data/peptide_changes.matched.tsv
echo "finished"

read -n 1 -s -r -p "Press any key to continue"