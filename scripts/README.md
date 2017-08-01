# 1. Update the Vseq images from the data directory

# 2. Extract the important columns from the "sticker outputs"
cd data

rm aa &
cut -f 1,16,26,36,50 heart_isotopWithQuant_stickerOUT.txt > aa &&
for A in $(ls -d */ | sed 's/\/\s*/\\|/g' | tr '\n' ' ' | sed 's/ //g' | sed s'/\\|$//g'); do echo "grep -e '$A' aa > heart_isotopWithQuant_stickerOUT.impCols.txt" && grep -e "$A" aa > heart_isotopWithQuant_stickerOUT.impCols.txt; done

rm aa &
cut -f 1,16,26,36,50 liver_isotopWithQuant_stickerOUT.data.txt > aa &&
for A in $(ls -d */ | sed 's/\/\s*/\\|/g' | tr '\n' ' ' | sed 's/ //g' | sed s'/\\|$//g'); do echo "grep -e '$A' aa > liver_isotopWithQuant_stickerOUT.data.impCols.txt" && grep -e "$A" aa > liver_isotopWithQuant_stickerOUT.data.impCols.txt; done

# 3. Create the Vseq data files into the input directory for the use in the website
rm data/*.json & ./scripts/create_vseq_data.py -i data -x data/heart_isotopWithQuant_stickerOUT.impCols.txt data/liver_isotopWithQuant_stickerOUT.data.impCols.txt
