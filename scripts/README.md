# 1. Update the Vseq images from the data directory.
RH_Heart_TMTHF_FR1/{peptide}_{scan}.png
...
RH_Heart_TMTHF_FR7/{peptide}_{scan}.png
RH_Liver_TMTHF_FR1/{peptide}_{scan}.png
..
RH_Liver_TMTHF_FR7/{peptide}_{scan}.png

# 2. Save into tabular text file the whole Excell documents: 
Supplementary_Table_Heart_changes-PTMs.xlsx -> Supplementary_Table_Heart_changes-PTMs.txt
Supplementary_Table_Liver_changes-PTMs.xlsx -> Supplementary_Table_Liver_changes-PTMs.txt

# 3. Create the Vseq data for the website
./scripts/create_vseq_data.sh

