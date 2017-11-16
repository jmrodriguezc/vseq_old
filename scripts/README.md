# 1. Save into tabular text file the whole Excell documents: 
Supplementary_Table_1_basalPeptidome.xlsx  -> Supplementary_Table_1_basalPeptidome.txt
Supplementary_Table_2_Changes.xlsx         -> Supplementary_Table_2_Changes.txt

# 2. Create the relationship between the Basal and Changes peptides
./create_relpep_data.sh

# 3. Create the data for the Vseq program from Comet-PTM results
./create_cometptm_for_vseq.sh \
    D:\data\Ratones_Heteroplasmicos_HF\muscle\comet_ptm \
    D:\data\Ratones_Heteroplasmicos_HF\muscle\data_for_vseq

## Note: in the case there are peptides that don't match
./create_cometptm_for_vseq.sh \
    D:\data\Ratones_Heteroplasmicos_HF\muscle\comet_ptm \
    D:\data\Ratones_Heteroplasmicos_HF\muscle\data_for_vseq \
    D:\projects\Vseq\server\data\peptide_changes.not_matched.txt

# 4. Execute Vseq program
cd D:\projects\Vseq\code\tags\Vseq_v4.6
# For the moment, we use the following batch file
Vseq.bat

# 5. Copy the Vseq images infor server/data directory.
cp -rp RH_Heart_TMTHF_FR1/{peptide}_{scan}.png server/data/.
...
RH_Heart_TMTHF_FR7/{peptide}_{scan}.png
RH_Liver_TMTHF_FR1/{peptide}_{scan}.png
..
RH_Liver_TMTHF_FR7/{peptide}_{scan}.png

# WARNING: Don't create another folder within 'data' folder different than 'RH_...' !!!

# 6. Create the Vseq data for the website
./create_vseq_data.sh


