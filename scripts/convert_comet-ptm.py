import argparse, logging, os
import json, re
import subprocess
from pprint import pprint

__author__ = 'jmrodriguezc'

XCORR_THRESHOLD = 1.11

def main(args):
    ''' Main function'''    

    logging.info('init changed comet-ptm file')
    outcont = ('FirstScan\tnum\tCharge\texp_neutral_mass\tcalc_neutral_mass\texp_mz\tMonoisotopicMass\tScoreValue\t'
 'ScoreID\tsp_score\tq_score\tions_matched\tions_total\tSequence\tpeptide\tRetentionTime\tmodifications\tdelta_mods\t'
 'b_series_delta_mods\ty_series_delta_mods\tdelta_jumps\tdelta_peptide\tprev_aa\tnext_aa\tProteinDescriptions\tduplicate_protein_count\n')

    with open(args.outfile, 'w') as outfile:
        outfile.write(outcont)
    
    low = lambda s: s[:1].lower() + s[1:] if s else ''
    callback = lambda pat: pat.group(1).lower()

    logging.info('read comet-ptm file')
    f = open(args.infile,'r')
    inlines = f.readlines()
    f.close()

    if len(inlines) > 1:
        logging.info('read comet-ptm file')
        for line in inlines[2:]:
            line = re.sub(r"\n*$", "", line)
            cols = line.split("\t")
            exp_neutral_mass = cols[3]
            evalue = cols[6]
            xcorr  = cols[7]
            delta_cn = cols[8]
            plain_peptide = cols[13]
            # discard the scans with threshold in the xcorr
            if float(xcorr) >= XCORR_THRESHOLD:
                # obtain the MonoisotopicMass: exp_neutral_mass + 1.0078 (six decimals)
                mono_isotopic_mass = format(float(exp_neutral_mass) + 1.0078, '.6f')
                # delta_cn -> 9 (ScoreID)
                delta_cn = 9
                # change pepetide sequence:            
                plain_peptide = low(plain_peptide) # first char in lowercase
                plain_peptide = re.sub(r'([K])', callback, plain_peptide) # all K -> k
                # replace values
                # replace e-value by MonoisotopicMass
                cols[6] = mono_isotopic_mass
                cols[8] = delta_cn
                cols[13] = plain_peptide
                # join and prints cols
                s = "\t".join(str(x) for x in cols) + "\n"
                with open(args.outfile, 'a') as outfile:
                    outfile.writelines(s)
    else:
        print("empty comet-ptm file")
    
    # logging.info('read comet-ptm file')
    # with open(args.infile, 'r') as infile:
    #     # CometVersion 2017.01 rev. 1	/data_proteomic/comet/data/Mitochondria_conPlastic/RH_muscle_TMTHF_FR6	03/21/2017, 10:55:36 AM	/data_proteomic/comet/data/Mitochondria_conPlastic/ID_uniprot_MusMusculus_dic2016_curated.fasta
    #     # scan	num	charge	exp_neutral_mass	calc_neutral_mass	exp_mz	e-value	xcorr	delta_cn	sp_score	q_score	ions_matched	ions_total	plain_peptide	peptide	retention_time	modifications	delta_mods	b_series_delta_mods	y_series_delta_mods	delta_jumps	delta_peptide	prev_aa	next_aa	protein	duplicate_protein_count
    #     next(infile)
    #     next(infile)
    #     for line in infile:
    #         line = re.sub(r"\n*$", "", line)
    #         cols = line.split("\t")
    #         exp_neutral_mass = cols[3]
    #         evalue = cols[6]
    #         xcorr  = cols[7]
    #         delta_cn = cols[8]
    #         plain_peptide = cols[13]
    #         # discard the scans with threshold in the xcorr
    #         if float(xcorr) >= XCORR_THRESHOLD:
    #             # obtain the MonoisotopicMass: exp_neutral_mass + 1.0078 (six decimals)
    #             mono_isotopic_mass = format(float(exp_neutral_mass) + 1.0078, '.6f')
    #             # delta_cn -> 9 (ScoreID)
    #             delta_cn = 9
    #             # change pepetide sequence:            
    #             plain_peptide = low(plain_peptide) # first char in lowercase
    #             plain_peptide = re.sub(r'([K])', callback, plain_peptide) # all K -> k
    #             # replace values
    #             # replace e-value by MonoisotopicMass
    #             cols[6] = mono_isotopic_mass
    #             cols[8] = delta_cn
    #             cols[13] = plain_peptide
    #             # join and prints cols
    #             s = "\t".join(str(x) for x in cols) + "\n"
    #             with open(args.outfile, 'a') as outfile:
    #                 outfile.writelines(s)            
    

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Convert comet-ptm files to be used by Vseq',
        epilog='''
        Example:
            convert_comet-ptm.py -i ~/d/data/Ratones_Heteroplasmicos_HF/muscle/comet_ptm/RH_muscle_TMTHF_FR1.txt -o ~/d/data/Ratones_Heteroplasmicos_HF/muscle/comet_ptm/RH_muscle_TMTHF_FR1.tsv
        ''')
    parser.add_argument('-i',  '--infile', required=True, help='comet-ptm file')
    parser.add_argument('-o',  '--outfile', required=True, help='modified comet-ptm file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info('start '+os.path.basename(__file__))
    main(args)
    logging.info('end '+os.path.basename(__file__))
