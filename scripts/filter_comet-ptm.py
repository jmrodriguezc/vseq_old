import argparse, logging, os
import json, re
import subprocess
from pprint import pprint

__author__ = 'jmrodriguezc'

XCORR_THRESHOLD = 1.11

def main(args):
    ''' Main function''' 

    b = os.path.basename(args.infile)
    raw = os.path.splitext(b)[0]
    logging.info('init raw: '+raw)

    logging.info('read filtered file')
    filt_scans = {}
    filt_scans[raw] = {}
    with open(args.infilter, 'r') as infile:
        for line in infile:
            line = re.sub(r"\n*$", "", line)
            cols = line.split("\t")
            vfile = cols[0]
            r = cols[1]
            scan = cols[3]
            if r == raw:
                filt_scans[raw][scan] = vfile
    logging.debug(filt_scans)
    
    logging.info('read comet-ptm file')
    f = open(args.infile,'r')
    inlines = f.readlines()
    f.close()

    header_line = inlines[:2]
    with open(args.outfile, 'w') as outfile:
        outfile.writelines(header_line)

    for line in inlines[2:]:
        line = re.sub(r"\n*$", "", line)
        cols = line.split("\t")
        scan = cols[0]
        # print only the scans that are within the Supplementary Material
        if scan in filt_scans[raw]:
            with open(args.outfile, 'a') as outfile:
                outfile.writelines(line+"\n")
    

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Convert comet-ptm files to be used by Vseq',
        epilog='''
        Example:
            filter_paper_data_comet-ptm.py -i ~/d/data/Ratones_Heteroplasmicos_HF/muscle/comet_ptm/RH_muscle_TMTHF_FR1.txt -f ~/server/data/peptide_changes.not_matched.txt -o ~/d/data/Ratones_Heteroplasmicos_HF/muscle/comet_ptm/RH_muscle_TMTHF_FR1.filter.tsv
        ''')
    parser.add_argument('-i',  '--infile', required=True, help='comet-ptm result')
    parser.add_argument('-f',  '--infilter', required=True, help='filter file from Supplementary data')
    parser.add_argument('-o',  '--outfile', required=True, help='filtered comet-ptm file')
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
