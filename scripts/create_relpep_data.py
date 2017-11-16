#!/usr/bin/env python
import argparse, logging, os
import json, re
import subprocess
from pprint import pprint

__author__ = 'jmrodriguezc'


def extract_basal_data(stfile):
    '''Extract the information of BASAL Peptidomes'''
    dat = {}
    # Scan    FileName        CorXcor Corr_Seq-mass   Assigned Modification   MOUSE TISSUE
    with open(stfile, 'r') as outfile:
        next(outfile) # jump first line
        for line in outfile:
            line = re.sub(r"\n*$", "", line)
            cols = line.split("\t")
            scn = cols[0]
            raw = cols[1]
            cxr = cols[2]
            seq = cols[3]
            mod = cols[4]
            tis = cols[5].lower()
            if not tis in dat:
                dat[tis] = {}
            prog = re.compile('^(RH\_[^\_]*\_TMTHF\_FR[0-9]{1})')
            if prog.match(raw):
                r = prog.match(raw).group(1)
                stkey = r + '_' + scn
                if not seq in dat[tis]:
                    dat[tis][seq] = [{
                        'scn': scn,
                        'raw': r,
                        'cxr': cxr,
                        'mod': mod
                    }]
                else:
                    dat[tis][seq].append({
                        'scn': scn,
                        'raw': r,
                        'cxr': cxr,
                        'mod': mod
                    })
    return dat

def extract_changes_data(spfile):
    '''Extract the information of PTM CHANGES'''
    dat = {}
    # FastaProteinDescription, Sequence, compl/control, p-value (Comp vs Cont), heteropl/control, p-value (Heterop vs Cont), Assigned Modification, tissue
    with open(spfile, 'r') as outfile:
        next(outfile) # jump first line
        for line in outfile:
            line = re.sub(r"\n*$", "", line)
            cols = line.split("\t")
            dsc = cols[0]
            seq = cols[1]
            zco = cols[2]
            pco = cols[3]
            zhe = cols[4]
            phe = cols[5]
            mod = cols[6]
            tis = cols[7].lower()
            dsc = re.sub(r"\"*", "", dsc)
            r = {
                'dsc': dsc,
                'zco': zco,
                'pco': pco,
                'zhe': zhe,
                'phe': phe,
                'mod': mod
            }
            if not tis in dat:
                dat[tis] = {}
            if not seq in dat[tis]:
                dat[tis][seq] = [r]
            else:
                dat[tis][seq].append(r)
    return dat

def create_relpep_data(dbasals, dchanges):
    ''' Extract the vseq data for a given tissue '''
    out = ""
    for tis in dchanges:
        dchange = dchanges[tis]
        if not tis in dbasals:
            logging.error(tis+' does not exits in the basal file')
        dbasal = dbasals[tis]
        for res in dchange:
            if res in dbasal:
                meta_info = 'ONE_IN_PTM'
                # when in changes file there are more than one result for the same peptide
                if len(dchange[res]) > 2:
                    meta_info = 'DUP_IN_PTM'
                for dc in dchange[res]:
                    changes_info = ''
                    basal_info = ''
                    d = {}
                    for db in dbasal[res]:
                        if len(d) == 0:
                            d = db
                        if db['cxr'] > d['cxr']:
                            d = db
                    changes_info = dc['zco']+"\t"+dc['pco']+"\t"+dc['zhe']+"\t"+dc['phe']+"\t"+dc['mod']+"\t"+dc['dsc']
                    if len(d) == 0:
                        meta_info += '_NO_BASAL_SC'
                        basal_info = 'NA'+"\t"+'NA'
                    else:
                        if not d['mod'] == dc['mod']:
                            meta_info += '_DIF_MOD_LABEL'
                        basal_info = d['cxr']+"\t"+d['raw']+"\t"+d['scn']
                    out += meta_info+"\t"+tis+"\t"+res+"\t"+basal_info+"\t"+changes_info+"\n"
            else:
                out += 'NA'+"\t"+tis+"\t"+res+"\t"+"\n"
    return out


def main(args):
    ''' Main function'''
   
    logging.info('extract info from basal peptides')
    dbasal = extract_basal_data(args.ibfile)
    logging.debug(dbasal)

    logging.info('extract info from ptm peptides')
    dchanges = extract_changes_data(args.icfile)
    logging.debug(dchanges)

    logging.info('create relationship between the BASAL and CHANGES peptides')
    outcont = create_relpep_data(dbasal, dchanges)

    logging.info('print the file')
    with open(args.outfile, 'w') as outfile:
        outfile.write(outcont)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Create the relationship of Basal and Changes peptides',
        epilog='''
        Example:
            create_relpep_data.py -ib data/Supplementary_Table_1_Heart_Basal_Peptidome.impCols.txt  -ic data/Supplementary_Table_4_Heart_changes.impCols.txt
        ''')
    parser.add_argument('-ib', '--ibfile', required=True, help='File with the basal peptides')
    parser.add_argument('-ic', '--icfile', required=True, help='File with the ptm peptides')
    parser.add_argument('-o',  '--outfile', required=True, help='Output file')
    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')
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

    main(args)