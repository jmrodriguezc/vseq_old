#!/usr/bin/env python
import argparse, logging, os
import json, re
import subprocess
from pprint import pprint

__author__ = 'jmrodriguezc'

datapath = None
notmatchedfname = 'not_matched_peptides.txt'

def extract_vseq_main_data(datapath):
    ''' Extract the main data of vseq for all tissue'''    
    out = {}
    if os.path.isdir(datapath):
        for d in os.listdir(datapath):
            dpath = datapath+'/'+d
            if os.path.isdir(dpath):
                tissue = re.match('^RH\_([^\_]*)', d).group(1).lower()
                trep = {
                    'raw':  d,
                    'path': dpath,
                }
                if tissue in out:
                    out[tissue].append(trep)
                else:
                    out[tissue] = [trep]
    return out

def extract_relpep_data(stfile):
    '''Extract information from the relationship file'''
    dat = {}
    with open(stfile, 'r') as outfile:
        for line in outfile:
            line = re.sub(r"\n*$", "", line)
            cols = line.split("\t")
            met = cols[0]
            tis = cols[1]
            seq = cols[2]
            cxr = cols[3]
            raw = cols[4]
            scn = cols[5]
            mod = cols[6]
            dsc = cols[7]
            if not tis in dat:
                dat[tis] = {}
            prog = re.compile('^(RH\_[^\_]*\_TMTHF\_FR[0-9]{1})')
            if prog.match(raw):
                r = prog.match(raw).group(1)
                stkey = r + '_' + scn
                dat[tis][stkey] = {
                    'seq': seq,
                    'cxr': cxr,
                    'raw': raw,
                    'scn': scn,
                    'mod': mod,
                    'dsc': dsc
                }

    return dat

def extract_vseq_data(tissue, drelpep, tpaths, vseq):
    ''' Extract the vseq data for a given tissue '''
    global datapath
    vseqt = {
        'data':[]
    }
    vfilter = {}

    low = lambda s: s[:1].lower() + s[1:] if s else ''
    callback = lambda pat: pat.group(1).lower()    
    for stkey in drelpep:
        dpep = drelpep[stkey]
        seq = dpep['seq']
        raw = dpep['raw']
        scn = dpep['scn']
        cxr = dpep['cxr']
        mod = dpep['mod']
        dsc = dpep['dsc']
        # change pepetide sequence:            
        pep = re.sub(r"\[[^\]]*\]", "", seq)
        pep = low(pep) # first char in lowercase
        pep = re.sub(r'([K])', callback, pep) # all K -> k
        res = '-'
        prog2 = re.compile('^([a-zA-Z]+)_([^\$]*)')
        if prog2.match(mod):
            res = prog2.match(mod).group(1)
            mod = prog2.match(mod).group(2)
        prt = dsc
        p = re.match('^\>(sp|tr)\|([^\|]*)\|',dsc).group(2)
        prt = re.sub(r"^\>[^\s]*\s*", "", prt)
        prt = re.sub(r"\s*PE=\d*\s*SV=\d*\s*", "", prt)
        prt += ' ('+p+')'                    

        vsid = raw+'/'+pep+'_'+scn
        vsfile = datapath+'/'+vsid+'.png'
        if os.path.isfile(vsfile):
            vrep = {
                'raw':          raw,
                'scan':         scn,
                'modification': mod,
                'residue':      res,
                'protein':      prt,
                'pdesc':        dsc,
                'pepmass':      seq,
                'peptide':      pep,
                'cxcorr':       cxr,
                'vsfile':       vsfile
            }
            vseqt['data'].append(vrep) # save all
            if seq in vfilter:
                if cxr > vfilter[seq]['cxcorr']:
                    vfilter[seq] = vrep    
            else:
                vfilter[seq] = vrep
        else:
            t = vsfile + "\t" + raw + "\t" + pep + "\t" + scn + "\t" + seq + "\n"
            nfile  = datapath + '/' + notmatchedfname
            with open(nfile, 'a') as outfile:
                outfile.writelines(t)

    tfname = 'vseq-'+tissue.lower()+'.json'
    tfile = datapath+'/'+'vseq-'+tissue.lower()+'.json'
    vseq.append({
        'name': tissue.title(),
        'data': tfname
    })

    # create vseq report with the ppetides with the best cxcorr
    for vk2,vrep in vfilter.items():
        vseqt['data'].append(vrep)

    with open(tfile, 'w') as outfile:
        json.dump(vseqt, outfile, indent=1)


def main(args):
    ''' Main function'''
    global datapath

    logging.info('extract the main data of vseq for all tissues')
    datapath = args.indir
    logging.debug(datapath)
    dpaths = extract_vseq_main_data(datapath)
    logging.debug(dpaths)
    
    logging.info('filter file from the supplementary material')
    drelpep = extract_relpep_data(args.ffile)
    logging.debug(drelpep)

    logging.info('create vseq data for each tissues')
    outvseq = []
    for tissue,tpaths in dpaths.items():
        extract_vseq_data(tissue, drelpep[tissue], tpaths, outvseq)

    logging.info('create file of vseq data for all tissues')
    outvseqfile  = datapath+'/vseq.json'
    with open(outvseqfile, 'w') as outfile:
        json.dump(outvseq, outfile, indent=1)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Create the Vseq data files into the input directory for the use in the website',
        epilog='''
        Example:
            create_vseq_data.py -i data -x data/heart_isotopWithQuant_stickerOUT.impCols.txt data/liver_isotopWithQuant_stickerOUT.data.impCols.txt
        ''')
    parser.add_argument('-i',  '--indir',  required=True, help='Directory with Vseq images whose directories are organize')
    parser.add_argument('-f',  '--ffile', required=True, help='List of files with the relationshiops between Basal and Changes peptides')
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

    main(args)