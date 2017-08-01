#!/usr/bin/env python3.6
import argparse, logging, os
import json, re
import subprocess
from pprint import pprint

__author__ = 'jmrodriguezc'

datapath = None

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

def extract_sticker_data(stfiles):
    '''Extract the complement information from the sticker output files'''
    stick = {}
    # p1 = subprocess.Popen(["cut", "-f", "1,16,26,36,50", args.stfils[1]], stdout=subprocess.PIPE)
    # p2 = subprocess.Popen(["grep", "-e", "RH_Liver_TMTHF_FR1\|RH_Liver_TMTHF_FR2"], stdin=p1.stdout,stdout=subprocess.PIPE)
    # p1.stdout.close()
    # for line in p2.stdout:
    #     print(line.decode('utf-8').split("\t"))
    for stfile in stfiles:
        stfname = os.path.basename(stfile)
        tissue = re.match('^([^\_]*)\_', stfname).group(1).lower()
        stick[tissue] = {
            'file': stfile,
            'data': {}
        }
        with open(stfile, 'r') as outfile:
            for line in outfile:
                line = re.sub(r"\n*$", "", line)
                cols = line.split("\t")
                scn = cols[0]
                raw = cols[1]
                dsc = cols[2]
                res = cols[3]
                mod = cols[4]
                dsc = re.sub(r"\"*", "", dsc)
                prog = re.compile('^(RH\_[^\_]*\_TMTHF\_FR[0-9]{1})')
                if prog.match(raw):
                    r = prog.match(raw).group(1)
                    stkey = r + '_' + scn
                    stick[tissue]['data'][stkey] = {
                        'scn': scn,
                        'raw': raw,
                        'dsc': dsc,
                        'res': res,
                        'mod': mod
                    }

    return stick

def extract_vseq_data(tissue, tpaths, dstick, vseq):
    ''' Extract the vseq data for a given tissue '''
    global datapath
    vseqt = {
        'data':[]
    }
    for tpath in tpaths:
        vpath = tpath['path']
        for vf in os.listdir(vpath):
            prog = re.compile('^([^\_]*)\_([^\.]*)\.png$')
            if prog.match(vf):
                vpep = prog.match(vf).group(1)
                vscn = prog.match(vf).group(2)                
                vraw = tpath['raw']
                
                vkey = vraw + '_' + vscn
                if vkey in dstick['data']:
                    vmod = dstick['data'][vkey]['mod']
                    prog2 = re.compile('^([a-zA-Z]+)_([^\$]*)')
                    if prog2.match(vmod):
                        vres = prog2.match(vmod).group(1)
                        vmod = prog2.match(vmod).group(2)
                    else:
                        vmod = dstick['data'][vkey]['mod']
                        vres = '-'
                    # vmod = dstick['data'][vkey]['mod']
                    # vres = '-'
                    vdsc = dstick['data'][vkey]['dsc']
                    vprt = vdsc
                    p = re.match('^\>(sp|tr)\|([^\|]*)\|',vdsc).group(2)
                    vprt = re.sub(r"^\>[^\s]*\s*", "", vprt)
                    vprt = re.sub(r"\s*PE=\d*\s*SV=\d*\s*", "", vprt)
                    vprt += ' ('+p+')'                    
                    vrep = {
                        'raw':          vraw,
                        'modification': vmod,
                        'residue':      vres,
                        'protein':      vprt,
                        'pdesc':        vdsc,
                        'peptide':      vpep,
                        'vseq': {
                            'name':    vscn,
                            'path':    vraw,
                            'imgfile': vf
                        }
                    }
                    vseqt['data'].append(vrep)

                # prog2 = re.compile('^(\w+)\_([^\$]*)$')
                # if vkey in dstick['data'] and prog2.match(vmod) :
                #     #vres = re.match('\w*(\w{1})\[\-*\d*\.\d*\]',dstick['data'][vkey]['res']).group(1)
                #     #vres = dstick['data'][vkey]['res']
                #     vres = prog2.match(vmod).group(1)
                #     vmod = prog2.match(vmod).group(2)
                # else:
                #     vmod = dstick['data'][vkey]['mod']
                #     vres = '-'
                # if vkey in dstick['data']:
                # else:
                #     vdsc = '-'
                #     vprt = vdsc
                

    tfname = 'vseq-'+tissue.lower()+'.json'
    tfile = datapath+'/'+'vseq-'+tissue.lower()+'.json'
    vseq.append({
        'name': tissue.title(),
        'data': tfname
    })    
    with open(tfile, 'w') as outfile:
        json.dump(vseqt, outfile, indent=1)


def main(args):
    ''' Main function'''
    global datapath

    logging.info('extract the main data of vseq for all tissues')
    datapath = os.getcwd()+'/'+args.indir
    logging.debug(datapath)
    dpaths = extract_vseq_main_data(datapath)
    logging.debug(dpaths)
    
    logging.info('extract complement info from the sticker outputs')
    dstick = extract_sticker_data(args.stfils)
    logging.debug(dstick)

    logging.info('create vseq data for each tissues')
    outvseq = []
    for tissue,tpaths in dpaths.items():
        extract_vseq_data(tissue,tpaths, dstick[tissue], outvseq)

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
    parser.add_argument('-i', '--indir', required=True, help='Directory with Vseq images whose directories are organize')
    parser.add_argument('-x', '--stfils', required=True, nargs='+', help='List of files with the sticker results')
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