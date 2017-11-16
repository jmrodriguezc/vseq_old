#!/bin/bash

# init vars
indir="$1"
outdir="$2"
ffile="$3"
scriptdir="D:\projects\Vseq\server\scripts"

run_cmd() {
  cmd="$1"
  echo "** $cmd"
  $cmd
  echo ""
}

# extract all 
for cometptmfile in $(find $indir -type f -name "*.txt"); do
  s=${cometptmfile##*/}
  outfname=${s%.txt}
  infile=$cometptmfile
  outfile=$outdir/$outfname.filter.tsv
  intmp=$cometptmfile
  outtmp=$outdir/$outfname.tmp

  # Note: filter the scans from the peptide changes in the paper
  if [[ $ffile != "" ]]; then
    run_cmd "python $scriptdir/filter_comet-ptm.py -i $intmp -f $ffile -o $outtmp"
    infile=$outtmp
  fi

  # convert the file
  run_cmd "python $scriptdir/convert_comet-ptm.py -i $infile -o $outfile"

  # delete tmp files
  run_cmd "rm $outtmp"
done

read -n 1 -s -r -p "Press any key to continue"