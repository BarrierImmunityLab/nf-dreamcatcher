#!/bin/bash 

## this is ONLY necessary for STARsolo Unmapped reads fastq output
## the whitelist is looked up by file name in $WHITELISTS_DIR (the whitelists bundled in the container image) first,
## since the absolute path recorded in Log.out is usually not visible inside this container; the recorded path is only a fallback

STARSOLO_DIR=$1
CPUS=${2:-1}
WHITELISTS_DIR=${WHITELISTS_DIR:-/opt/whitelists}
UNCORR=Unmapped_filt_uncorr.R1.fastq

if [[ $STARSOLO_DIR == "" ]]
then
	>&2 echo "Usage: [WHITELISTS_DIR=<dir with whitelists>] ./get_barcode_corrections.sh <starsolo_dir> [cpus]" 
	exit 1
fi 

echo -e "get_barcode_corrections.sh: generating a FASTA file from unique cellular barcodes that need to be corrected against the whitelist.." 
awk 'NR%4==1' $UNCORR | perl -ne '@t=split/_/; print "$t[-2]\n"' | sort | uniq | awk '{print ">"$1; print $1}' > barcodes_to_correct.fa

WL=`grep "soloCBwhitelist" $STARSOLO_DIR/Log.out  | grep RE-DEFINED | awk '{print $2}'`
echo "get_barcode_corrections.sh: found whitelist $WL; making a fasta file .." 
if [[ -s $WHITELISTS_DIR/$(basename $WL) ]]
then
    WL=$WHITELISTS_DIR/$(basename $WL)
elif [[ ! -s $WL ]]
then
    >&2 echo "get_barcode_corrections.sh: ERROR - whitelist $(basename $WL) was found neither in $WHITELISTS_DIR nor at the path recorded in Log.out!"
    exit 1
fi
cat $WL | awk '{print ">"$1; print $1}' > whitelist.fa
bowtie2-build --threads $CPUS whitelist.fa whitelist &> /dev/null 

echo -e "get_barcode_corrections.sh: mapping barcodes to the whitelist .." 
bowtie2 -af -p$CPUS --no-hd --very-sensitive -U barcodes_to_correct.fa -x whitelist > barcodes_vs_whitelist.sam 2> barcode_bowtie2.log
perl -ne '@t=split/\t/; print "$t[0]\t$t[2]\t$1\n" if ($t[5] eq "16M" && ($t[1] == 0 || $t[1] == 256) && m/\tNM:i:([0,1])\t/)' barcodes_vs_whitelist.sam > barcodes_vs_whitelist.tsv
rm barcodes_vs_whitelist.sam whitelist*.bt2
