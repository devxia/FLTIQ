#!/bin/bash
# Description: FLTIQ-Full-Length Transcript Identification and Quantification Pipeline
# Author: Binbin Xia

usage() {
    echo ""
    echo "Usage: $0 <OPTIONS>"
    echo "Required Parameters:"
    echo "-i <sqanti_classification>      Path to the SQANTI classification"
    echo "-g <sqanti_gtf>                 Path to the SQANTI output GTF annotation"
    echo "-a <known_annotation>           GTF annotation file from database"
    echo "-d <output_dir>                 Path to the output directory."
    echo ""
    exit 1
}

while getopts ":i:g:a:d" i;
do
    case "${i}" in
    i)
        sqanti_classification=$OPTARG;;
    g)
        sqanti_gtf=$OPTARG;;
    a)
        known_annotation=$OPTARG;;
    d)
        output_dir=$OPTARG;;
    esac
done

# Parse input and set defaults
if [[ "$sqanti_classification" == "" || "$sqanti_gtf" == "" || "$known_annotation" == "" || "$output_dir" == "" ]] ; then
    usage
fi




# This bash script looks for the run_merge.py script in its current working directory, if it does not exist then exits
fltiq_script="$(readlink -f $(dirname $0))/run_merge.py"
if [ ! -f "$fltiq_script" ]; then
    echo "Cleaner python script $fltiq_script does not exist."
    exit 1
fi

# Run script with required parameters
python $fltiq_script \
--sqanti_classification=$sqanti_classification \
--sqanti_gtf=$sqanti_gtf \
--known_annotation=$known_annotation \
--output_dir=$output_dir
