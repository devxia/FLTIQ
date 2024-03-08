#!/bin/bash
# Description: FLTIQ-Full-Length Transcript Identification and Quantification Pipeline
# Author: Binbin Xia

usage() {
    echo ""
    echo "Usage: $0 <OPTIONS>"
    echo "Required Parameters:"
    echo "-d <data_dir>           Path to directory of supporting data"
    echo "-o <output_dir>         Path to a directory that will store the results."
    echo "-p <cpu_number>         Set number of CPUs"
    echo "-s <sqanti_dir>         Path to the SQANTI directory."
    echo "-n <proj_name>          Project will be subdirectory name store the results."
    echo "-r <reference>          Reference genome (Fasta format)."
    echo "-a <annotation>         GTF annotation file."
    echo "-f <fl_reads>           Raw reads in fasta or fastq format. This argument accepts multiple (comma/space separated) files."
    echo "-i <hisat2_idx>         Path to hisat2 index."
    echo "-1 <r1>                 Path to NGS R1."
    echo "-2 <r2>                 Path to NGS R2."
    echo ""
    exit 1
}

while getopts ":d:o:p:s:n:r:a:f:i:1:2" i;
do
    case "${i}" in
    d)
        data_dir=$OPTARG;;
    o)
        output_dir=$OPTARG;;
    p)
        cpu_number=$OPTARG;;
    s)
        sqanti_dir=$OPTARG;;
    n)
        proj_name=$OPTARG;;
    r)
        reference=$OPTARG;;
    a)
        annotation=$OPTARG;;
    f)
        fl_reads=$OPTARG;;
    i)
        hisat2_idx=$OPTARG;;
    1)
        r1=$OPTARG;;
    2)
        r2=$OPTARG;;
    esac
done

# Parse input and set defaults
if [[ "$data_dir" == "" || "$output_dir" == "" || "$sqanti_dir" == "" || "$annotation" == "" ]] ; then
    usage
fi

if [[ "$cpu_number" == "" ]] ; then
    cpu_number="4"
fi



# This bash script looks for the run_FLTIQ.py script in its current working directory, if it does not exist then exits
fltiq_script="$(readlink -f $(dirname $0))/run_FLTIQ.py"
if [ ! -f "$fltiq_script" ]; then
    echo "FLTIQ python script $fltiq_script does not exist."
    exit 1
fi

# Run AlphaFold with required parameters
python $fltiq_script \
--proj_name=$proj_name \
--sqanti_dir=$sqanti_dir \
--output_dir=$output_dir \
--reference=$reference \
--annotation=$annotation \
--is_native_rna=True \
--n_cpu=$n_cpu \
--hisat2_idx=$hisat2_idx \
--r1=$r1 \
--r2=$r2

