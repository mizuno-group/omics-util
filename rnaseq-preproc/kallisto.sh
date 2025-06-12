#!/bin/bash

########################
# description
# runner for kallisto
# Note
#   the version of kallisto and the downloaded index should be corresponding

# version
ver=1.0.0

########################
# history
# 240710 start writing (Tadahaya Mizuno)

########################
# preparation
# function for help
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTION] index_path fastq_path1 fastq_path2...
  -h          Display help
  -b VALUE    Number of bootstrap samples, default 100
  -t VALUE    Number of threads, default 8
EOM

  exit 2
}

# get filename
get_filename() {
    tmp=$1
    echo ${tmp%%.*}
}

# get absolute path
realpath() {
  case "$1" in /*) ;; *) printf '%s/' "$PWD";; esac; echo "$1"
}

# option check
n_boot=100
n_threads=8
while getopts b:t:hv opt; do
  case "$opt" in
    h)
      usage
      exit
      ;;
    b)
      n_boot=$OPTARG
      ;;
    t)
      n_threads=$OPTARG
      ;;
    v)
      echo "v$ver"
      exit
      ;;
    \?)
      echo '!! Unexpected argument !!'
      exit 1
      ;;
  esac
done
shift $((OPTIND - 1))

# url argument check
# index file check
if [ "`echo $1 | grep index`" ]; then
  : # do nothing
else
  echo "!! Unexpected argument: give index file !!"
fi
# fastq file check
if [ $# -eq 2 ]; then
  pe=false
elif [ $# -eq 3 ]; then
  pe=true
else
  echo "!! Unexpected argument: give 1 or 2 fastq files !!"
  exit 1
fi

########################
# main function

# path handling
idx=`realpath $1` # full path
full1=`realpath $2` # full path
work_dir=`dirname ${full1}` # full path

# kallisto
f1=`basename "${full1}"`
n1=`get_filename ${f1}`
outdir=${work_dir}/KALLISTO_${n1}

if "${pe}"; then
  full2=`realpath $3`
  f2=`basename "${full2}"`
  n2=`get_filename ${f2}`
  kallisto quant \
    -i ${idx} -o ${outdir} \
    -b ${n_boot} -t ${n_threads} \
    ${full1} ${full2}
else
  kallisto quant \
    -i ${idx} -o ${outdir} \
    -b ${n_boot} -t ${n_threads} \
    ${full1}
fi