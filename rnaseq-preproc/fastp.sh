#!/bin/bash

########################
# description
# runner for fastp
# export the trimmed fastq files to the parent directory of the given fastq files

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
Usage: $(basename "$0") [OPTION] fastq_path1 fastq_path2...
  -h          Display help
  -t VALUE    Add a tag for the output, default trim
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

# make tmp_dir under the given
make_tmp() {
  tmp_path=${1}/tmp_dir
  if [ -e ${tmp_path} ]; then
      rm -rf ${tmp_path}
  fi
  mkdir ${tmp_path}
}

# url argument check
if [ $# -eq 1 ]; then
  pe=false
elif [ $# -eq 2 ]; then
  pe=true
else
  echo "!! Unexpected argument: give 1 or 2 fastq files !!"
  exit 1
fi

# option check
tag=TRIM
while getopts t:hv opt; do
  case "$opt" in
    h)
      usage
      exit
      ;;
    t)
      tag=$OPTARG
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

########################
# main function

# path handling
full1=`realpath $1` # full path
work_dir=`dirname ${full1}` # full path

# fastp
f1=`basename "${full1}"`
n1=`get_filename ${f1}`
out1=${work_dir}/${tag}_${f1}
html1=${work_dir}/report_${n1}.html
json1=${work_dir}/report_${n1}.json

if "${pe}"; then
  full2=`realpath $2`
  f2=`basename "${full2}"`
  n2=`get_filename ${f2}`
  out2=${work_dir}/${tag}_${f2}
  fastp \
    --detect_adapter_for_pe \
    -i ${full1} -I ${full2} \
    -o ${out1} -O ${out2} \
    -h ${html1} -j ${json1} \
    -3 -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -f 1 -F 1
else
  fastp \
    --detect_adapter_for_pe \
    -i ${full1} -o ${out1} \
    -h ${html1} -j ${json1} \
    -3 -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -f 1 -F 1
fi