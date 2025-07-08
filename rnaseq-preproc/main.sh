#!/bin/bash

########################
# description
# main runner for fastp + kalisto

# version
ver=1.0.0

########################
# history
# 240711 start writing (Tadahaya Mizuno)

########################
# preparation
# function for help
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTION] index_path fastq_dir...
  -h          Display help
  -o VALUE    Output directory, default current
  -b VALUE    Number of bootstrap samples, default 100
  -t VALUE    Number of threads, default 2
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
  tmp_path=${1}/TMPDIR
  if [ -e ${tmp_path} ]; then
      rm -rf ${tmp_path}
  fi
  mkdir ${tmp_path}
}

# option check
outdir_name=""
res_only=false
n_boot=100
n_threads=2 # 250326 8 → 2に変更
while getopts o:b:t:hv opt; do
  case "$opt" in
    h)
      usage
      exit
      ;;
    o)
      outdir_name=$OPTARG
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

# argument check
# index file check
if [ "`echo $1 | grep index`" ]; then
  : # do nothing
else
  echo "!! Unexpected argument: give index file !!"
fi
# fastq path check
if [ $# -eq 2 ]; then
  :
else
  echo "!! Unexpected argument: give index file path and fastq dir path !!"
  exit 1
fi

########################
# main function
# display time
sta=`date +%s`
echo "> start fastp + kallisto"

# path handling
# parent
#   |- workdir
#   |- outdir
index_path=`realpath $1` # full path
work_dir=`realpath $2` # full path
parent=`dirname ${work_dir}` # full path
if [ -z "$outdir_name" ]; then
  outdir="${parent}/RESULT"
else
  outdir="${parent}/${outdir_name}"
fi
# make RESULT directory if it does not exist
if [ ! -e ${outdir} ]; then
  mkdir ${outdir}
fi
# make DONE directory if it does not exist
if [ ! -e "${work_dir}/DONE" ]; then
  mkdir "${work_dir}/DONE"
fi
path_script=`realpath $0`
dir_script=`dirname ${path_script}`
path_fastp="${dir_script}/fastp.sh"
path_kallisto="${dir_script}/kallisto.sh"

# get fastq file list
q1=()
q2=()

for f1 in "${work_dir}/"*_1.*; do
  # extract characteristic part of the file name
  echo "検出されたファイル: $f1" 
  basename1=$(basename "$f1")
  base="${basename1%_1.*}"

  # check if the pair file exists
  f2=$(find "${work_dir}/" -type f -name "${base}_2.*")
  echo "ペアのファイル: $f2"

  if [[ -n "$f2" ]]; then
    q1+=("$f1")
    q2+=("$f2")
  fi
done
l1=${#q1[@]}
l2=${#q2[@]}

echo "q1 count: ${#q1[@]}"
echo "q2 count: ${#q2[@]}"

# make scripts executable
chmod +x ${path_fastp}
chmod +x ${path_kallisto}

# main loop
if [ ${l1} == 0 ] && [ ${l2} == 0 ]; then
  echo "!! No fastq files were found !!"
  exit 1
elif [ ${l2} == 0 ]; then
  echo ">> single-end"
  for ix in ${!q1[@]}; do
    echo "--- iter "$ix" ---"
    # fastp
    echo ">> fastp"
    bash ${path_fastp} ${q1[ix]}
    # kallisto
    sleep 5
    # get fastq files starting with TRIM_
    echo ">> kallisto"
    tmp1=`find ${work_dir} -maxdepth 1 -name "TRIM_*"`
    bash ${path_kallisto} -b ${n_boot} -t ${n_threads} ${index_path} ${tmp1}
    # summarize the result
    sleep 5
    # obtain the path of the result
    # the path is the first one in the list
    res_path=`find ${work_dir} -maxdepth 1 -name "KALLISTO_*" -print -quit`
    # move report files to the above
    mv "${work_dir}/report_"* ${res_path}
    # rename the res_path by removing the prefix
    sleep 60 # waiting for the completion of the trasfer
    fname=`basename ${q1[ix]}`
    fname2=`get_filename ${fname}`
    mv "${res_path}" "${work_dir}/${fname2}"
    # move the result to the outdir
    sleep 60 # waiting for the completion of the trasfer
    # to avoid outdir overwrite
    if [ -z "$outdir_name" ]; then
      outdir="${parent}/RESULT"
    else
      outdir="${parent}/${outdir_name}"
    fi
    mv "${work_dir}/${fname2}" "${outdir}"
    # remove the intermediate files
    rm -rf "${work_dir}/TRIM_"*
    # move the completed input to the DONE directory
    mv "${q1[ix]}" "${work_dir}/DONE"
    echo ">> iter ""$ix"" done"
    echo "sleep for the completion of RAM release..."
    sleep 120 # waiting for the completion of RAM release
  done
elif [ ${l1} == ${l2} ]; then
  echo ">> pair-end"
  for ix in ${!q1[@]}; do
    echo "--- iter "$ix" ---"
    # fastp
    echo ">> fastp"
    bash ${path_fastp} ${q1[ix]} ${q2[ix]} 
    # kallisto
    sleep 5
    # get fastq files starting with TRIM_
    echo ">> kallisto"
    tmp1=`find ${work_dir} -maxdepth 1 -name "TRIM_*_1.*"`
    tmp2=`find ${work_dir} -maxdepth 1 -name "TRIM_*_2.*"`
    bash ${path_kallisto} -b ${n_boot} -t ${n_threads} ${index_path} ${tmp1} ${tmp2}
    # summarize the result
    sleep 5
    # obtain the path of the result
    # the path is the first one in the list
    res_path=`find ${work_dir} -maxdepth 1 -name "KALLISTO_*" -print -quit`
    echo "[DEBUG] res_path = $res_path"
    # move report files to the above
    mv "${work_dir}/report_"* ${res_path}
    # rename the res_path by removing the prefix
    sleep 60 # waiting for the completion of the trasfer
    fname=`basename ${q1[ix]}`
    fname2=`get_filename ${fname}`
    mv "${res_path}" "${work_dir}/${fname2}"
    # move the result to the outdir
    sleep 60 # waiting for the completion of the trasfer
    # to avoid outdir overwrite
    if [ -z "$outdir_name" ]; then
      outdir="${parent}/RESULT"
    else
      outdir="${parent}/${outdir_name}"
    fi
    mv "${work_dir}/${fname2}" "${outdir}"
    # remove the intermediate files
    rm -rf "${work_dir}/TRIM_"*
    # move the completed input to the DONE directory
    mv "${q1[ix]}" "${work_dir}/DONE"
    mv "${q2[ix]}" "${work_dir}/DONE"
    echo ">> iter ""$ix"" done"
    echo "sleep for the completion of RAM release..."
    sleep 300 # waiting for the completion of RAM release
  done
else
  echo "The number of ends were mismatched !!"
  exit 1
fi

# display time
end=`date +%s`
pt=`expr ${end} - ${sta}`
hr=`expr ${pt} / 3600`
pt=`expr ${pt} % 3600`
mi=`expr ${pt} / 60`
se=`expr ${pt} % 60`
echo "> end"
echo "--- Elapsed Time ---"
echo "${hr} h ${mi} m ${se} s"
echo "--------------------"