#!/bin/ksh
#
# set -xe

fileslst=$(ls *$1*.ctl)

m=0

for files in $fileslst; do
  nam1=$(echo $files | awk 'BEGIN{FS="ctl"}{print $1}')
  
  files2=${nam1}nc

  echo $files
  echo $files2
  if [ $m -eq 0 ]; then
  echo correct? y/n
  read n
  m=1
  fi
  if [[ $n = "y" ]]; then   
  cdo -f nc import_binary $files $files2
  fi
done
