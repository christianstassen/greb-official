#!/bin/ksh
#
# set -xe

fileslst=$1

m=0

for files in $fileslst; do
  nam1=$(echo $files | awk 'BEGIN{FS="ctl"}{print $1}')

  files2=${nam1}nc

  echo $files
  echo $files2
  if [ $m -eq 0 ]; then
  m=1
  fi
  cdo -f nc import_binary $files $files2
done
