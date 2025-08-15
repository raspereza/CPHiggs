#!/bin/bash
# $1 - channel
# $2 - suffix
chan=$1
suffix=$2

n=$#
if [[ $n -ne 2 ]]; then
    echo usage : HaddSamplesRun3.bash [CHANNEL] [SUFFIX]
    echo CHANNEL = [mt,et,mm,ee]
    echo SUFFIX = [x,x_xtrig,x_ipcut1,etc...]
    exit
fi


folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection

rm ${folder}/${chan}_Run3_${suffix}.root
hadd ${folder}/${chan}_Run3_${suffix}.root ${folder}/${chan}_Run3_2022_${suffix}.root ${folder}/${chan}_Run3_2023_${suffix}.root
