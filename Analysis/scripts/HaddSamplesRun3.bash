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

rm ${folder}/${chan}_${suffix}.root
#hadd ${folder}/${chan}_${era}_${suffix}.root ${folder}/dy_${chan}_${era}_${suffix}.root ${folder}/top_${chan}_${era}_${suffix}.root ${folder}/vv_${chan}_${era}_${suffix}.root ${folder}/wjets_${chan}_${era}_${suffix}.root ${folder}/data_${chan}_${era}_${suffix}.root
hadd ${folder}/${chan}_${suffix}.root ${folder}/${chan}_Run3_2022_${suffix}.root ${folder}/${chan}_Run3_2023_${suffix}.root
