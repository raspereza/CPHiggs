#!/bin/bash
# $1 - era
# $2 - channel
# $3 - suffix
era=$1
chan=$2
suffix=$3

n=$#
if [[ $n -ne 3 ]]; then
    echo usage : RunSelectionSeq.bash [ERA] [CHANNEL] [SUFFIX]
    echo ERA = [Run3_2022, Run3_2023]
    echo CHANNEL = [mt,et,mm,ee]
    echo SUFFIX = [x,x_xtrig,x_ipcut1,etc...]
    exit
fi


folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection

rm ${folder}/${chan}_${era}_${suffix}.root
#hadd ${folder}/${chan}_${era}_${suffix}.root ${folder}/dy_${chan}_${era}_${suffix}.root ${folder}/top_${chan}_${era}_${suffix}.root ${folder}/vv_${chan}_${era}_${suffix}.root ${folder}/wjets_${chan}_${era}_${suffix}.root ${folder}/data_${chan}_${era}_${suffix}.root
hadd ${folder}/${chan}_${era}_${suffix}.root ${folder}/*_${chan}_${era}*_${suffix}.root
