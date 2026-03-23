#!/bin/bash
# $1 - channel
# $2 - suffix
# $3 - analysis
n=$#
if [[ $n -ne 3 ]]; then
    echo usage : HaddSamples.bash [CHANNEL] [SUFFIX] [ANALYSIS]
    echo CHANNEL = [mt,et,mm,ee]
    echo SUFFIX = [x,x_xtrig,x_ipcut1,etc...]
    echo ANALYSIS = [baseline,ipSig,datacardsPhiCP,jetFakes]
    exit
fi

chan=${1}
suffix=${2}
analysis=${3}
folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection/${analysis}

rm ${folder}/${chan}_Run3_${suffix}.root
hadd ${folder}/${chan}_Run3_${suffix}.root ${folder}/*_${chan}_Run3_202*_${suffix}.root

