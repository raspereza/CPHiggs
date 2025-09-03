#!/bin/bash
# $1 - channel
# $2 - suffix
# $3 - analysis
n=$#
if [[ $n -ne 3 ]]; then
    echo usage : HaddSamples.bash [CHANNEL] [SUFFIX] [ANALYSIS]
    echo CHANNEL = [mt,et,mm,ee]
    echo SUFFIX = [x,x_xtrig,x_ipcut1,etc...]
    echo ANALYSIS = [baseline,ipSig,datacardsPhiCP]
    exit
fi

chan=${1}
suffix=${2}
analysis=${3}
folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection/${analysis}

for era in Run3_2022preEE Run3_2022postEE Run3_2023preBPix Run3_2023postBPix
do
    rm ${folder}/${chan}_${era}_${suffix}.root
    hadd ${folder}/${chan}_${era}_${suffix}.root ${folder}/*_${chan}_${era}_${suffix}.root
done

rm ${folder}/${chan}_Run3_2022_${suffix}.root
hadd ${folder}/${chan}_Run3_2022_${suffix}.root ${folder}/${chan}_Run3_2022preEE_${suffix}.root ${folder}/${chan}_Run3_2022postEE_${suffix}.root

rm ${folder}/${chan}_Run3_2023_${suffix}.root
hadd ${folder}/${chan}_Run3_2023_${suffix}.root ${folder}/${chan}_Run3_2023preBPix_${suffix}.root ${folder}/${chan}_Run3_2023postBPix_${suffix}.root

rm ${folder}/${chan}_Run3_${suffix}.root
hadd ${folder}/${chan}_Run3_${suffix}.root ${folder}/${chan}_Run3_2022_${suffix}.root ${folder}/${chan}_Run3_2023_${suffix}.root
