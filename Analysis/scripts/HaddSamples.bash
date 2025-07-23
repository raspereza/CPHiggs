#!/bin/bash
# $1 - chan
# $2 - era
chan=$1
era=$2
folder=/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection
rm ${folder}/${chan}_${era}.root
hadd ${folder}/${chan}_${era}.root ${folder}/dy_${chan}_${era}.root ${folder}/top_${chan}_${era}.root ${folder}/vv_${chan}_${era}.root ${folder}/wjets_${chan}_${era}.root ${folder}/data_${chan}_${era}.root
