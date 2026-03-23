#!/bin/bash

n=$#
if [[ $n -ne 4 ]]; then
    echo usage : RunSelectionSequence.bash [CHANNEL] [ANALYSIS] [IPCUT] [FF]
    echo CHANNEL = [mt,et,mm,ee]
    echo ANALYSIS = [baseline,ipSig,datacardsPhiCP,jetFakes]
    echo IPCUT = [0,1]
    echo FF = [0,1]
    exit
fi

chan=${1}
analysisType=${2}
ipcut=${3}
ff=${4}

for era in Run3_2022preEE Run3_2022postEE Run3_2023preBPix Run3_2023postBPix
do
    
    for sample in data ztt_0j ztt_1j ztt_2j zll_0j zll_1j zll_2j zll_incl wjets top_2l2v top_lv2q vv st # ggH_sm ggH_ps ggH_mm qqH HWplus HWminus ZH 
    do
	./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${ipcut} ${ff}
    done
    
    if [ ${era} == 'Run3_2022preEE' ]
    then
	for sample in zll_ext top_2l2v_ext top_lv2q_ext
	do	
	    ./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${ipcut} ${ff}
	done
    fi
    
    if [ ${era} == 'Run3_2022postEE' ]
    then
	for sample in zll_ext top_2l2v_ext top_lv2q_ext
	do
	    ./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${ipcut} ${ff}
	done
    fi
done


