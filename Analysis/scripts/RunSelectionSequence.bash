#!/bin/bash

n=$#
if [[ $n -ne 8 ]]; then
    echo usage : RunSelectionSeq.bash [CHANNEL] [ANALYSIS] [XTRIG] [MTCUT] [BVETO] [IPCUT1] [PROMPT_SF] [TAU_SF]
    echo CHANNEL = [mt,et,mm,ee]
    echo ANALYSIS = [baseline,ipSig,datacardsPhiCP]
    echo XTRIG = [0,1]
    echo MTCUT = [0,1]
    echo BVETO = [0,1]
    echo IPCUT1 = [0,1]
    echo PROMPT_SF = [0,1]
    echo TAU_SF = [0,1]
    exit
fi

chan=${1}
analysisType=${2}
xtrig=${3}
mtcut=${4}
bveto=${5}
ipcut1=${6}
promptSF=${7}
tauSF=${8}

for era in Run3_2022preEE Run3_2022postEE Run3_2023preBPix Run3_2023postBPix
do

    for sample in data ztt_0j ztt_1j ztt_2j zll_0j zll_1j zll_2j zll_incl wjets top vv ggH_sm ggH_ps ggH_mm qqH HWplus HWminus ZH # dy_incl dy_1j dy_2j dy_3j dy_4j
    do
	./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${mtcut} ${bveto} ${ipcut1} ${promptSF} ${tauSF}
    done
    
    if [ ${era} == 'Run3_2022preEE' ]
    then
	for sample in zll_ext # dy_ext
	do	
	    ./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${mtcut} ${bveto} ${ipcut1} ${promptSF} ${tauSF}
	done
    fi

    if [ ${era} == 'Run3_2022postEE' ]
    then
	for sample in zll_ext # dy_ext
	do
	    ./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${mtcut} ${bveto} ${ipcut1} $${promptSF} ${tauSF}
	done
    fi
    
done


