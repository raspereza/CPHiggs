#!/bin/bash

n=$#
if [[ $n -ne 8 ]]; then
    echo usage : RunSelectionSeq.bash [ERA] [CHANNEL] [ANALYSIS] [XTRIG] [IPCUT1] [IPCUT2] [PROMPT_SF] [TAU_SF]
    echo ERA = [Run3_2022preEE, Run3_2022postEE, Run3_2023preBPix, Run3_2023postBPix]
    echo CHANNEL = [mt,et,mm,ee]
    echo ANALYSIS = [baseline,ipSig,jetFakes,datacards,phiCP]
    echo XTRIG = [0,1]
    echo IPCUT1 = [0,1]
    echo IPCUT2 = [0,1]
    echo PROMPT_SF = [0,1]
    echo TAU_SF = [0,1]
    exit
fi

era=$1 
chan=$2
analysisType=$3
xtrig=$4
ipcut1=$5
ipcut2=$6
promptSF=$7
tauSF=$8

for sample in data top vv wjets ztt_0j ztt_1j ztt_2j zll_0j zll_1j zll_2j zll_incl dy_1j dy_2j dy_3j dy_4j dy_incl zll_powheg 
do
    ./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${ipcut1} ${ipcut2} ${promptSF} ${tauSF}
done

if [ ${chan} == 'mt']
then
    for sample in ztt_powheg
    do
	./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${ipcut1} ${ipcut2} ${promptSF} ${tauSF}
    done
fi

if [ ${chan} == 'et']
then
    for sample in ztt_powheg
    do
	./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${ipcut1} ${ipcut2} ${promptSF} ${tauSF}
    done
fi


if [ ${era} == 'Run3_2022preEE' ]
then
    for sample in dy_ext zll_ext
    do	
	./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${ipcut1} ${ipcut2} ${promptSF} ${tauSF}
    done
fi

if [ ${era} == 'Run3_2022postEE' ]
then
    for sample in dy_ext zll_ext
    do
	./scripts/RunSelection.bash ${era} ${chan} ${sample} ${analysisType} ${xtrig} ${ipcut1} ${ipcut2} ${promptSF} ${tauSF}
    done
fi
