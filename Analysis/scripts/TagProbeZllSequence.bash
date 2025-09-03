#!/bin/bash
# $1 - channel
n=$#
if [[ $n -ne 2 ]]; then
    echo usage : TagProbeZllSequence.bash [CHANNEL] [GENERATOR]
    echo CHANNEL = [mm,ee]
    echo GENERATOR = [amcatnlo,MG,powheg]
    exit
fi

channel=${1}
generator=${2}

for era in Run3_2022 Run3_2023
do
    for binEta in 1 2 3
    do
	./scripts/TagProbeZll.py --era ${era} --channel ${channel} --binEta ${binEta} --generator ${generator}
    done
done
