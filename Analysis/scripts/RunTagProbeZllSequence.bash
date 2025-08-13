#!/bin/bash

n=$#
if [[ $n -ne 3 ]]; then
    echo usage : RunTagProbeZllSequence.bash [ERA] [CHANNEL] [GENERATOR]
    echo ERA = [Run3_2022, Run3_2023]
    echo CHANNEL = [mm,ee]
    echo GENERATOR = [amcatnlo,MG,powheg]
    exit
fi

era=$1 
channel=$2
generator=$3

for binEta in 1 2 3 
do
    ./scripts/TagProbeZll.py --era ${era} --channel ${channel} --generator ${generator} --binEta ${binEta}
done

for binEta in 1 2 3 
do
    ./scripts/TagProbeZll.py --era ${era} --channel ${channel} --generator ${generator}	--binEta ${binEta} --secondLep
done


 
