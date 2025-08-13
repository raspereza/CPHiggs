#!/bin/bash
era=$1
chan=$2
generator=$3
applyIP1=$4
applyIP2=$5
applySF=$6

n=$#
if [[ $n -ne 6 ]]; then
    echo usage : PlotZllSequence.bash [ERA] [CHANNEL] [IP_CUT1] [IP_CUT2] [IP_SF] [GENERATOR]
    echo ERA = [Run3_2022, Run3_2023]
    echo CHANNEL = [mm,ee]
    echo GENERATOR = [amcatnlo,MG,powheg]
    echo IP_CUT1 = [0,1] 
    echo IP_CUT2 = [0,1] 
    echo IP_SF = [0,1] 
    exit
fi

./scripts/PlotZll.py --era ${era} --channel ${chan} --variable m_vis --nbins 50 --xmin 50 --xmax 150 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable met --nbins 30 --xmin 0 --xmax 150 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable pt_1 --nbins 30 --xmin 0 --xmax 150  --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable pt_2 --nbins 18 --xmin 0 --xmax 90 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5  --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable ipsig_1 --nbins 50 --xmin 0 --xmax 10  --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable ipsig_2 --nbins 50 --xmin 0 --xmax 10 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
#./scripts/PlotZll.py --era ${era} --channel ${chan} --variable n_jets --nbins 11 --xmin -0.5 --xmax 10.5 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable met --nbins 30 --xmin 0 --xmax 150 --applyIP1 ${applyIP1} --applyIP2 ${applyIP2} --applySF ${applySF} --generator ${generator}

