#!/bin/bash
era=$1
chan=$2
applyIP=$3
applySF=$4
isMG=$5

n=$#
if [[ $n -ne 5 ]]; then
    echo usage : PlotZtautauSequence.bash [ERA] [CHANNEL] [IP_CUT] [IP_SF] [IS_MG]
    exit
fi

./scripts/PlotZll.py --era ${era} --channel ${chan} --variable m_vis --nbins 50 --xmin 50 --xmax 150 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable met --nbins 30 --xmin 0 --xmax 150 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable pt_1 --nbins 30 --xmin 0 --xmax 150  --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable pt_2 --nbins 18 --xmin 0 --xmax 90 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5  --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable ipsig_1 --nbins 50 --xmin 0 --xmax 10  --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}
./scripts/PlotZll.py --era ${era} --channel ${chan} --variable ipsig_2 --nbins 50 --xmin 0 --xmax 10 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF} --isMG ${isMG}

