#!/bin/bash
era=$1
applyIP=$2
applySF=$3

for chan in mm ee
do
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable m_vis --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable pt_1 --nbins 24 --xmin 0 --xmax 120  --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable pt_2 --nbins 24 --xmin 0 --xmax 120 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5  --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable ipsig_1 --nbins 30 --xmin 0 --xmax 6  --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
    ./scripts/PlotZll.py --era ${era} --channel ${chan} --variable ipsig_2 --nbins 30 --xmin 0 --xmax 6 --applyIP1 ${applyIP} --applyIP2 ${applyIP} --applySF ${applySF}
done
