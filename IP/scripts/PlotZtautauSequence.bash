#!/bin/bash
era=$1
chan=$2
ip=$3
sf=$4
xtrig=$5

./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable m_vis --nbins 32 --xmin 40 --xmax 200 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable mt_1 --nbins 40 --xmin 0 --xmax 200 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable pt_1 --nbins 24 --xmin 0 --xmax 120 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable pt_2 --nbins 24 --xmin 0 --xmax 120 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable ipsig_1 --nbins 30 --xmin 0 --xmax 6 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable ipsig_2 --nbins 30 --xmin 0 --xmax 6 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable met --nbins 40 --xmin 0 --xmax 200 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable n_jets --nbins 7 --xmin -0.5 --xmax 6.5 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable n_bjets --nbins 5 --xmin -0.5 --xmax 4.5 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} 
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable jpt_1 --nbins 17 --xmin 30 --xmax 200 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable jpt_2 --nbins 17 --xmin 30 --xmax 200 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable mjj --nbins 20 --xmin 0 --xmax 1000 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable jdeta --nbins 16 --xmin 0 --xmax 8 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig}
