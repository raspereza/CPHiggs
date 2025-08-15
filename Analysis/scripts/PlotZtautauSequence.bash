#!/bin/bash
era=$1
chan=$2
generator=$3
xtrig=$4
ip=$5
sf=$6
calibrWJ=$7

n=$#
if [[ $n -ne 7 ]]; then
    echo usage : PlotZtautauSequence.bash [ERA] [CHANNEL] [GENERATOR] [CROSS_TRIGGER] [IP_CUT] [IP_SF] [CALIBRATE_WJ]
    echo ERA = [Run3_2022, Run3_2022preEE, Run3_2022postEE, Run3_2023, Run3_2023preBPix, Run3_2023postBPix]
    echo CHANNEL = [mt,et]
    echo GENERATOR = [amcatnlo,MG,powheg]
    echo CROSS_TRIGGER = [0,1]
    echo IP_CUT = [0,1]
    echo IP_SF = [0,1]
    echo CALIBRATE_WJ [0,1]
    exit
fi

#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable m_vis --nbins 40 --xmin 0 --xmax 200 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable mt_1 --nbins 50 --xmin 0 --xmax 250 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable pt_1 --nbins 30 --xmin 0 --xmax 150 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable pt_2 --nbins 30 --xmin 0 --xmax 150 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable ipsig_1 --nbins 50 --xmin 0 --xmax 10 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable ipsig_2 --nbins 50 --xmin 0 --xmax 10 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable met --nbins 40 --xmin 0 --xmax 200 --applyIP ${ip} --applySF ${sf}  --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable n_jets --nbins 7 --xmin -0.5 --xmax 6.5 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable n_bjets --nbins 5 --xmin -0.5 --xmax 4.5 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ} 
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable jpt_1 --nbins 27 --xmin 30 --xmax 300 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable jpt_2 --nbins 27 --xmin 30 --xmax 300 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable mjj --nbins 20 --xmin 0 --xmax 1000 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
#./scripts/PlotZtautau.py --era ${era} --channel ${chan} --variable jdeta --nbins 16 --xmin 0 --xmax 8 --applyIP ${ip} --applySF ${sf} --useCrossTrigger ${xtrig} --generator ${generator} --calibrWJ ${calibrWJ}
