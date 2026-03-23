#!/bin/bash

n=$#
if [[ $n -ne 1 ]]; then
    echo usage : PlotZtautauSequence.bash [ERA] [CHANNEL] [MT_REGION]
    echo CHANNEL = [mt,et]
    exit
fi

chan=$1


./scripts/PlotZtautau.py --channel ${chan} --variable bdt_ditau --nbins 14 --xmin 0.3 --xmax 1.0 
./scripts/PlotZtautau.py --channel ${chan} --variable bdt_fakes --nbins 14 --xmin 0.3 --xmax 1.0 
./scripts/PlotZtautau.py --channel ${chan} --variable bdt_signal --nbins 14 --xmin 0.3 --xmax 1.0 
./scripts/PlotZtautau.py --channel ${chan} --variable m_vis --nbins 40 --xmin 0 --xmax 200 
#./scripts/PlotZtautau.py --channel ${chan} --variable dR --nbins 50 --xmin 0 --xmax 5 
./scripts/PlotZtautau.py --channel ${chan} --variable mt_1 --nbins 30 --xmin 0 --xmax 60 
./scripts/PlotZtautau.py --channel ${chan} --variable pt_1 --nbins 30 --xmin 0 --xmax 150 
./scripts/PlotZtautau.py --channel ${chan} --variable pt_2 --nbins 30 --xmin 0 --xmax 150 
./scripts/PlotZtautau.py --channel ${chan} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5 
./scripts/PlotZtautau.py --channel ${chan} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 
#./scripts/PlotZtautau.py --channel ${chan} --variable ipsig_1 --nbins 50 --xmin 0 --xmax 10 
#./scripts/PlotZtautau.py --channel ${chan} --variable ipsig_2 --nbins 50 --xmin 0 --xmax 10 
./scripts/PlotZtautau.py --channel ${chan} --variable met --nbins 40 --xmin 0 --xmax 200 
./scripts/PlotZtautau.py --channel ${chan} --variable n_jets --nbins 7 --xmin -0.5 --xmax 6.5 
#./scripts/PlotZtautau.py --channel ${chan} --variable n_bjets --nbins 5 --xmin -0.5 --xmax 4.5 
./scripts/PlotZtautau.py --channel ${chan} --variable jpt_1 --nbins 27 --xmin 30 --xmax 300 
./scripts/PlotZtautau.py --channel ${chan} --variable jpt_2 --nbins 27 --xmin 30 --xmax 300 
./scripts/PlotZtautau.py --channel ${chan} --variable mjj --nbins 20 --xmin 0 --xmax 1000 
./scripts/PlotZtautau.py --channel ${chan} --variable jdeta --nbins 16 --xmin 0 --xmax 8 
