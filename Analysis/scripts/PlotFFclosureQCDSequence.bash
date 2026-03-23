#!/bin/bash

n=$#
if [[ $n -ne 3 ]]; then
    echo usage : PlotFFclosureQCDSequence.bash [CHANNEL] [REGION] [FF]
    echo CHANNEL = [mt,et]
    echo REGION = [lowmt_os_antiiso,lowmt_ss_antiiso]
    echo FF = [ss_antiiso,os_antiiso]
    exit
fi

chan=$1
region=$2
ff=$3

./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable bdt_signal --nbins 14 --xmin 0.3 --xmax 1.0 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable bdt_fakes --nbins 14 --xmin 0.3 --xmax 1.0
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable bdt_ditau --nbins 14 --xmin 0.3 --xmax 1.0
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable dR --nbins 50 --xmin 0 --xmax 5 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable m_vis --nbins 40 --xmin 0 --xmax 200 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable mt_1 --nbins 14 --xmin 0 --xmax 70 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable pt_1 --nbins 30 --xmin 0 --xmax 150 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable pt_2 --nbins 30 --xmin 0 --xmax 150 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable eta_1 --nbins 50 --xmin -2.5 --xmax 2.5 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable eta_2 --nbins 50 --xmin -2.5 --xmax 2.5 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable met --nbins 40 --xmin 0 --xmax 200 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable n_jets --nbins 7 --xmin -0.5 --xmax 6.5 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable jpt_1 --nbins 27 --xmin 30 --xmax 300 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable jpt_2 --nbins 27 --xmin 30 --xmax 300 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable mjj --nbins 20 --xmin 0 --xmax 1000 
./scripts/PlotFFclosureQCD.py --ff ${ff} --channel ${chan} --region ${region} --variable jdeta --nbins 16 --xmin 0 --xmax 8 
