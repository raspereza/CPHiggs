#! /usr/bin/env python3

import yaml
from yaml.loader import SafeLoader
import CPHiggs.IP.utils as utils

############
#   MAIN   #
############

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-e', '--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-s', '--sample', dest='sample', default='DYto2L_M_50_madgraphMLM')
    args = parser.parse_args()

    yaml_file = utils.tupleFolder+'/params/'+args.era+'.yaml'
    with open(yaml_file,'r') as f:
        data = list(yaml.load_all(f,Loader=SafeLoader))
        lumi = data[0]['lumi']
        print('lumi = %f'%(lumi))
        sample = args.sample
        xsec = data[0][args.sample]['xs']
        nevt = data[0][args.sample]['eff']
        print('%s : xsec=%f   nevt=%f'%(args.sample,xsec,nevt))
