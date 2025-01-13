#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Selection of samples for measurement of the IPSig cut SFs
import ROOT
import math
import yaml
from yaml.loader import SafeLoader
from array import array

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils
import CPHiggs.IP.analysisIP as analysis
import os

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et','ee','mm'])
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',action='store_true')
    parser.add_argument('-mtCut','--mtCut',dest='mtCut',action='store_true')
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    useCrossTrigger = args.useCrossTrigger
    applyMTCut = args.mtCut
    
    mTcut = 9999999.
    if applyMTCut: mTcut = 40.
    
    ptbins = utils.ptbins[chan]
    etabins = utils.etabins[chan]

    print('ptbins ',ptbins)
    print('etabins ',etabins)
    
    nbinsPt = len(ptbins) - 1
    nbinsEta = len(etabins) - 1
    
    histPtBins = ROOT.TH1D('histPtBins','',nbinsPt,array('d',list(ptbins)))
    histEtaBins = ROOT.TH1D('histEtaBins','',nbinsEta,array('d',list(etabins)))
    
    dyNames,topNames,vvNames,wjetsNames = utils.samplesEra(era,
                                                           utils.dy_samples,
                                                           utils.top_samples,
                                                           utils.vv_samples,
                                                           utils.wjets_samples)
    
    # definition of cuts

    # mt channel 
    ptLep1Cut = 21.
    etaLep1Cut = 2.4
    ptLep2Cut = 20.
    etaLep2Cut = 2.3
    ptSingleLepTrigger = 26.
    etaSingleLepTrigger = 2.4
    ptLepCrossTrigger = 21.
    ptTauCrossTrigger = 32.

    # et channel
    if chan=='et':
        ptLep1Cut = 25.
        etaLep1Cut = 2.5
        ptSingleLepTrigger = 31.
        etaSingleLepTrigger = 2.5
        ptLepCrossTrigger = 25.
        ptTauTauCrossTrigger = 35.

    # ee channel
    if chan=='ee':
        ptSingleLepTrigger = 31.
        etaSingleLepTrigger = 2.5
        ptLep1Cut = 25.
        etaLep1Cut = 2.5
        ptLep2Cut = 25.
        etaLep2Cut = 2.5

    # mm channel
    if chan=='mm':
        ptLep1Cut = 21.
        etaLep1Cut = 2.4
        ptLep2Cut = 21.
        etaLep2Cut = 2.4

    cuts = analysis.AnalysisCuts(mtCut=mTcut,
                                 etaLep1Cut=etaLep1Cut,
                                 ptLep1Cut=ptLep1Cut,
                                 etaLep2Cut=etaLep2Cut,
                                 ptLep2Cut=ptLep2Cut,
                                 ptSingleLepTrigger=ptSingleLepTrigger,
                                 etaSingleLepTrigger=etaSingleLepTrigger,
                                 ptLepCrossTrigger=ptLepCrossTrigger,
                                 ptTauCrossTrigger=ptTauCrossTrigger,
                                 useCrossTrigger=useCrossTrigger)

    yaml_file = utils.tupleFolder+'/params/'+era+'.yaml'
    if not os.path.isfile(yaml_file):
        print('file %s not found'%(yaml_file))
        exit()

    metafile = open(yaml_file,'r')
    metadata = list(yaml.load_all(metafile,Loader=SafeLoader))
    lumi = metadata[0]['lumi']
    
    tupleFolder = utils.tupleFolder
    print('')
    print('initializing data samples >>>')
    dataNames = utils.muonSamples[era]
    if args.channel=='et' or args.channel=='ee': dataNames =  utils.elecSamples[era]
    dataSamples = {}
    for name in dataNames:
        dataSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,1.0,True)
        dataSamples[name].SetConfig(cuts,histPtBins,histEtaBins)
        
    print('')
    print('initializing DY samples >>')
    dySamples = {} 
    for name in dyNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        dySamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        dySamples[name].SetConfig(cuts,histPtBins,histEtaBins)
        
    print('')
    print('initializing Top samples >>')
    topSamples = {} 
    for name in topNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        topSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        topSamples[name].SetConfig(cuts,histPtBins,histEtaBins)
        
    print('')
    print('initializing VV samples >>')
    vvSamples = {} 
    for name in vvNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        vvSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        vvSamples[name].SetConfig(cuts,histPtBins,histEtaBins)
        
    print('')
    print('initializing WJets samples >>')
    wjetsSamples = {}
    for name in wjetsNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        wjetsSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        wjetsSamples[name].SetConfig(cuts,histPtBins,histEtaBins)

    print('')
    # running selection
    hists = {}
    hists['data']  = analysis.RunSamplesTuple(dataSamples,'data')
    hists['dy']    = analysis.RunSamplesTuple(dySamples,'dy')
    hists['top']   = analysis.RunSamplesTuple(topSamples,'top')
    hists['vv']    = analysis.RunSamplesTuple(vvSamples,'vv')
    hists['wjets'] = analysis.RunSamplesTuple(wjetsSamples,'wjets')

    suffix = ''
    if useCrossTrigger:
        if applyMTCut:
            suffix = '_xtrig_mtcut'
        else:
            suffix = '_xtrig'
    else:
        if applyMTCut:
            suffix = '_mtcut'
    
    outputFileName = 'selection/%s_%s%s.root'%(chan,era,suffix)
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    for hist in hists:
        for histname in hists[hist]:
            hists[hist][histname].Write(histname) 
    outputFile.Close()
