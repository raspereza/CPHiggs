#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Selection of samples for measurement of the SF related to IPSig cut
import ROOT
import math
import yaml
from yaml.loader import SafeLoader
from array import array

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils
import CPHiggs.IP.analysisIP as analysis
from CPHiggs.IP.ScaleFactor import ScaleFactor
import os

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et','ee','mm'])
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',action='store_true')
    parser.add_argument('-applyMTCut','--mtCut',dest='applyMTCut',action='store_true')
    parser.add_argument('-applyMVisCut','--applyMVisCut',dest='applyMVisCut',action='store_true')
    parser.add_argument('-applyIPSigLep1Cut','--applyIPSigLep1Cut',dest='applyIPSigLep1Cut',action='store_true')
    parser.add_argument('-applyIPSigLep2Cut','--applyIPSigLep2Cut',dest='applyIPSigLep2Cut',action='store_true')
    parser.add_argument('-applyIPSigPromptLepSF','--applyIPSigPromptLepSF',dest='applyIPSigPromptLepSF',action='store_true')
    parser.add_argument('-applyIPSigTauLepSF','--applyIPSigTauLepSF',dest='applyIPSigTauLepSF',action='store_true')
    parser.add_argument('-ipSigCut','--ipSigCut',dest='ipSigCut',type=float,default=1.0)
    parser.add_argument('-runOnData','--runOnData',dest='runOnData',action='store_true')
    parser.add_argument('-runOnDY','--runOnDY',dest='runOnDY',action='store_true')
    parser.add_argument('-runOnSignal','--runOnSignal',dest='runOnSignal',action='store_true')
    args = parser.parse_args()

    basedir = '%s/src/CPHiggs/IP'%(os.getenv('CMSSW_BASE'))

    era = args.era
    chan = args.channel
    useCrossTrigger = args.useCrossTrigger
    applyMTCut = args.applyMTCut
    applyMVisCut = args.applyMVisCut
    applyIPSigLep1Cut = args.applyIPSigLep1Cut
    applyIPSigLep2Cut = args.applyIPSigLep2Cut
    applyIPSigPromptLepSF = args.applyIPSigPromptLepSF
    applyIPSigTauLepSF = args.applyIPSigTauLepSF
    runOnData = args.runOnData
    runOnDY = args.runOnDY
    runOnSignal = args.runOnSignal
    
    mTcut = 999999.
    if applyMTCut: mTcut = 70.

    mvisCut = 999999.
    if applyMVisCut: mvisCut = 80.
    
    ptbins = utils.ptbins[chan]
    etabins = utils.etabins[chan]

    print('')
    print('ptbins ',ptbins)
    print('etabins ',etabins)
    print('')
    
    nbinsPt = len(ptbins) - 1
    nbinsEta = len(etabins) - 1
    
    histPtBins = ROOT.TH1D('histPtBins','',nbinsPt,array('d',list(ptbins)))
    histEtaBins = ROOT.TH1D('histEtaBins','',nbinsEta,array('d',list(etabins)))

    dyNames,topNames,vvNames,wjetsNames = utils.samplesEra(era,
                                                           utils.dy_samples,
                                                           utils.top_samples,
                                                           utils.vv_samples,
                                                           utils.wjets_samples)

    oddNames = utils.higgs_samples
    evenNames = utils.higgs_samples
    
    ipSigPromptLepSF = None
    ipSigTauLepSF = None
    folderSF = os.getenv('CMSSW_BASE')+'/src/CPHiggs/IP/ScaleFactors'

    suffixPromptLep = 'PromptMu'
    suffixTauLep = 'TauMu'
    if chan=='ee' or chan=='et':
        suffixPromptLep = 'PromptE'
        suffixTauLep = 'TauE'
        
    fileNamePromptLepSF = '%s/SF_%s_%s.root'%(folderSF,suffixPromptLep,era)
    if applyIPSigPromptLepSF:
        ipSigPromptLepSF = ScaleFactor(filename=fileNamePromptLepSF)

    labels = {
        'Run3_2022': 'Run3_2022All',
        'Run3_2022EE': 'Run3_2022All',
        'Run3_2023': 'Run3_2023All',
        'Run3_2023BPix': 'Run3_2023All'
    }
        
    fileNameTauLepSF = '%s/SF_%s_%s.root'%(folderSF,suffixTauLep,labels[era])
    if applyIPSigTauLepSF:
        ipSigTauLepSF = ScaleFactor(filename=fileNameTauLepSF)
        
    
    ######################
    # definition of cuts #
    ######################
    # mt channel 
    ptLep1Cut = 21.
    etaLep1Cut = 2.4
    ptLep2Cut = 20.
    etaLep2Cut = 2.3
    ptSingleLepTrigger = 26.
    etaSingleLepTrigger = 2.4
    ptLepCrossTrigger = 21.
    ptTauCrossTrigger = 32.
    lepMomScale = 0.002
    lepTauMomScale = 0.01
    ipSigLepCut = args.ipSigCut
    
    # et channel
    if chan=='et':
        ptLep1Cut = 25.
        etaLep1Cut = 2.1
        ptSingleLepTrigger = 32.
        etaSingleLepTrigger = 2.1
        ptLepCrossTrigger = 25.
        ptTauCrossTrigger = 35.
        lepMomScale = 0.025
        lepTauMomScale = 0.04

    # ee channel
    if chan=='ee':
        ptSingleLepTrigger = 32.
        etaSingleLepTrigger = 2.1
        ptLep1Cut = 25.
        etaLep1Cut = 2.1
        ptLep2Cut = 25.
        etaLep2Cut = 2.1
        lepMomScale = 0.025
        lepTauMomScale = 0.04

    # mm channel
    if chan=='mm':
        ptLep1Cut = 25.
        etaLep1Cut = 2.4
        ptLep2Cut = 25.
        etaLep2Cut = 2.4

    cuts = analysis.AnalysisCuts(mtCut=mTcut,
                                 mvisCut=mvisCut,
                                 etaLep1Cut=etaLep1Cut,
                                 ptLep1Cut=ptLep1Cut,
                                 etaLep2Cut=etaLep2Cut,
                                 ptLep2Cut=ptLep2Cut,
                                 ptSingleLepTrigger=ptSingleLepTrigger,
                                 etaSingleLepTrigger=etaSingleLepTrigger,
                                 ptLepCrossTrigger=ptLepCrossTrigger,
                                 ptTauCrossTrigger=ptTauCrossTrigger,
                                 useCrossTrigger=useCrossTrigger,
                                 applyIPSigLep1Cut=applyIPSigLep1Cut,
                                 applyIPSigLep2Cut=applyIPSigLep2Cut,
                                 ipsigLepCut=ipSigLepCut)

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
        dySamples[name].SetConfig(cuts,histPtBins,histEtaBins,
                                  applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                  applyIPSigTauLepSF=applyIPSigTauLepSF,
                                  ipSigPromptLepSF=ipSigPromptLepSF,
                                  ipSigTauLepSF=ipSigTauLepSF)
        
    print('')
    print('initializing Top samples >>')
    topSamples = {} 
    for name in topNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        topSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        topSamples[name].SetConfig(cuts,histPtBins,histEtaBins,
                                   applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                   applyIPSigTauLepSF=applyIPSigTauLepSF,
                                   ipSigPromptLepSF=ipSigPromptLepSF,
                                   ipSigTauLepSF=ipSigTauLepSF)
        
    print('')
    print('initializing VV samples >>')
    vvSamples = {} 
    for name in vvNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        vvSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        vvSamples[name].SetConfig(cuts,histPtBins,histEtaBins,
                                  applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                  applyIPSigTauLepSF=applyIPSigTauLepSF,
                                  ipSigPromptLepSF=ipSigPromptLepSF,
                                  ipSigTauLepSF=ipSigTauLepSF)
        
    print('')
    print('initializing WJets samples >>')
    wjetsSamples = {}
    for name in wjetsNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        wjetsSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        wjetsSamples[name].SetConfig(cuts,histPtBins,histEtaBins,
                                     applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                     applyIPSigTauLepSF=applyIPSigTauLepSF,
                                     ipSigPromptLepSF=ipSigPromptLepSF,
                                     ipSigTauLepSF=ipSigTauLepSF)

    print('')

    print('')
    print('initializing CP-even Higgs samples >>')
    evenSamples = {}
    for name in evenNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        evenSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        evenSamples[name].SetConfig(cuts,histPtBins,histEtaBins,
                                    applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                    applyIPSigTauLepSF=applyIPSigTauLepSF,
                                    ipSigPromptLepSF=ipSigPromptLepSF,
                                    ipSigTauLepSF=ipSigTauLepSF,
                                    applyWeightCP=1)

    print('')

    print('')
    print('initializing CP-odd Higgs samples >>')
    oddSamples = {}
    for name in oddNames:
        norm = lumi*metadata[0][name]['xs']/metadata[0][name]['eff']
        oddSamples[name] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False)
        oddSamples[name].SetConfig(cuts,histPtBins,histEtaBins,
                                   applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                   applyIPSigTauLepSF=applyIPSigTauLepSF,
                                   ipSigPromptLepSF=ipSigPromptLepSF,
                                   ipSigTauLepSF=ipSigTauLepSF,
                                   applyWeightCP=2)

    print('')

    #######################
    ## running selection ##
    #######################
    hists = {}
    if runOnData:
        hists['data']  = analysis.RunSamplesTuple(dataSamples,'data')
    elif runOnDY:
        hists['dy']    = analysis.RunSamplesTuple(dySamples,'dy')
    elif runOnSignal:
        hists['even'] = analysis.RunSamplesTuple(evenSamples,'even')
        hists['odd']  = analysis.RunSamplesTuple(oddSamples,'odd')
    else:
        hists['data']  = analysis.RunSamplesTuple(dataSamples,'data')
        hists['dy']    = analysis.RunSamplesTuple(dySamples,'dy')
        hists['top']   = analysis.RunSamplesTuple(topSamples,'top')
        hists['vv']    = analysis.RunSamplesTuple(vvSamples,'vv')
        hists['wjets'] = analysis.RunSamplesTuple(wjetsSamples,'wjets')
    
    suffix_mt = ''
    suffix_mvis = ''
    suffix_xtrig = ''
    suffix_ip1 = ''
    suffix_ip2 = ''
    suffix_prompt = ''
    suffix_tau = ''
    if useCrossTrigger:
        suffix_xtrig = '_xtrig'
    if applyMVisCut:
        suffix_mvis = '_mvis'
    if applyMTCut:
        suffix_mt = '_mtcut'
    if applyIPSigLep1Cut:
        suffix_ip1 = '_ipcut1'
    if applyIPSigLep2Cut:
        suffix_ip2 = '_ipcut2'
    if applyIPSigPromptLepSF:
        suffix_prompt = '_promptSF'
    if applyIPSigTauLepSF:
        suffix_tau = '_tauSF'
    
    suffix = suffix_mt+suffix_mvis+suffix_xtrig+suffix_ip1+suffix_ip2+suffix_prompt+suffix_tau 
    outputFileName = '%s/selection/%s_%s%s.root'%(basedir,chan,era,suffix)
    if runOnData:
        outputFileName = '%s/selection/data_%s_%s%s.root'%(basedir,chan,era,suffix)
    if runOnDY:
        outputFileName = '%s/selection/dy_%s_%s%s.root'%(basedir,chan,era,suffix)
    if runOnSignal:
        outputFileName = '%s/selection/signal_%s_%s%s.root'%(basedir,chan,era,suffix)
        
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    histPtBins.Write('ptBins')
    histEtaBins.Write('etaBins')
    for hist in hists:
        for histname in hists[hist]:
            hists[hist][histname].Write(histname) 
    outputFile.Close()
