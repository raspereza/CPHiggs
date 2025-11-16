#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Selection of samples for measurement of the SF related to IPSig cut
import ROOT
import math
import yaml
from yaml.loader import SafeLoader
from array import array
import correctionlib

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils
import CPHiggs.Analysis.analysis as analysis
from IPcorrectionsRun3.IPsignificance.ScaleFactor import ScaleFactor
from IPcorrectionsRun3.FakeFactors.FakeFactor import FakeFactor
import os

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3_2022preEE',choices=['Run3_2022preEE','Run3_2022postEE','Run3_2023preBPix','Run3_2023postBPix'])
    parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et','ee','mm'])
    parser.add_argument('-applyIPSigLepCut','--applyIPSigLepCut',dest='applyIPSigLepCut',type=int,default=1)
    parser.add_argument('-applyFakeFactor','--applyFakeFactor',dest='applyFakeFactor',type=int,default=1)
    parser.add_argument('-sample','--sample',dest='sample',default='data',choices=['data','ztt_0j','ztt_1j','ztt_2j','zll_0j','zll_1j','zll_2j','zll_incl','zll_ext','wjets','top','vv','st','top_2l2v','top_2l2v_ext','top_lv2q','top_lv2q_ext','ggH_sm','ggH_ps','ggH_mm','qqH','HWplus','HWminus','ZH'])
    parser.add_argument('-analysisType','--analysisType',dest='analysisType',default='baseline',choices=['baseline','ipSig','phiCP','datacardsPhiCP','jetFakes'])
    parser.add_argument('-ff_version','--ff_version',dest='ff_version',default='v3')
    args = parser.parse_args()

    eras = utils.periods[args.era]
    chan = args.channel
    sample = args.sample
    analysisType = args.analysisType

    useCrossTrigger = False
    applyBVeto = False
    applyIPSigLep1Cut = False
    applyIPSigLep2Cut = False
    applyIPSigPromptLepSF = False
    applyIPSigTauLepSF = False
    applyIPSigJsonSF = False
    applyFakeFactor = False
    if args.applyIPSigLepCut==1:
        applyIPSigLep1Cut = True
        applyIPSigPromptLepSF = True
        applyIPSigTauLepSF = True
        applyIPSigJsonSF = True
    if args.applyFakeFactor==1:
        applyFakeFactor = True

    ff_version = args.ff_version
    baseFolder = '%s/src/CPHiggs/Analysis'%(os.getenv('CMSSW_BASE'))
    tupleFolderPOWHEG = utils.tupleFolderPOWHEG
    tupleFolderV2 = utils.tupleFolderV2
    tupleFolderMG = utils.tupleFolderMG
    outputFolder = '%s'%(utils.outputFolder)
    folderSF = '%s/src/IPcorrectionsRun3/IPsignificance/data'%(os.getenv('CMSSW_BASE'))
    folderIPSigJsonSF = '%s/src/IPcorrectionsRun3/IPsignificance/JSON'%(os.getenv('CMSSW_BASE'))
    
    mvisUpperCut = 999999999.
    mvisLowerCut = 0.
    if chan=='mt':
        mvisLowerCut = 40.
    
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

    #########################
    # Loading scale factors #
    #########################
    # custom classes
    ipSigPromptLepSF = None
    ipSigTauLepSF = None
    # json corrector
    ipSigJsonSF = None
    
    suffixPromptLep = 'PromptMu'
    suffixTauLep = 'TauMu'
    suffixLep = 'muon'
    if chan=='ee' or chan=='et':
        suffixPromptLep = 'PromptE'
        suffixTauLep = 'TauE'
        suffixLep = 'electron'

    if applyIPSigJsonSF: # use json file ->
        fileNameIPSigJson = '%s/IP_Significance_Correction_Run3_2022-2023_%s.json'%(folderIPSigJsonSF,suffixLep)
        cset = correctionlib.CorrectionSet.from_file(fileNameIPSigJson)
        ipSigJsonSF = cset["ipsig_correction"]
        
    else: # use custom interface (ScaleFactor)
        if applyIPSigPromptLepSF:
            fileNamePromptLepSF = '%s/SF_%s_Run3_2022-2023.root'%(folderSF,suffixPromptLep)
            ipSigPromptLepSF = ScaleFactor(filename=fileNamePromptLepSF)
        if applyIPSigTauLepSF:
            fileNameTauLepSF = '%s/SF_%s_Run3_2022-2023.root'%(folderSF,suffixTauLep)
            ipSigTauLepSF = ScaleFactor(filename=fileNameTauLepSF)

    fakeFactorFolder = '%s/src/IPcorrectionsRun3/FakeFactors/data'%(os.getenv('CMSSW_BASE'))
    filenameFakeFactor = '%s/FF_Run3_%s_%s.root'%(fakeFactorFolder,chan,ff_version)
    fakeFactor = None
    if applyFakeFactor:
        fakeFactor = FakeFactor(filename=filenameFakeFactor)
            
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
    antiMu = 4
    antiE = 2
    
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
        antiMu = 4
        antiE = 6

    # ee channel
    if chan=='ee':
        ptSingleLepTrigger = 32.
        etaSingleLepTrigger = 2.1
        ptLep1Cut = 32.
        etaLep1Cut = 2.1
        ptLep2Cut = 25.
        etaLep2Cut = 2.1
        lepMomScale = 0.025
        lepTauMomScale = 0.04

    # mm channel
    if chan=='mm':
        ptLep1Cut = 26.
        etaLep1Cut = 2.4
        ptLep2Cut = 25.
        etaLep2Cut = 2.4

    cuts = analysis.AnalysisCuts(mvisUpperCut=mvisUpperCut,
                                 mvisLowerCut=mvisLowerCut,
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
                                 applyBVeto=applyBVeto,
                                 antiMu=antiMu,
                                 antiE=antiE)

    # reading lumi, cross section,
    # effective number of MC events and
    # filter efficiency
    metadata = {}
    for era in eras:
        yaml_file = baseFolder+'/params/'+era+'.yaml'
        if not os.path.isfile(yaml_file):
            print('file %s not found'%(yaml_file))
            exit()
        metafile = open(yaml_file,'r')
        metadata[era] = list(yaml.load_all(metafile,Loader=SafeLoader))

    Samples = {}
    if sample=='data':
        print('')
        print('initializing data samples >>>')
        dataNames = {}
        for era in eras:
            datasets = utils.muonSamples[era]
            if args.channel=='et' or args.channel=='ee': datasets =  utils.elecSamples[era]
            dataNames[era] = datasets

        for era in eras:
            for name in dataNames[era]:
                sampleName = name
                Samples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,1.0,True,analysisType=analysisType)
    else:
        print('')
        print('initializing %s samples >>'%(sample))
        for era in eras:
            for name in utils.samplesDict[sample][era]:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                Samples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)

    print('')
    print('++++++++++++++++++++++++++++++++++++++++++++')
    print('')
    #######################
    ## running selection ##
    #######################
    hists = {}
    applyWeightCP = False
    if sample in utils.signal_samples: applyWeightCP = True
    for s in Samples:
        if sample=='data':
            Samples[s].SetConfig(cuts,histPtBins,histEtaBins,
                                 applyFakeFactor=applyFakeFactor,
                                 fakeFactor=fakeFactor)
        else:
            Samples[s].SetConfig(cuts,histPtBins,histEtaBins,
                                 applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                 applyIPSigTauLepSF=applyIPSigTauLepSF,
                                 applyIPSigJsonSF=applyIPSigJsonSF,
                                 applyFakeFactor=applyFakeFactor,
                                 ipSigPromptLepSF=ipSigPromptLepSF,
                                 ipSigTauLepSF=ipSigTauLepSF,
                                 ipSigJsonSF=ipSigJsonSF,
                                 fakeFactor=fakeFactor,
                                 applyWeightCP=applyWeightCP)

    
    hists[sample]  = analysis.RunSamplesTuple(Samples,utils.samplesNameDict[sample])
    
    suffix = '_x'
    if useCrossTrigger:
        suffix += '_xtrig'
    if applyBVeto:
        suffix += '_bveto'
    if applyIPSigLep1Cut:
        suffix += '_ipcut1'
    if applyFakeFactor:
        suffix += '_ff_'+ff_version
    
    outputFileName = '%s/selection/%s/%s_%s_%s%s.root'%(outputFolder,analysisType,sample,chan,args.era,suffix)

    print('Saving histograms to file %s'%(outputFileName))
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    histPtBins.Write('ptBins')
    histEtaBins.Write('etaBins')
    for hist in hists:
        for histname in hists[hist]:
            hists[hist][histname].Write(histname) 
    outputFile.Close()
