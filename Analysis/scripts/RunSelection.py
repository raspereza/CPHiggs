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
    parser.add_argument('-era','--era',dest='era',default='Run3_2022postEE',choices=['Run3_2022preEE','Run3_2022postEE','Run3_2023preBPix','Run3_2023postBPix'])
    parser.add_argument('-channel','--channel',dest='channel',default='tt',choices=['mt','et','ee','mm','tt'])
    parser.add_argument('-applyIPSigLepCut','--applyIPSigLepCut',dest='applyIPSigLepCut',type=int,default=0)
    parser.add_argument('-applyFakeFactor','--applyFakeFactor',dest='applyFakeFactor',type=int,default=0)
    parser.add_argument('-sample','--sample',dest='sample',default='ggH_sm',choices=['data','ztt_0j','ztt_1j','ztt_2j','zll_0j','zll_1j','zll_2j','zll_incl','zll_ext','wjets','top','vv','st','top_2l2v','top_2l2v_ext','top_lv2q','top_lv2q_ext','ggH_sm','ggH_ps','ggH_mm','qqH','HWplus','HWminus','ZH'])
    parser.add_argument('-analysisType','--analysisType',dest='analysisType',default='phiCP',choices=['baseline','ipSig','phiCP','datacardsPhiCP','jetFakes'])
    parser.add_argument('-ff_version','--ff_version',dest='ff_version',default='ipcut',choices=['ipcut','noipcut'])
    parser.add_argument('-cpp_fit','--cpp_fit',dest='cpp_fit',action='store_true')
    parser.add_argument('-phi_scan','--phi_scan',dest='phi_scan',action='store_true')
    args = parser.parse_args()

    eras = utils.periods[args.era]
    chan = args.channel
    sample = args.sample
    analysisType = args.analysisType

    useCrossTrigger = False
    setCPPkinfit = args.cpp_fit
    phiScan = args.phi_scan
    applyBVeto = False
    applyIPSigLep1Cut = False
    applyIPSigLep2Cut = False
    applyFakeFactor = False
    if args.applyIPSigLepCut==1:
        applyIPSigLep1Cut = True
    if args.applyFakeFactor==1:
        applyFakeFactor = True

    ff_version = args.ff_version
    baseFolder = '%s/src/CPHiggs/Analysis'%(os.getenv('CMSSW_BASE'))

    tupleFolder = utils.tupleFolderV2
    if analysisType=='phiCP':
        tupleFolder = utils.tupleFolderPhiCP
    
    outputFolder = '%s'%(utils.outputFolder)
    ipSigCorrFolder = '%s/src/IPcorrectionsRun3/IPsignificance/JSON'%(os.getenv('CMSSW_BASE'))
    fakeFactorFolder = '%s/src/IPcorrectionsRun3/FakeFactors/JSON/Nov18/'%(os.getenv('CMSSW_BASE'))
    
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
    # Loading fake factors  #
    #########################

    fakeFactor = None
    if applyFakeFactor:
        filenameFakeFactor = '%s/FF_Run3_%s_%s.json'%(fakeFactorFolder,chan,ff_version)
        print('Opening file %s'%(filenameFakeFactor))
        cset = correctionlib.CorrectionSet.from_file(filenameFakeFactor)
        fakeFactor = cset['Fake_factors']
            
    ######################
    # definition of cuts #
    ######################
    # mt channel 
    ptLep1Cut = 21.
    etaLep1Cut = 2.4
    ptLep2Cut = 20.
    etaLep2Cut = 2.5
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
                Samples[sampleName] = analysis.analysisSample(tupleFolder,era,chan,name,1.0,True,analysisType=analysisType)
    else:
        print('')
        print('initializing %s samples >>'%(sample))
        for era in eras:
            for name in utils.samplesDict[sample][era]:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                Samples[sampleName] = analysis.analysisSample(tupleFolder,era,chan,name,norm,False,analysisType=analysisType)

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
                                 applyFakeFactor=applyFakeFactor,
                                 fakeFactor=fakeFactor,
                                 applyWeightCP=applyWeightCP,
                                 cppKinFit=setCPPkinfit,
                                 phiScan=phiScan)

    
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
    if analysisType=='phiCP':
        if setCPPkinfit:
            outputFileName = '%s/selection/%s/%s_%s_%s_cpp.root'%(outputFolder,analysisType,sample,chan,args.era)
        else:
            outputFileName = '%s/selection/%s/%s_%s_%s_python.root'%(outputFolder,analysisType,sample,chan,args.era)
            

    print('Saving histograms to file %s'%(outputFileName))
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    histPtBins.Write('ptBins')
    histEtaBins.Write('etaBins')
    for hist in hists:
        for histname in hists[hist]:
            hists[hist][histname].Write(histname) 
    outputFile.Close()
