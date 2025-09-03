#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Selection of samples for measurement of the SF related to IPSig cut
import ROOT
import math
import yaml
from yaml.loader import SafeLoader
from array import array

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils
import CPHiggs.Analysis.analysis as analysis
from CPHiggs.Analysis.ScaleFactor import ScaleFactor
import os

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3_2022preEE',choices=['Run3_2022preEE','Run3_2022postEE','Run3_2023preBPix','Run3_2023postBPix'])
    parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et','ee','mm'])
    parser.add_argument('-applyMTCut','--applyMTCut',dest='applyMTCut',type=int,default=0)
    parser.add_argument('-applyMVisCut','--applyMVisCut',dest='applyMVisCut',type=int,default=0)
    parser.add_argument('-useCrossTrigger','--useCrossTrigger',dest='useCrossTrigger',type=int,default=0)
    parser.add_argument('-bVeto','--bVeto',dest='bVeto',type=int,default=0)
    parser.add_argument('-applyIPSigLep1Cut','--applyIPSigLep1Cut',dest='applyIPSigLep1Cut',type=int,default=0)
    parser.add_argument('-applyIPSigLep2Cut','--applyIPSigLep2Cut',dest='applyIPSigLep2Cut',type=int,default=0)
    parser.add_argument('-applyIPSigPromptLepSF','--applyIPSigPromptLepSF',dest='applyIPSigPromptLepSF',type=int,default=0)
    parser.add_argument('-applyIPSigTauLepSF','--applyIPSigTauLepSF',dest='applyIPSigTauLepSF',type=int,default=0)
    parser.add_argument('-sample','--sample',dest='sample',default='data',choices=['data','ztt_0j','ztt_1j','ztt_2j','zll_0j','zll_1j','zll_2j','zll_incl','zll_ext','wjets','top','vv','even','odd','dy_incl','dy_ext','dy_1j','dy_2j','dy_3j','dy_4j','zll_powheg','ztt_powheg','ggH_sm','ggH_ps','ggH_mm','qqH','HWplus','HWminus','ZH'])
    parser.add_argument('-analysisType','--analysisType',dest='analysisType',default='baseline',choices=['baseline','ipSig','phiCP','datacardsPhiCP'])
    args = parser.parse_args()

    eras_2022 = ['Run3_2022','Run3_2022preEE','Run3_2022postEE']
    
    eras = utils.periods[args.era]
    chan = args.channel
    sample = args.sample
    applyMTCut = False
    if args.applyMTCut==1:
        applyMTCut = True
    applyMVisCut = False
    if args.applyMVisCut==1:
        applyMVisCut = False
    analysisType = args.analysisType
    applyBVeto = False
    if args.bVeto==1:
        applyBVeto = True
    
    useCrossTrigger = False
    applyIPSigLep1Cut = False
    applyIPSigLep2Cut = False
    applyIPSigPromptLepSF = False
    applyIPSigTauLepSF = False
    if args.useCrossTrigger==1:
        useCrossTrigger = True
    if args.applyIPSigLep1Cut==1:
        applyIPSigLep1Cut = True
    if args.applyIPSigLep2Cut==1:
        applyIPSigLep2Cut = True
    if args.applyIPSigPromptLepSF==1:
        applyIPSigPromptLepSF = True
    if args.applyIPSigTauLepSF==1:
        applyIPSigTauLepSF = True
    
    baseFolder = '%s/src/CPHiggs/Analysis'%(os.getenv('CMSSW_BASE'))
    tupleFolderPOWHEG = utils.tupleFolderPOWHEG
    tupleFolderV2 = utils.tupleFolderV2
    tupleFolderMG = utils.tupleFolderMG
    outputFolder = '%s'%(utils.outputFolder)
    folderSF = utils.outputFolder+'/ScaleFactors'
    
    mTcut = 999999.
    if applyMTCut: mTcut = 70.

    mvisUpperCut = 999999.
    if applyMVisCut: mvisUpperCut = 80.

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
    ipSigPromptLepSF = None
    ipSigTauLepSF = None

    suffixPromptLep = 'PromptMu'
    suffixTauLep = 'TauMu'
    if chan=='ee' or chan=='et':
        suffixPromptLep = 'PromptE'
        suffixTauLep = 'TauE'

    sf_map = {'Run3_2022' : 'Run3_2022',
              'Run3_2022preEE' : 'Run3_2022',
              'Run3_2022postEE' : 'Run3_2022',
              'Run3_2023' : 'Run3_2023',
              'Run3_2023preBPix' : 'Run3_2023',
              'Run3_2023postBPix' : 'Run3_2023',
              }
        
    fileNamePromptLepSF = '%s/SF_%s_%s.root'%(folderSF,suffixPromptLep,sf_map[args.era])
    if applyIPSigPromptLepSF:
        ipSigPromptLepSF = ScaleFactor(filename=fileNamePromptLepSF,label='promptSF')
    
    fileNameTauLepSF = '%s/SF_%s_%s.root'%(folderSF,suffixTauLep,sf_map[args.era])
    if applyIPSigTauLepSF:
        ipSigTauLepSF = ScaleFactor(filename=fileNameTauLepSF,label='tauSF')
    
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
    antiJet = 7
    ipSigLepCut = 1.00 # don't change
    ipSigTauCut = 1.25 # don't change
    
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
        antiJet = 7

    # ee channel
    if chan=='ee':
        ptSingleLepTrigger = 31.
        etaSingleLepTrigger = 2.1
        ptLep1Cut = 25.
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

    cuts = analysis.AnalysisCuts(mtCut=mTcut,
                                 mvisUpperCut=mvisUpperCut,
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
                                 ipsigLepCut=ipSigLepCut,
                                 ipsigTauCut=ipSigTauCut,
                                 applyBVeto=applyBVeto,
                                 antiMu=antiMu,
                                 antiE=antiE,
                                 antiJet=antiJet)

    # reading lumi, cross section,
    # effective number of MC events and 
    metadata = {}
    for era in eras:
        yaml_file = baseFolder+'/params/'+era+'.yaml'
        if not os.path.isfile(yaml_file):
            print('file %s not found'%(yaml_file))
            exit()
        metafile = open(yaml_file,'r')
        metadata[era] = list(yaml.load_all(metafile,Loader=SafeLoader))

        
    Samples = {}
    print('')
    print('initializing data samples >>>')
    dataNames = {}
    for era in eras:
        datasets = utils.muonSamples[era]
        if args.channel=='et' or args.channel=='ee': datasets =  utils.elecSamples[era]
        dataNames[era] = datasets

    dataSamples = {}
    for era in eras:
        for name in dataNames[era]:
            sampleName = name
            dataSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,1.0,True,analysisType=analysisType)
            dataSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins)
    Samples['data'] = dataSamples
    
    """
    # commenting out MadGraph
    print('')
    print('initializing dy_incl_MG samples >>')
    dyInclSamples = {} 
    for era in eras:
        for name in utils.dy_incl_MG_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            dyInclSamples[sampleName] = analysis.analysisSample(tupleFolderMG,era,chan,name,norm,False,analysisType=analysisType)
            dyInclSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['dy_incl'] = dyInclSamples

    print('')
    print('initializing dy_1j_MG samples >>')
    dy1JSamples = {} 
    for era in eras:
        for name in utils.dy_1j_MG_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            dy1JSamples[sampleName] = analysis.analysisSample(tupleFolderMG,era,chan,name,norm,False,analysisType=analysisType)
            dy1JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['dy_1j'] = dy1JSamples

    print('')
    print('initializing dy_2j_MG samples >>')
    dy2JSamples = {} 
    for era in eras:
        for name in utils.dy_2j_MG_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            dy2JSamples[sampleName] = analysis.analysisSample(tupleFolderMG,era,chan,name,norm,False,analysisType=analysisType)
            dy2JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['dy_2j'] = dy2JSamples

    print('')
    print('initializing dy_3j_MG samples >>')
    dy3JSamples = {} 
    for era in eras:
        for name in utils.dy_3j_MG_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            dy3JSamples[sampleName] = analysis.analysisSample(tupleFolderMG,era,chan,name,norm,False,analysisType=analysisType)
            dy3JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['dy_3j'] = dy3JSamples

    print('')
    print('initializing dy_4j_MG samples >>')
    dy4JSamples = {} 
    for era in eras:
        for name in utils.dy_4j_MG_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            dy4JSamples[sampleName] = analysis.analysisSample(tupleFolderMG,era,chan,name,norm,False,analysisType=analysisType)
            dy4JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['dy_4j'] = dy4JSamples
    """

    print('')
    print('initializing ztt_0j samples >>')
    ztt0JSamples = {} 
    for era in eras:
        for name in utils.ztt_0j_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            ztt0JSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            ztt0JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['ztt_0j'] = ztt0JSamples

    print('')
    print('initializing ztt_1j samples >>')
    ztt1JSamples = {} 
    for era in eras:
        for name in utils.ztt_1j_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            ztt1JSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            ztt1JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['ztt_1j'] = ztt1JSamples

    print('')
    print('initializing ztt_2j samples >>')
    ztt2JSamples = {} 
    for era in eras:
        for name in utils.ztt_2j_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            ztt2JSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            ztt2JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['ztt_2j'] = ztt2JSamples

    print('')
    print('initializing zll_0j samples >>')
    zll0JSamples = {} 
    for era in eras:
        for name in utils.zll_0j_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            zll0JSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            zll0JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['zll_0j'] = zll0JSamples
    
    print('')
    print('initializing zll_1j samples >>')
    zll1JSamples = {} 
    for era in eras:
        for name in utils.zll_1j_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            zll1JSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            zll1JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['zll_1j'] = zll1JSamples
    
    print('')
    print('initializing zll_2j samples >>')
    zll2JSamples = {} 
    for era in eras:
        for name in utils.zll_2j_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            zll2JSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            zll2JSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['zll_2j'] = zll2JSamples

    print('')
    print('initializing zll_incl samples >>')
    zllInclSamples = {} 
    for era in eras:
        for name in utils.zll_incl_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            zllInclSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            zllInclSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                 applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                 applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                 ipSigPromptLepSF=ipSigPromptLepSF,
                                                 ipSigTauLepSF=ipSigTauLepSF)
    Samples['zll_incl'] = zllInclSamples

    zllExtSamples = {}
    dyExtSamples = {}
    if args.era in eras_2022:
        print('')
        print('initializing zll_ext samples >>')
        for era in eras:
            for name in utils.zll_ext_samples[era]:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                zllExtSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                zllExtSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                 applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                 applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                 ipSigPromptLepSF=ipSigPromptLepSF,
                                                 ipSigTauLepSF=ipSigTauLepSF)
        Samples['zll_ext'] = zllExtSamples
        """
        # commenting out MadGraph
        print('')
        print('initializing dy_ext samples >>')
        for era in eras:
            for name in utils.dy_ext_MG_samples[era]:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                dyExtSamples[sampleName] = analysis.analysisSample(tupleFolderMG,era,chan,name,norm,False,analysisType=analysisType)
                dyExtSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                    applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                    applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                    ipSigPromptLepSF=ipSigPromptLepSF,
                                                    ipSigTauLepSF=ipSigTauLepSF)
        Samples['dy_ext'] = dyExtSamples
        """
    
    """
    # commenting out POWHEG
    if chan in ['mm','mt']:
        print('')
        print('initializing zll_powheg samples >>')
        zll_powheg_samples = utils.zmm_powheg_samples[era]
        zllPowhegSamples = {}
        for era in eras:
            for name in zll_powheg_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                zllPowhegSamples[sampleName] = analysis.analysisSample(tupleFolderPOWHEG,era,chan,name,norm,False,analysisType=analysisType)
                zllPowhegSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                       applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                       applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                       ipSigPromptLepSF=ipSigPromptLepSF,
                                                       ipSigTauLepSF=ipSigTauLepSF)        
        Samples['zll_powheg'] = zllPowhegSamples
    
        if chan=='mt':
            print('')
            print('initializing Ztt_powheg samples >>')
            zttPowhegSamples = {}
            for era in eras:
                for name in utils.ztt_powheg_samples[era]:
                    norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                    sampleName = name+'_'+era
                    zttPowhegSamples[sampleName] = analysis.analysisSample(tupleFolderPOWHEG,era,chan,name,norm,False,analysisType=analysisType)
                    zttPowhegSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                           applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                           applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                           ipSigPromptLepSF=ipSigPromptLepSF,
                                                           ipSigTauLepSF=ipSigTauLepSF)        
            Samples['ztt_powheg'] = zttPowhegSamples

    """

    # tt~ background samples
    print('')
    print('initializing Top samples >>')
    topSamples = {}
    for era in eras:
        for name in utils.top_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            topSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            topSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                             applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                             applyIPSigTauLepSF=applyIPSigTauLepSF,
                                             ipSigPromptLepSF=ipSigPromptLepSF,
                                             ipSigTauLepSF=ipSigTauLepSF)        
    Samples['top'] = topSamples

    # VV and single-top background samples
    print('')
    print('initializing VV samples >>')
    vvSamples = {} 
    for era in eras:
        for name in utils.vv_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            vvSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            vvSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                            applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                            applyIPSigTauLepSF=applyIPSigTauLepSF,
                                            ipSigPromptLepSF=ipSigPromptLepSF,
                                            ipSigTauLepSF=ipSigTauLepSF)
    Samples['vv'] = vvSamples
        
    # WJets background samples
    print('')
    print('initializing WJets samples >>')
    wjetsSamples = {}
    for era in eras:
        for name in utils.wjets_samples[era]:
            norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
            sampleName = name+'_'+era
            wjetsSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
            wjetsSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                               applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                               applyIPSigTauLepSF=applyIPSigTauLepSF,
                                               ipSigPromptLepSF=ipSigPromptLepSF,
                                               ipSigTauLepSF=ipSigTauLepSF)
    Samples['wjets'] = wjetsSamples

    # Signal samples (phiCP and datacardsPhiCP) ->
    if analysisType in ['phiCP','datacardsPhiCP']:    
        print('')
        print('initializing CP-even ggH samples >>')
        ggHevenSamples = {}
        for era in eras:
            for name in utils.ggH_even_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                ggHevenSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                ggHevenSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                     applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                     applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                     ipSigPromptLepSF=ipSigPromptLepSF,
                                                     ipSigTauLepSF=ipSigTauLepSF,
                                                     applyWeightCP=True)
        Samples['ggH_sm'] = ggHevenSamples

        print('')
        print('initializing CP-odd Higgs samples >>')
        ggHoddSamples = {}
        for era in eras:
            for name in utils.ggH_odd_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                ggHoddSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                ggHoddSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                    applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                    applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                    ipSigPromptLepSF=ipSigPromptLepSF,
                                                    ipSigTauLepSF=ipSigTauLepSF,
                                                    applyWeightCP=True)
        Samples['ggH_ps'] = ggHoddSamples

        print('')
        print('initializing CP-maxmix Higgs samples >>')
        ggHmaxmixSamples = {}
        for era in eras:
            for name in utils.ggH_maxmix_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                ggHmaxmixSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                ggHmaxmixSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                       applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                       applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                       ipSigPromptLepSF=ipSigPromptLepSF,
                                                       ipSigTauLepSF=ipSigTauLepSF,
                                                       applyWeightCP=True)
        Samples['ggH_mm'] = ggHmaxmixSamples
    
        print('')
        print('initializing qqH samples >>')
        qqHSamples = {}
        for era in eras:
            for name in utils.qqH_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                qqHSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                qqHSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                             applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                             applyIPSigTauLepSF=applyIPSigTauLepSF,
                                             ipSigPromptLepSF=ipSigPromptLepSF,
                                             ipSigTauLepSF=ipSigTauLepSF,
                                             applyWeightCP=True)
        Samples['qqH'] = qqHSamples


        print('')
        print('initializing HWplus samples >>')
        HWplusSamples = {}
        for era in eras:
            for name in utils.HWplus_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                HWplusSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                HWplusSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                    applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                    applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                    ipSigPromptLepSF=ipSigPromptLepSF,
                                                    ipSigTauLepSF=ipSigTauLepSF,
                                                    applyWeightCP=True)
        Samples['HWplus'] = HWplusSamples

        print('')
        print('initializing HWminus samples >>')
        HWminusSamples = {}
        for era in eras:
            for name in utils.HWminus_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                HWminusSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                HWminusSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                     applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                     applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                     ipSigPromptLepSF=ipSigPromptLepSF,
                                                     ipSigTauLepSF=ipSigTauLepSF,
                                                     applyWeightCP=True)
        Samples['HWminus'] = HWminusSamples

        print('')
        print('initializing ZH samples >>')
        ZHSamples = {}
        for era in eras:
            for name in utils.ZH_samples:
                norm = metadata[era][0]['lumi']*metadata[era][0][name]['xs']*metadata[era][0][name]['filter_efficiency']/metadata[era][0][name]['eff']
                sampleName = name+'_'+era
                ZHSamples[sampleName] = analysis.analysisSample(tupleFolderV2,era,chan,name,norm,False,analysisType=analysisType)
                ZHSamples[sampleName].SetConfig(cuts,histPtBins,histEtaBins,
                                                applyIPSigPromptLepSF=applyIPSigPromptLepSF,
                                                applyIPSigTauLepSF=applyIPSigTauLepSF,
                                                ipSigPromptLepSF=ipSigPromptLepSF,
                                                ipSigTauLepSF=ipSigTauLepSF,
                                                applyWeightCP=True)
        Samples['ZH'] = ZHSamples



    print('')
    print('++++++++++++++++++++++++++++++++++++++++++++')
    print('')
    #######################
    ## running selection ##
    #######################
    hists = {}
    hists[sample]  = analysis.RunSamplesTuple(Samples[sample],sample)
    
    suffix_mt = ''
    suffix_bveto = ''
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
    if applyBVeto:
        suffix_bveto = '_bveto'
    if applyIPSigLep1Cut:
        suffix_ip1 = '_ipcut1'
    if applyIPSigLep2Cut:
        suffix_ip2 = '_ipcut2'
    if applyIPSigPromptLepSF:
        suffix_prompt = '_promptSF'
    if applyIPSigTauLepSF:
        suffix_tau = '_tauSF'
    
    suffix = '_x'+suffix_mt+suffix_mvis+suffix_bveto+suffix_xtrig+suffix_ip1+suffix_ip2+suffix_prompt+suffix_tau 
    outputFileName = '%s/selection/%s/%s_%s_%s%s.root'%(outputFolder,analysisType,sample,chan,args.era,suffix)
        
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    histPtBins.Write('ptBins')
    histEtaBins.Write('etaBins')
    for hist in hists:
        for histname in hists[hist]:
            hists[hist][histname].Write(histname) 
    outputFile.Close()
