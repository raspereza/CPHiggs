#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Creation of datacards
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.utils as utils
import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

samples = ['data','top','vv','zll','ztt','wjets']
bkg_samples = ['top','vv','zll','ztt','wjets']
signal_samples = ['ggH_sm','qqH','HWplus','HWminus','ZH']
cp_hypotheses = ['sm','ps','mm']
histNames = {'chi2_kinfit': [100,0.,10.], 
             'bdt_signal': [7,0.3,1.0],
             'bdt_ditau': [7,0.3,1.0],
             'bdt_fakes': [7,0.3,1.0],
             'pt_1': [40,0.,200.],
             'pt_2': [40,0.,200.],
             'eta_1': [25,-2.5,2.5],
             'eta_2': [25,-2.5,2.5],
             'm_vis': [48, 0.,240.],
             'm_FastMTT': [25,0,250]
}
options = {
    'default': '',
    'kinfit': '_KinFit',
    'genpv': '_GenPV',
}

decayModes = {
    'pi_a1_3pr' : [8,0,360],
    'rho_a1_3pr' : [8,0,360],
    'a1_1pr_a1_3pr' : [4,0,360],
    'a1_3pr_a1_3pr' : [8,0,360],
}
    
bdt_boundaries = [0.3,0.7,0.8,0.9,1.05]

phicp_nbins = 8
phicp_min = 0.
phicp_max = 360.

signs = ['os','ss']
decayModes = ['pi_a1_3pr','rho_a1_3pr','a1_1pr_a1_3pr','a1_3pr_a1_3pr']

def ExtractChi2Fractions(f):

    hists = {}
    for sample in samples: 
        for sign in signs:
            nameHist = '%s_chi2_kinfit_%s_iso_all_all'%(sample,sign)
            name = '%s_%s'%(sample,sign)
            hists[name]=f.Get(nameHist)

    for sample in signal_samples:
        for 

def ExtractHistos(f,**kwargs):

    option = kwargs.get('option','default')

    SuffixCP = options[option]
    histos = {}
    
    # data and background samples 
    for sample in samples: 
        for histName in histNames:
            for sign in signs:
                nameHist = '%s_%s_%s_iso_all_all'%(sample,histName,sign)
                name = '%s_%s_%s'%(sample,histName,sign)
                nb = histNames[histName][0]
                xmin = histNames[histName][1]
                xmax = histNames[histName][2]
                bins = utils.createBins(nb,xmin,xmax)
                histos[name] = utils.rebinHisto(f.Get(nameHist),bins,'_rebinned')
        for dm in decayModes:
            for sign in signs:
                nameHist = '%s_phicp_%s%s_%s_iso_all_all'%(sample,dm,SuffixCP,sign)
                name = '%s_phicp_%s_%s'%(sample,dm,sign)
                nb = decayModes[dm][0]
                xmin = decayModes[dm][1]
                xmax = decayModes[dm][2]
                bins = utils.createBins(nb,xmin,xmax)
                histos[name] = utils.rebinHisto(f.Get(nameHist),bins,'_rebinned')
                
    # signal samples
    for sample in signal_samples:
        for histName in histNames:
            for sign in signs:
                for hypothesis in cp_hypotheses:
                    nameHist = '%s_%s_%s_iso_all_all_%s'%(sample,histName,sign,hypothesis)
                    name = '%s_%s_%s_%s'%(sample,hypothesis,histName,sign)
                    nb = histNames[histName][0]
                    xmin = histNames[histName][1]
                    xmax = histNames[histName][2]
                    bins = utils.createBins(nb,xmin,xmax)
                    histos[name] = utils.rebinHisto(f.Get(nameHist),bins,'_rebinned')
        for dm in decayModes:
            for sign in signs:
                for hypothesis in cp_hypotheses:
                    nameHist = '%s_phicp_%s%s_%s_iso_all_all_%s'%(sample,dm,SuffixCP,sign,hypothesis)
                    name = '%s_%s_phicp_%s_%s'%(sample,dm,hypothesis,sign)
                    nb = decayModes[dm][0]
                    xmin = decayModes[dm][1]
                    xmax = decayModes[dm][2]
                    bins = utils.createBins(nb,xmin,xmax)
                    histos[name] = utils.rebinHisto(f.Get(nameHist),bins,'_rebinned')

    return histos

def ExtractDistributions(hists,var,**kwargs):

    var = kwargs.get('var','chi2_kinfit')
    isAsimov = kwargs.get('isAsimov',True)
    qcdScale = kwargs.get('qcdScale',1.2)
    lumiScale = kwargs.get('lumiScale',1.0)

    outhists = {}
    outhists['TTT'] = hists['top_%s_os'%(var)].Clone('TTT')
    outhists['ZTT'] = hists['ztt_%s_os'%(var)].Clone('ZTT')
    outhists['ZL'] = hists['ztt_%s_os'%(var)].Clone('ZL')
    
    outhists['VVT'] = hists['vv_%s_os'%(var)].Clone('VVT')
    outhists['VVT'].Add(outhists['VVT'],hists['st_%s_os'%(var)],1.,1.)
    outhists['JetFakes'] = hists['wjets_%s_os'].Clone('JetFakes')

    histQCD = hists['data_%s_ss'%(var)]
    for s in bkg_samples:
        histQCD.Add(histQCD,hists['%s_%s_ss'%(s,var)],1.,-1.)

    histQCD->Scale(qcdScale)
    outhists['JetFakes'].Add(outhists['JetFakes'],histQCD) 

    outhists['TTT'].Scale(lumiScale)
    outhists['VVT'].Scale(lumiScale)
    outhists['ZTT'].Scale(lumiScale)
    outhists['JetFakes'].Scale(lumiScale)
    outhists['ZL'].Scale(lumiScale)
        
    for h in hypotheses:
        outhists['ggH_%s_prod_sm_htt125'%(h)] = hists['ggH_sm_%s_os_%s'%(var,h)].Clone('ggH_%s_prod_sm_htt125'%(h))
        outhists['qqH_%s_htt125'%(h)] = hists['qqH_%s_os_%s'%(var,h)].Clone('qqH_%s_htt125'%(h))
        outhists['']

    for decay_hypothesis in decay_hypotheses:
        histograms['WH_%s_htt125'%(decay_hypothesis)] = histograms['HWplus_%s_htt125'%(decay_hypothesis)]
        histograms['WH_%s_htt125'%(decay_hypothesis)].Add(histograms['WH_%s_htt125'%(decay_hypothesis)],
                                                          histograms['HWminus_%s_htt125'%(decay_hypothesis)],1.,1.)
        
    return outhists
                
                
    
if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022postEE',choices=['Run3','Run3_2022postEE'])
    parser.add_argument('-qcdNorm','--qcdNorm',dest='qcdNorm',type=float,default=1.2)
    parser.add_argument('-kinfit','--kinfit',dest='kinfit',action='store_true')
    parser.add_argument('-chi2cut','--chi2cut',dest='chi2cut',action='store_true')
    parser.add_argument('-rescaleL','--rescaleL',dest='rescaleL',action='store_true')
    parser.add_argument('-option','--option',dest='option',default='default',choices=['default','kinfit','genpv'])
    args = parser.parse_args()

    era = args.era
    chan = 'tt'
    kinfit = args.kinfit
    chi2cut = args.chi2cut
    rescaleL = args.rescaleL

    lumi_2022 = 7980.4
    lumi_2022EE = 26671.7
    lumi_2023 = 18063.0
    lumi_2023BPix = 9693.0
    lumi_total = lumi_2022 + lumi_2022EE + lumi_2023 + lumi_2023BPix
    lumiScale = lumi_total/lumi_2022EE
    if era=='Run3_2022postEE':
        lumiScale = 1.0
    

    suffix = 'x_ipcut1'
    basedir = utils.outputFolder
    inputFileName = '%s/selection/datacardsPhiCP/%s_%s_%s.root'%(basedir,chan,era,suffix)
    if os.path.isfile(inputFileName):
        print('')
        print('Loading ROOT file %s'%(inputFileName))
    else:
        print('')
        print('ROOT file %s is not found'%(inputFileName))
        print('Quitting')
        exit()

    inputFile = ROOT.TFile(inputFileName,'read')

    # defining new binning
    nbins_new = utils.nbins_bdt_signal[chan]
    bdt_bins['signal'] = DefineBinning(inputFile,nbins_new)
    print(bdt_bins['signal'])
    exit()
    
    binsMT = [0.,70.,150.]
    hists = utils.extractHistos(inputFile)
    CalibrateBkgs(hists_mt)

    if dryRun:
        exit()

    list_of_histos = ['ZTT','ZL','JetFakes','TTT','VVT',]
    groups_bdt = {}
    groups_phicp = {}
    for cat in cats:
        varname = 'bdt_%s'%(cat)
        print('%s -> '%(varname))
        groups_bdt[cat] = ExtractDatacards(inputFile,
                                           varname,
                                           bdt_bins[cat],
                                           True,
                                           era=era,
                                           scaleWJ=scaleWJ,
                                           scaleDY=scaleDY,
                                           scaleQCD=scaleQCD)
        x_data = groups_bdt['data_obs'].GetSumOfWeights()
        

    sigbins = bdt_bins['signal']
    nbins = len(sigbins)
    for dm in modes_bins:
        varname = 'phicp_vs_bdt_'+dm
        for ibin in range(0,nbins-1):
            sig_bdt_bin = '%1i'%(ibin)
            name_phicp = '%s_bdtBin%s'%(dm,sig_bdt_bin)
            print('%s -> '%(name_phicp),sigbins[ibin],sigbins[ibin+1])            
            phicp_bins = utils.createBins(modes_bins[dm],phicp_min,phicp_max)
            groups_phicp[name_phicp] = ExtractDatacards(inputFile,
                                                        varname,
                                                        phicp_bins,
                                                        False,
                                                        era=era,
                                                        scaleWJ=scaleWJ,
                                                        scaleDY=scaleDY,
                                                        scaleQCD=scaleQCD,
                                                        bdtbin=sig_bdt_bin,
                                                        xmin=sigbins[ibin],
                                                        xmax=sigbins[ibin+1])
#            for hist in groups_phicp[name_phicp]:
#                print(hist,groups_phicp[name_phicp][hist])
            print('')

    suffix_out = ''
    if bVeto:
        suffix_out += '_bveto'
    if useCrossTrigger:
        suffix_out += '_xtrig'
    outputfileName = utils.outputFolder+'/datacardsPhiCP/'+chan+'_'+era+suffix_out+'.root'
    outputfile = ROOT.TFile(outputfileName,'recreate')
    for cat in cats:
        outputfile.cd('')
        name_dir = '%s_mva_%s'%(chan,cats[cat])
        outputfile.mkdir(name_dir)
        outputfile.cd(name_dir)
        for hist in groups_bdt[cat]:
            groups_bdt[cat][hist].Write(hist)

    for dm in modes_map:
        for ibin in range(0,nbins-1):
            sig_bdt_bin = '%1i'%(ibin)
            name_phicp = '%s_bdtBin%s'%(dm,sig_bdt_bin)
            name_dir = '%s_mva_higgs_cat%s_%s%s'%(chan,sig_bdt_bin,chan_lep[chan],modes_map[dm])
            outputfile.cd('')
            outputfile.mkdir(name_dir)
            outputfile.cd(name_dir)
            for hist in groups_phicp[name_phicp]:
                groups_phicp[name_phicp][hist].Write(hist)

    outputfile.Close()
