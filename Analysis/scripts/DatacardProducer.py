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

def DefineBinning(inputFile,nbins_new):

    hist = inputFile.Get('ggH_sm_bdt_signal_os_iso_all_hadtau_sm').Clone('h_sm')
    h_qqH_sm = inputFile.Get('qqH_bdt_signal_os_iso_all_hadtau_sm')
    h_ggH_ps = inputFile.Get('ggH_sm_bdt_signal_os_iso_all_hadtau_ps')
    h_qqH_ps = inputFile.Get('qqH_bdt_signal_os_iso_all_hadtau_ps')
    hist.Add(hist,h_qqH_sm,1.,1.)
    hist.Add(hist,h_ggH_ps,1.,1.)
    hist.Add(hist,h_qqH_ps,1.,1.)
    nbins = hist.GetNbinsX()

    yield_sum = 0.0 
    for ib in range(1,nbins+1):
        yield_sum += hist.GetBinContent(ib)
        
    nx = float(nbins_new)+0.5
     
    average_yield = yield_sum/nx

    bins_new = []
    yields = []
    bins_new.append(0.33)
    ibin = 1
    tot_yield = 0.
    for ib in range(1,nbins_new):
        go = True
        nevents = 0.0
        while go:
            nevents += hist.GetBinContent(ibin)
            nevents_up = nevents + hist.GetBinContent(ibin+1)
            condition = nevents<average_yield and nevents_up>average_yield
            if condition:
                go = False
                diff = abs(average_yield-nevents)
                diff_up = abs(nevents_up-average_yield)
                condition2 = diff<diff_up
                if condition2:
                    bins_new.append(hist.GetXaxis().GetBinLowEdge(ibin+1))
                    yields.append(nevents)
                    tot_yield += nevents
                else:
                    bins_new.append(hist.GetXaxis().GetBinLowEdge(ibin+2))
                    yields.append(nevents_up)
                    tot_yield += nevents_up
                    ibin += 1
                nevents = 0.0
            ibin += 1
    
    bins_new.append(1.0)
    yields.append(yield_sum-tot_yield)
    print('')
    print('Average yield = %3.1f'%(average_yield))
    print('NBins = %1i'%(nbins_new))
    totsum_check = 0.0
    sum_prod_average = average_yield*nx
    for ib in range(0,nbins_new):
        print('%1i -> [%4.2f,%4.2f] : %3.1f'%(ib,bins_new[ib],bins_new[ib+1],yields[ib]))
        totsum_check += yields[ib]
    print('Total yield = %4.2f  %4.2f  %4.2f'%(yield_sum,totsum_check,sum_prod_average))
    print('')
    return bins_new

def CalibrateBkgs(hists):
    
    hist_data = hists['data_mt_1_os_iso_all'].Clone('h_data')

    hist_top  = hists['top_mt_1_os_iso_all'].Clone('h_top')
    hist_vv   = hists['vv_mt_1_os_iso_all'].Clone('h_vv')
    hist_dy   = hists['dy_mt_1_os_iso_all'].Clone('h_dy')
    hist_wj   = hists['wjets_mt_1_os_iso_all'].Clone('h_wj')

    hist_data_lowmt  = hist_data.Clone('hist_data_lowmt')
    hist_data_highmt = hist_data.Clone('hist_data_highmt') 
    
    hist_data_ss = hists['data_mt_1_ss_iso_all']
    hist_qcd     = hist_data_ss.Clone('h_qcd')
    hist_top_ss  = hists['top_mt_1_ss_iso_all']
    hist_vv_ss   = hists['vv_mt_1_ss_iso_all']
    hist_dy_ss   = hists['dy_mt_1_ss_iso_all']
    hist_wj_ss   = hists['wjets_mt_1_ss_iso_all']
    hist_qcd.Add(hist_qcd,hist_top_ss,1.,-1.)
    hist_qcd.Add(hist_qcd,hist_vv_ss,1.,-1.)
    hist_qcd.Add(hist_qcd,hist_dy_ss,1.,-1.)
    hist_qcd.Add(hist_qcd,hist_wj_ss,1.,-1.)

    hist_data_lowmt.Add(hist_data_lowmt,hist_top,1.,-1.)
    hist_data_lowmt.Add(hist_data_lowmt,hist_vv,1.,-1.)
    hist_data_lowmt.Add(hist_data_lowmt,hist_wj,1.,-1.)
    hist_data_lowmt.Add(hist_data_lowmt,hist_qcd,1.,-1.)

    hist_data_highmt.Add(hist_data_highmt,hist_top,1.,-1.)
    hist_data_highmt.Add(hist_data_highmt,hist_vv,1.,-1.)
    hist_data_highmt.Add(hist_data_highmt,hist_dy,1.,-1.)
    hist_data_highmt.Add(hist_data_highmt,hist_qcd,1.,-1.)
    
    hist_data_os_antiiso = hists['data_mt_1_os_antiiso_all'].Clone('h_data_os')
    hist_top_os_antiiso  = hists['top_mt_1_os_antiiso_all'].Clone('h_top_os')
    hist_vv_os_antiiso   = hists['vv_mt_1_os_antiiso_all'].Clone('h_vv_os')
    hist_dy_os_antiiso   = hists['dy_mt_1_os_antiiso_all'].Clone('h_dy_os')
    hist_wj_os_antiiso   = hists['wjets_mt_1_os_antiiso_all'].Clone('h_wj_os')    

    hist_data_ss_antiiso = hists['data_mt_1_ss_antiiso_all'].Clone('h_data_ss')
    hist_top_ss_antiiso  = hists['top_mt_1_ss_antiiso_all'].Clone('h_top_ss')
    hist_vv_ss_antiiso   = hists['vv_mt_1_ss_antiiso_all'].Clone('h_vv_ss')
    hist_dy_ss_antiiso   = hists['dy_mt_1_ss_antiiso_all'].Clone('h_dy_ss')
    hist_wj_ss_antiiso   = hists['wjets_mt_1_ss_antiiso_all'].Clone('h_wj_ss')

    hist_data_os_antiiso.Add(hist_data_os_antiiso,hist_top_os_antiiso,1.,-1.)
    hist_data_os_antiiso.Add(hist_data_os_antiiso,hist_vv_os_antiiso,1.,-1.)
    hist_data_os_antiiso.Add(hist_data_os_antiiso,hist_dy_os_antiiso,1.,-1.)
    hist_data_os_antiiso.Add(hist_data_os_antiiso,hist_wj_os_antiiso,1.,-1.)
    
    hist_data_ss_antiiso.Add(hist_data_ss_antiiso,hist_top_ss_antiiso,1.,-1.)
    hist_data_ss_antiiso.Add(hist_data_ss_antiiso,hist_vv_ss_antiiso,1.,-1.)
    hist_data_ss_antiiso.Add(hist_data_ss_antiiso,hist_dy_ss_antiiso,1.,-1.)
    hist_data_ss_antiiso.Add(hist_data_ss_antiiso,hist_wj_ss_antiiso,1.,-1.)
    
    # normalization of the DY background
    data_lowmt = hist_data_lowmt.GetBinContent(1)
    data_lowmt_err = hist_data_lowmt.GetBinError(1)
    data_lowmt_rel = data_lowmt_err/data_lowmt
    dy_lowmt = hist_dy.GetBinContent(1)
    dy_lowmt_err = hist_dy.GetBinError(1)
    dy_lowmt_rel = dy_lowmt_err/dy_lowmt
    scaleDY = data_lowmt/dy_lowmt
    scaleDY_rel = math.sqrt(data_lowmt_rel*data_lowmt_rel+dy_lowmt_rel*dy_lowmt_rel)
    scaleDY_err = scaleDY * scaleDY_rel
    
    # normalization of the WJ background
    data_highmt = hist_data_highmt.GetBinContent(2)
    data_highmt_err = hist_data_highmt.GetBinError(2)
    data_highmt_rel = data_highmt_err/data_highmt
    wj_highmt = hist_wj.GetBinContent(2)
    wj_highmt_err = hist_wj.GetBinError(2)
    wj_highmt_rel = wj_highmt_err/wj_highmt
    scaleWJ = data_highmt/wj_highmt
    scaleWJ_rel = math.sqrt(data_highmt_rel*data_highmt_rel+wj_highmt_rel*wj_highmt_rel)
    scaleWJ_err = scaleWJ * scaleWJ_rel

    # QCD SS->OS extrapolation factor
    data_os = hist_data_os_antiiso.GetBinContent(1)
    data_os_err = hist_data_os_antiiso.GetBinError(1)
    data_os_rel = data_os_err/data_os
    data_ss = hist_data_ss_antiiso.GetBinContent(1)
    data_ss_err = hist_data_ss_antiiso.GetBinError(1)
    data_ss_rel = data_ss_err/data_ss
    scaleQCD = data_os/data_ss
    scaleQCD_rel = math.sqrt(data_os_rel*data_os_rel+data_ss_rel*data_ss_rel)
    scaleQCD_err = scaleQCD * scaleQCD_rel
    
    nbins = hist_data.GetNbinsX()
    xbins = []
    for ib in range(1,nbins+2):
        xbins.append(hist_data.GetXaxis().GetBinLowEdge(ib))

    print('')
    print('mT distribution used for calibration')
    print('nbins=%1i'%(nbins),xbins)
    print('')
    print('Calibration coefficients ->')
    print('WJ  = %5.3f +/- %5.3f'%(scaleWJ,scaleWJ_err))
    print('DY  = %5.3f +/- %5.3f'%(scaleDY,scaleDY_err))
    print('QCD = %5.3f +/- %5.3f'%(scaleQCD,scaleQCD_err))
    print('')
            
def ExtractDatacards(f,var,bins,isBDT,**kwargs):

    xmin = kwargs.get('xmin',0.)
    xmax = kwargs.get('xmax',1.)
    bdtbin = kwargs.get('bdtbin','')
    scaleWJ = kwargs.get('scaleWJ',1.0)
    scaleDY = kwargs.get('scaleDY',1.0)
    scaleQCD = kwargs.get('scaleQCD',1.0)
    era = kwargs.get('era','Run3')
    
    signal_samples = ['ggH_sm','ggH_ps','ggH_mm','qqH','HWminus','HWplus','ZH']
    bkg_samples = ['zll_incl','zll_ext','zll_0j','zll_1j','zll_2j','ztt_0j','ztt_1j','ztt_2j','top','vv','wjets']
    nowjets_samples = ['zll_incl','zll_ext','zll_0j','zll_1j','zll_2j','ztt_0j','ztt_1j','ztt_2j','top','vv']
    dy_samples = ['zll_incl','zll_ext','zll_0j','zll_1j','zll_2j','ztt_0j','ztt_1j','ztt_2j']
    if era in ['Run3_2023','Run3_2023preBPix','Run3_2023postBPix']:
        bkg_samples.remove('zll_ext')
        nowjets_samples.remove('zll_ext')
        dy_samples.remove('zll_ext')

    sign_labels = ['os','ss']
    iso_labels = ['iso','antiiso']
    typetau_labels = ['faketau','hadtau','leptau']

    decay_hypotheses = ['sm','ps','mm','flat']
    prod_hypotheses = ['sm','ps','mm']
    procs = ['qqH','HWplus','HWminus','ZH','ggH']

    hists = {}
    for sample in bkg_samples: 
        for sign in sign_labels:
            for iso in iso_labels:
                for typetau in typetau_labels:
                    name = '%s_%s_%s_%s_all_%s'%(sample,var,sign,iso,typetau)
                    histo = f.Get(name)
                    if isBDT:
                        hists[name] = utils.rebinHisto(histo,bins,'%s_rebinned'%(var))
                    else:
                        hists[name] = utils.slice2DHisto(histo,xmin,xmax,bins,'%s_%s_rebinned'%(var,bdtbin))
                

    for sample in signal_samples: 
        for sign in sign_labels:
            for typetau in typetau_labels:
                name = '%s_%s_%s_iso_all_%s'%(sample,var,sign,typetau)
                histo = f.Get(name)
                if isBDT:
                    hists[name] = utils.rebinHisto(histo,bins,'%s_rebinned'%(var))
                else:
                    hists[name] = utils.slice2DHisto(histo,xmin,xmax,bins,'%s_%s_rebinned'%(var,bdtbin))
                for hypothesis in prod_hypotheses:
                    name = '%s_%s_%s_iso_all_%s_%s'%(sample,var,sign,typetau,hypothesis)
                    histo = f.Get(name)
                    if isBDT:
                        hists[name] = utils.rebinHisto(histo,bins,'%s_rebinned'%(var))
                    else:
                        hists[name] = utils.slice2DHisto(histo,xmin,xmax,bins,'%s_%s_rebinned'%(var,bdtbin))

    histograms = {}
    name_data = 'data_%s_os_iso_all_all'%(var)
    histo_data = f.Get(name_data)
    if isBDT:
        histograms['data_obs'] = utils.rebinHisto(histo_data,bins,'%s_rebinned'%(var))
    else:
        histograms['data_obs'] = utils.slice2DHisto(histo_data,xmin,xmax,bins,'%s_%s_rebinned'%(var,bdtbin))

    basehist = hists['zll_incl_'+var+'_os_iso_all_hadtau']

    # JetFakes ->

    # QCD estimate from SS region ->
    name_data_ss = 'data_%s_ss_iso_all_all'%(var) 
    histo_data_ss = f.Get(name_data_ss)
    if isBDT:
        histograms['JetFakes'] = utils.rebinHisto(histo_data_ss,bins,'%s_rebinned'%(var))
    else:
        histograms['JetFakes'] = utils.slice2DHisto(histo_data_ss,xmin,xmax,bins,'%s_%s_rebinned'%(var,bdtbin))
    mc_ss = utils.ReplicaHist(basehist,'MC_ss_'+var)
    for sample in bkg_samples:
        for typetau in typetau_labels:
            name = '%s_%s_ss_iso_all_%s'%(sample,var,typetau)
            mc_ss.Add(mc_ss,hists[name],1.,1.)
    histograms['JetFakes'].Add(histograms['JetFakes'],mc_ss,1.,-1.)
    histograms['JetFakes'].Scale(scaleQCD)
    
    # summing up non-QCD fakes
    for sample in nowjets_samples:
        name = '%s_%s_os_iso_all_faketau'%(sample,var)
        histograms['JetFakes'].Add(histograms['JetFakes'],hists[name],1.,1.)
    for typetau in typetau_labels:    
        name = 'wjets_%s_os_iso_all_%s'%(var,typetau)
        histograms['JetFakes'].Add(histograms['JetFakes'],hists[name],1.,scaleWJ)
    
    histograms['ZL'] = utils.ReplicaHist(basehist,'ZL_'+var+'_'+bdtbin)
    histograms['ZTT'] = utils.ReplicaHist(basehist,'ZTT_'+var+'_'+bdtbin)
    histograms['VVT'] = utils.ReplicaHist(basehist,'VVT_'+var+'_'+bdtbin)
    histograms['TTT'] = utils.ReplicaHist(basehist,'TTT_'+var+'_'+bdtbin)

    histograms['ZL'].Scale(scaleDY)
    histograms['ZTT'].Scale(scaleDY)
    
    for sample in dy_samples:
        name = '%s_%s_os_iso_all_leptau'%(sample,var)
        histograms['ZL'].Add(histograms['ZL'],hists[name],1.,scaleDY)
        name = '%s_%s_os_iso_all_hadtau'%(sample,var)
        histograms['ZTT'].Add(histograms['ZTT'],hists[name],1.,scaleDY)

    for tautype in ['leptau','hadtau']:
        name = 'vv_%s_os_iso_all_%s'%(var,tautype)
        histograms['VVT'].Add(histograms['VVT'],hists[name],1.,1.)
        name = 'top_%s_os_iso_all_%s'%(var,tautype)
        histograms['TTT'].Add(histograms['TTT'],hists[name],1.,1.)

    
    for proc in procs:
        for decay_hypothesis in decay_hypotheses:
            decay_name = '_%s'%(decay_hypothesis)
            if decay_hypothesis=='flat':
                decay_name = ''
            if proc=='ggH':
                for prod_hypothesis in prod_hypotheses:
                    name = 'ggH_%s_%s_os_iso_all_hadtau%s'%(prod_hypothesis,var,decay_name)
                    template_name = 'ggH_%s_prod_%s_htt125'%(decay_hypothesis,prod_hypothesis)                    
                    histograms[template_name] = hists[name]
            else:
                name = '%s_%s_os_iso_all_hadtau%s'%(proc,var,decay_name)
                template_name = '%s_%s_htt125'%(proc,decay_hypothesis)
                histograms[template_name] = hists[name]

    for decay_hypothesis in decay_hypotheses:
        histograms['WH_%s_htt125'%(decay_hypothesis)] = histograms['HWplus_%s_htt125'%(decay_hypothesis)]
        histograms['WH_%s_htt125'%(decay_hypothesis)].Add(histograms['WH_%s_htt125'%(decay_hypothesis)],
                                                          histograms['HWminus_%s_htt125'%(decay_hypothesis)],1.,1.)
        
    return histograms
                
                
    
if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3',choices=['Run3_2022preEE','Run3_2022postEE','Run3_2022','Run3_2023preBPix','Run3_2023postBPix','Run3_2023','Run3'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',action='store_true')
    parser.add_argument('-bVeto','--bVeto', dest='bVeto',action='store_true')
    parser.add_argument('-wjNorm','--wjNorm',dest='wjNorm',type=float,default=1.0)
    parser.add_argument('-dyNorm','--dyNorm',dest='dyNorm',type=float,default=1.0)
    parser.add_argument('-qcdNorm','--qcdNorm',dest='qcdNorm',type=float,default=1.0)
    parser.add_argument('-runCalibr','--runCalibr',dest='runCalibr',action='store_true')
    parser.add_argument('-dryRun','--dryRun',dest='dryRun',action='store_true')
    args = parser.parse_args()

    era = args.era
    useCrossTrigger = args.useCrossTrigger
    bVeto = args.bVeto
    chan = args.channel
    runCalibr = args.runCalibr
    dryRun = args.dryRun
    
    scaleDY = args.dyNorm
    scaleWJ = args.wjNorm
    scaleQCD = args.qcdNorm
    
    bdt_et_bins = {
        'ditau' : [0.33,0.50,0.60,0.70,0.80,0.90,1.0],
        'fakes' : [0.33,0.50,0.60,0.70,0.80,0.90,1.0],
    }
    bdt_mt_bins = {
        'ditau' : [0.33,0.50,0.60,0.70,0.80,0.9,1.0],
        'fakes' : [0.33,0.50,0.60,0.70,0.80,0.9,1.0],
    }

    bdt_bins = bdt_mt_bins
    if chan=='et':
        bdt_bins = bdt_et_bins
    
    cats = {
        'signal': 'higgs',
        'ditau' : 'tau',
        'fakes' : 'fake'
    }

    chan_lep = {
        'mt' : 'mu',
        'et' : 'e',
    }

    modes_bins = {
        'lep_pi' : 8,
        'lep_rho' : 10,
        'lep_a1_1pr' : 8,
        'lep_a1_3pr' : 10,
    }
    
    modes_map = {
        'lep_pi' : 'pi',
        'lep_rho' : 'rho',
        'lep_a1_1pr' : 'a11pr',
        'lep_a1_3pr' : 'a1',
    }
    
    phicp_min = 0.
    phicp_max = 360.

    suffix_xtrig = ''
    suffix_bveto = ''
    if useCrossTrigger:
        suffix_xtrig = '_xtrig'
    if bVeto:
        suffix_bveto = '_bveto'
    
    suffix = '_x_mtcut'+suffix_bveto+suffix_xtrig+'_ipcut1_promptSF_tauSF'
    basedir = utils.outputFolder
    inputFileName = '%s/selection/datacardsPhiCP/%s_%s%s.root'%(basedir,chan,era,suffix)
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
    hists_mt = utils.extractHistos(inputFile,'mt_1',binsMT,'amcatnlo',era,chan)
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
