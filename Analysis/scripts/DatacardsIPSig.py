#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.utils as utils

def CreateCards(hists,**kwargs):

    era = kwargs.get('era','Run3_2022')
    chan = kwargs.get('channel','mt')
    binPt = kwargs.get('binPt','1')
    binEta = kwargs.get('binEta','1')
    region = kwargs.get('region','pass')

    binPtEta = binPt+'_'+binEta
    # histograms
    h_data = hists['data_m_vis_'+region+'_'+binPtEta+'_os_iso_all'].Clone('h_data')

    # genuine taus ->
    h_ztt = hists['dy_m_vis_'+region+'_'+binPtEta+'_os_iso_tau'].Clone('h_ztt')
    h_top = hists['top_m_vis_'+region+'_'+binPtEta+'_os_iso_tau'].Clone('h_top_tau')    
    h_vv = hists['vv_m_vis_'+region+'_'+binPtEta+'_os_iso_tau'].Clone('h_vv_tau')    

    h_wj = hists['wjets_m_vis_pass_'+binPtEta+'_os_iso_all'].Clone('h_wj')
    
    # prompt leptons ->
    h_zll = hists['dy_m_vis_pass_'+binPtEta+'_os_iso_lep'].Clone('h_zll')
    h_zll.Add(h_zll,hists['dy_m_vis_pass_'+binPtEta+'_os_iso_had'],1.,1.)

    h_top_lep = hists['top_m_vis_pass_'+binPtEta+'_os_iso_lep'].Clone('h_top_lep')
    h_top_lep.Add(h_top_lep,hists['top_m_vis_pass_'+binPtEta+'_os_iso_had'],1.,1.)
    
    h_vv_lep = hists['vv_m_vis_pass_'+binPtEta+'_os_iso_lep'].Clone('h_vv_lep')
    h_vv_lep.Add(h_vv_lep,hists['vv_m_vis_pass_'+binPtEta+'_os_iso_had'],1.,1.)


    if region=='fail':

        h_zll_incl = hists['dy_m_vis_incl_'+binPtEta+'_os_iso_lep']
        h_zll_incl.Add(h_zll_incl,hists['dy_m_vis_incl_'+binPtEta+'_os_iso_had'],1.,1.)

        h_top_incl = hists['top_m_vis_incl_'+binPtEta+'_os_iso_lep']
        h_top_incl.Add(h_top_incl,hists['top_m_vis_incl_'+binPtEta+'_os_iso_had'],1.,1.)

        h_vv_incl = hists['vv_m_vis_incl_'+binPtEta+'_os_iso_lep']
        h_vv_incl.Add(h_vv_incl,hists['vv_m_vis_incl_'+binPtEta+'_os_iso_had'],1.,1.)
        
        h_wj_incl = hists['wjets_m_vis_incl_'+binPtEta+'_os_iso_all']
        
        h_zll.Add(h_zll_incl,h_zll,1.,-1.)
        h_top_lep.Add(h_top_incl,h_top_lep,1.,-1.)
        h_vv_lep.Add(h_vv_incl,h_vv_lep,1.,-1.)
        h_wj.Add(h_wj_incl,h_wj,1.,-1.)

    # qcd estimation
    h_qcd = hists['data_m_vis_'+region+'_'+binPtEta+'_ss_iso_all'].Clone('h_qcd')
    for mc_sample in ['dy','top','vv','wjets']:
        if region=='pass':
            h_qcd.Add(h_qcd,hists[mc_sample+'_m_vis_pass_'+binPtEta+'_ss_iso_all'],1.,-1.)
        else:
            h_qcd.Add(h_qcd,hists[mc_sample+'_m_vis_pass_'+binPtEta+'_ss_iso_all'],1.,1.)
            h_qcd.Add(h_qcd,hists[mc_sample+'_m_vis_incl_'+binPtEta+'_ss_iso_all'],1.,-1.)


    utils.removeNegativeBins(h_ztt)
    utils.removeNegativeBins(h_zll)
    utils.removeNegativeBins(h_top)
    utils.removeNegativeBins(h_vv)
    utils.removeNegativeBins(h_wj)
    utils.removeNegativeBins(h_top_lep)
    utils.removeNegativeBins(h_vv_lep)
    utils.removeNegativeBins(h_qcd)

    x_data = h_data.GetSumOfWeights()
    x_ztt  = h_ztt.GetSumOfWeights()
    x_zll  = h_zll.GetSumOfWeights()
    x_ttt  = h_top.GetSumOfWeights()
    x_vvt  = h_vv.GetSumOfWeights()
    x_wj   = h_wj.GetSumOfWeights()
    x_ttl  = h_top_lep.GetSumOfWeights()
    x_vvl  = h_vv_lep.GetSumOfWeights()
    x_qcd  = h_qcd.GetSumOfWeights()
    x_tot  = x_ztt+x_zll+x_ttt+x_vvt+x_wj+x_ttl+x_vvl+x_qcd

    print('')
    print('%s  %s  :  binPt%s_binEta%s  :  %s'%(era,chan,binPt,binEta,region))
    print('ZTT   = %7.0f'%(x_ztt))
    print('ZLL   = %7.0f'%(x_zll))
    print('TTT   = %7.0f'%(x_ttt))
    print('TTL   = %7.0f'%(x_ttl))
    print('VVT   = %7.0f'%(x_vvt))
    print('VVL   = %7.0f'%(x_vvl))
    print('WJ    = %7.0f'%(x_wj))
    print('QCD   = %7.0f'%(x_qcd))
    print('-----------------------------')
    print('Total = %7.0f'%(x_tot))
    print('Data  = %7.0f'%(x_data))
    print('')
    
    basedir = utils.outputFolder+'/datacards'

    suffix = 'binPt'+binPt+'_binEta'+binEta
    outputFileName = '%s/%s_%s_%s_%s.root'%(basedir,chan,era,suffix,region)
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    outputFile.mkdir(region)
    outputFile.cd(region)
    h_data.Write('data_obs')
    h_ztt.Write('ZTT_'+region)
    h_top.Write('TTT_'+region)
    h_vv.Write('VVT_'+region)
    h_qcd.Write('QCD')
    h_zll.Write('ZLL')
    h_top_lep.Write('TTL')
    h_vv_lep.Write('VVL')
    h_wj.Write('WJ')
    outputFile.Close()
    
    cardsFileName = '%s/%s_%s_%s_%s.txt'%(basedir,chan,era,suffix,region)
    f = open(cardsFileName,"w")
    f.write("imax 1  number of channels\n")
    f.write("jmax *  number of backgrounds\n")
    f.write("kmax *  number of nuisance parameters\n")
    f.write("---------------------------\n")
    f.write("observation %3.1f\n"%(x_data))
          
    f.write("---------------------------\n")
    f.write("shapes * * "+outputFileName+"  %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n"%(region,region))
    f.write("---------------------------\n")
    f.write("bin %s %s %s %s %s %s %s %s\n"%(region,region,region,region,region,region,region,region))
    f.write("process  TTT_%s  VVT_%s  ZTT_%s  TTL  VVL  ZLL  WJ  QCD\n"%(region,region,region))
    f.write("process              -2    -1     0     1     2     3     4     5\n")
    f.write("rate    %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n"%(x_ttt,x_vvt,x_ztt,x_ttl,x_vvl,x_zll,x_wj,x_qcd))
    f.write("---------------------------\n")
    f.write("lumi         lnN   1.02  1.02  1.02  1.02  1.02  1.02     -     -\n")
    f.write("xsec_tt      lnN   1.05     -     -  1.05     -     -     -     -\n")
    f.write("xsec_vv      lnN      -  1.07     -     -  1.07     -     -     -\n")
    f.write("xsec_dy      lnN      -     -  1.03     -     -  1.03     -     -\n")
    f.write("tau_id       lnN   1.05  1.05  1.05  1.05  1.05     -     -     -\n")
    f.write("lep_fakes    lnN     -      -     -     -     -  1.30     -     -\n")
    f.write("ipsig_prompt lnN     -      -     -  1.03  1.03  1.03  1.03     -\n")
    f.write("qcd_norm     lnN     -      -     -     -     -     -     -  1.30\n")
    f.write("jet_fakes  rateParam  %s  WJ   1.0  [0.5,1.5]\n"%(region))
    f.write("* autoMCStats 0\n")
    
if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022',choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='m_vis')
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',action='store_true')
    parser.add_argument('-applyMTCut','--applyMTCut', dest='applyMTCut', action='store_true')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=21)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=40.)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=250.)
    parser.add_argument('-generator','--generator',dest='generator',default='amcatnlo',choices=['amcatnlo','MG'])
    
    args = parser.parse_args()
    era = args.era
    applyMTCut = args.applyMTCut
    useCrossTrigger = args.useCrossTrigger
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    generator = args.generator
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    suffix_mt = ''
    suffix_xtrig = ''
    suffix_prompt = ''
    if useCrossTrigger:
        suffix_xtrig = '_xtrig'
    if applyMTCut:
        suffix_mt = '_mtcut'

    suffix = '_x'+suffix_mt+suffix_xtrig+'_promptSF'
    basedir = utils.outputFolder
    inputFileName = '%s/selection/%s_%s%s.root'%(basedir,chan,era,suffix)
    if os.path.isfile(inputFileName):
        print('')
        print('Loading ROOT file %s'%(inputFileName))
    else:
        print('')
        print('ROOT file %s is not found'%(inputFileName))
        print('Quitting')
        exit()
    inputFile = ROOT.TFile(inputFileName,'read')
    isSecond = False
    hists = utils.extractTagProbeHistos(inputFile,bins,generator,era,chan,isSecond)
    histPtBins = inputFile.Get('ptBins')
    histEtaBins = inputFile.Get('etaBins')
    nbinsPt = histPtBins.GetNbinsX()
    nbinsEta = histEtaBins.GetNbinsX()

    ptBins = []
    for iPt in range(1,nbinsPt+2):
        ptBins.append(histPtBins.GetBinLowEdge(iPt))

    etaBins = []
    for iEta in range(1,nbinsEta+2):
        etaBins.append(histEtaBins.GetBinLowEdge(iEta))

    print('')
    print('ptBins',ptBins)
    print('etaBins',etaBins)
    print('')

    region_labels = ['pass','fail']
    for iEta in range(1,nbinsEta+1):
        binEta = '%1i'%(iEta)
        for iPt in range(1,nbinsPt+1):
            binPt = '%1i'%(iPt)
            for region in region_labels:
                CreateCards(hists,era=era,channel=chan,binPt=binPt,binEta=binEta,region=region)
            
