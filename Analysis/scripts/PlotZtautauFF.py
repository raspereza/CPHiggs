#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def ExtractHistos(f,var,bins,region):

    hists = {}
    for sample in utils.samples:
        for typ in utils.type_labels:
            nameInput='%s_%s_%s_%s'%(sample,var,region,typ)
            name='%s_%s'%(sample,typ)
            hists[name] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')
            for ff in utils.ff_labels:
                nameInput='%s_%s_%s_%s_%s'%(sample,var,region,ff,typ)
                name='%s_%s_%s'%(sample,ff,typ)
                hists[name] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')

    return hists
    

def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    suffixIP = kwargs.get('suffixIP','ipcut')
    region = kwargs.get('region','lowmt_os_iso')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.501)
    ymax = kwargs.get('ymax',1.499)
    
    mc_samples = ['ztt','zll','top','vv','wjets']
    h_mc = {} # sample, {sig,ar,qcd,wj}, {jet, lepton} 
    # histograms
    
    h_data = hists['data_all'].Clone('h_data')
    for mc in mc_samples:
        # MC with gen_match_2 = [1,5]
        name = '%s_lepton'%(mc)
        h_mc[name] = hists['%s_lep'%(mc)].Clone('h_'+name)
        h_mc[name].Add(h_mc[name],hists['%s_tau'%(mc)],1.,1.)
        # MC with gen_match_2 = 0
        name = '%s_had'%(mc)
        h_mc[name] = hists['%s_had'%(mc)].Clone('h_'+name)

        for ff in utils.ff_labels:
            
            name = '%s_%s_lepton'%(mc,ff)
            h_mc[name] = hists['%s_%s_lep'%(mc,ff)].Clone('h_'+name)
            h_mc[name].Add(h_mc[name],hists['%s_%s_tau'%(mc,ff)],1.,1.)
            
            name = '%s_%s_had'%(mc,ff)
            h_mc[name] = hists['%s_%s_had'%(mc,ff)].Clone('h_'+name)
            
    h_data_qcd = hists['data_qcd_closure_all'].Clone('h_data_qcd')
    h_data_ar = hists['data_ar_all'].Clone('h_data_ar')
    h_data_ewk = hists['data_wj_all'].Clone('h_data_ewk')
    h_data_top = hists['data_mc_top_all'].Clone('h_data_top')
    
    for mc in mc_samples:
        hist = h_mc['%s_ar_lepton'%(mc)]
        h_data_ar.Add(h_data_ar,hist,1.,-1.)
        hist = h_mc['%s_qcd_closure_lepton'%(mc)]
        h_data_qcd.Add(h_data_qcd,hist,1.,-1.)
        hist = h_mc['%s_wj_lepton'%(mc)]
        h_data_ewk.Add(h_data_ewk,hist,1.,-1.)
        hist = h_mc['%s_mc_top_lepton'%(mc)]
        h_data_top.Add(h_data_top,hist,1.,-1.)

    ewk_ar = h_mc['wjets_ar_had']
    for mc in ['zll','ztt','vv']:
        ewk_ar.Add(ewk_ar,h_mc['%s_ar_had'%(mc)],1.,1.)        
    top_ar = h_mc['top_ar_had']
        
    nbins = h_data_ar.GetNbinsX() 
    for ib in range(1,nbins+1):
        x_data_ar = max(h_data_ar.GetBinContent(ib),1.)
        x_ewk_ar = ewk_ar.GetBinContent(ib)
        x_top_ar = top_ar.GetBinContent(ib)
        if x_ewk_ar<0: x_ewk_ar = 0.
        f_ewk = x_ewk_ar/x_data_ar
        if f_ewk>0.99: f_ewk = 0.99
        f_top = x_top_ar/x_data_ar
        f_mc = f_top + f_ewk
        if f_mc>0.99:
            f_top = 0.99 - f_ewk
            f_mc = 0.99
        f_qcd = 1.0-f_mc
        xlower = h_data_ar.GetXaxis().GetBinLowEdge(ib)
        xupper = h_data_ar.GetYaxis().GetBinLowEdge(ib+1)
        f_tot = f_top+f_ewk+f_qcd
        print('[%5.2f  %5.2f]   %5.3f  %5.3f  %5.3f  %5.3f'%(xlower,xupper,f_qcd,f_ewk,f_top,f_tot))
        x_ff_qcd = h_data_qcd.GetBinContent(ib)
        x_ff_ewk = h_data_ewk.GetBinContent(ib)
        x_ff_top = h_data_top.GetBinContent(ib) 
        x_ff = f_ewk*x_ff_ewk+f_qcd*x_ff_qcd#+f_top*x_ff_top
#        e_ff = ROOT.TMath.Sqrt(f_ewk*e_ff_ewk*f_ewk*e_ff_ewk+f_qcd*e_ff_qcd*f_qcd*e_ff_qcd)
        h_data_qcd.SetBinContent(ib,x_ff)
#        h_data_qcd.SetBinError(ib,e_ff)
        
    styles.InitData(h_data)

    h_ztt = h_mc['ztt_lepton'].Clone('ztt')
    h_zll = h_mc['zll_lepton'].Clone('zll')
    h_top = h_mc['top_lepton'].Clone('top')
    h_vv = h_mc['vv_lepton'].Clone('vv')
    h_wjets = h_mc['wjets_lepton'].Clone('wjets')
    h_qcd = h_data_qcd.Clone('qcd')
    
    xtitle = utils.XTitle[chan][var]
    styles.InitHist(h_ztt,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_vv,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    styles.InitHist(h_qcd,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)
    styles.InitHist(h_wjets,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)

    x_data = h_data.GetSumOfWeights()
    
    x_ztt   = h_ztt.GetSumOfWeights()
    x_zll   = h_zll.GetSumOfWeights() 
    x_top   = h_top.GetSumOfWeights()
    x_vv    = h_vv.GetSumOfWeights() 
    x_wjets = h_wjets.GetSumOfWeights()
    x_qcd   = h_qcd.GetSumOfWeights() 
    x_tot   = x_ztt + x_zll + x_top + x_vv + x_wjets + x_qcd

    print('')
    print('Yields ->')
    print('Ztautau    : %7.0f'%(x_ztt))
    print('Zll        : %7.0f'%(x_zll))
    print('TTbar      : %7.0f'%(x_top))
    print('VV+ST      : %7.0f'%(x_vv))
    print('WJets      : %7.0f'%(x_wjets))
    print('Fakes      : %7.0f'%(x_qcd))
    print('Total MC   : %7.0f'%(x_tot))
    print('Data       : %7.0f'%(x_data))
    print('')

    h_vv.Add(h_vv,h_top,1.,1.)
    h_vv.Add(h_vv,h_wjets,1.,1.)
    h_qcd.Add(h_qcd,h_vv,1.,1.)
    h_zll.Add(h_zll,h_qcd,1.,1.)
    h_ztt.Add(h_ztt,h_zll,1.,1.)

    nbins = h_ztt.GetNbinsX()
    h_tot = h_ztt.Clone("total")
    print('Total (cross-check) : %7.0f'%(h_tot.GetSumOfWeights()))
    print('')
    
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_ztt)
    utils.zeroBinErrors(h_zll)
    utils.zeroBinErrors(h_top)
    utils.zeroBinErrors(h_vv)
    utils.zeroBinErrors(h_wjets)
    utils.zeroBinErrors(h_qcd)

    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    
    h_data.GetYaxis().SetRangeUser(0.,1.2*YMax)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("events")
    h_ratio.GetYaxis().SetTitle("obs/exp")
    h_ratio.GetXaxis().SetTitle(xtitle)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    
    
    h_data.Draw('e1')
    h_ztt.Draw('hsame')
    h_zll.Draw('hsame')
    h_qcd.Draw('hsame')
    h_vv.Draw('hsame')
    h_top.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.65,0.4,0.85,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.042)
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_ztt,'Z#rightarrow#tau#tau','f')
    if chan=='mt': leg.AddEntry(h_zll,'Z#rightarrow#mu#mu','f')
    else: leg.AddEntry(h_zll,'Z#rightarrowee','f')
    leg.AddEntry(h_qcd,'j#rightarrow#tau fakes','f')
    leg.AddEntry(h_vv,'electroweak','f')
    leg.AddEntry(h_top,'t#bar{t}','f')
    if plotLegend: leg.Draw()

    styles.CMS_label(upper,era=era)

    upper.Draw("SAME")
    upper.RedrawAxis()
    upper.Modified()
    upper.Update()
    canvas.cd()

    # lower pad
    lower = ROOT.TPad("lower", "pad",0,0,1,0.30)
    lower.Draw()
    lower.cd()
    styles.InitLowerPad(lower)

    h_ratio.Draw('e1')
    h_tot_ratio.Draw('e2same')
    h_ratio.Draw('e1same')

    nbins = h_ratio.GetNbinsX()
    xmin = h_ratio.GetXaxis().GetBinLowEdge(1)    
    xmax = h_ratio.GetXaxis().GetBinLowEdge(nbins+1)
    line = ROOT.TLine(xmin,1.,xmax,1.)
    line.SetLineStyle(1)
    line.SetLineWidth(2)
    line.SetLineColor(4)
    line.Draw()
    lower.Modified()
    lower.RedrawAxis()

    canvas.cd()
    canvas.Modified()
    canvas.cd()
    canvas.SetSelected(canvas)
    canvas.Update()
    print('')
    outputFolder = '/eos/home-r/rasp/php-plots/plots/%s_FF/%s'%(chan,suffixIP)
    outputGraphics = '%s/%s.png'%(outputFolder,var)
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3',choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et'])
    parser.add_argument('-variable','--variable',dest='variable',default='m_vis')
    parser.add_argument('-nbins','--nbins', dest='nbins',type=int,default=50)
    parser.add_argument('-xmin','--xmin',dest='xmin',type=float,default=0.0)
    parser.add_argument('-xmax','--xmax',dest='xmax',type=float,default=250.)
    parser.add_argument('-ymin','--ymin',dest='ymin',type=float,default=0.501)
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.499)
    parser.add_argument('-suffix','--suffix',dest='suffix',default='x_ipcut1_ff_v3')
    parser.add_argument('-region','--region',dest='region',default='lowmt_os_iso')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    region = args.region
    suffix = args.suffix
    
    suffixIP = 'ipcut'
        
    plotLegend = True
    if var in ['eta_1','eta_2','bdt_ditau','bdt_fakes','bdt_signal']:
        plotLegend = False

    basedir = utils.outputFolder+'/selection/baseline'
    bins = utils.createBins(nbins,xmin,xmax)

    inputFileName = '%s/%s_%s_%s.root'%(basedir,chan,era,suffix)
    if os.path.isfile(inputFileName):
        print('')
        print('Loading ROOT file %s'%(inputFileName))
        print('')
    else:
        print('')
        print('ROOT file %s not found'%(inputFileName))
        print('quitting')
        print('')
    inputFile = ROOT.TFile(inputFileName,'read')
    hists = ExtractHistos(inputFile,var,bins,region)
    Plot(hists,era=era,var=var,channel=chan,
         suffixIP=suffixIP,ymin=ymin,ymax=ymax,
         plotLegend=plotLegend,region=region)
    

