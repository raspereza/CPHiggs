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
        for sign in utils.sign_labels:
            nameInput='%s_%s_%s_%s_iso_all'%(sample,var,region,sign)
            name='%s_%s'%(sample,sign)
            hists[name] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')

    return hists
    
def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    suffix = kwargs.get('suffix','ipcut')
    qcdNorm = kwargs.get('qcdNorm',1.20)
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.701)
    ymax = kwargs.get('ymax',1.299)
    blind_data = kwargs.get('blind_data',True)
    
    # histograms
    h_data = hists['data_os'].Clone('h_data')
    h_ztt = hists['ztt_os'].Clone('h_ztt')
    h_zll = hists['zll_os'].Clone('h_zll')
    h_top = hists['top_os'].Clone('h_top')
    h_vv = hists['vv_os'].Clone('h_vv')
    h_wjets = hists['wjets_os'].Clone('h_wjets')
    
    # qcd estimation
    h_qcd = hists['data_ss'].Clone('h_qcd')
    for mc_sample in ['zll','ztt','top','vv','wjets']:
        name = '%s_ss'%(mc_sample)
        h_qcd.Add(h_qcd,hists[name],1.,-1.) 

    h_qcd.Scale(qcdNorm)

    styles.InitData(h_data)

    xtitle = utils.XTitle[chan][var]
    styles.InitHist(h_ztt,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_vv,"","",ROOT.kCyan,1001)
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
    print('QCD        : %7.0f'%(x_qcd))
    print('Total MC   : %7.0f'%(x_tot))
    print('Data       : %7.0f'%(x_data))
    print('')

    h_qcd.Add(h_qcd,h_vv,1.,1.)
    h_wjets.Add(h_wjets,h_qcd,1.,1.)
    h_top.Add(h_top,h_wjets,1.,1.)
    h_zll.Add(h_zll,h_top,1.,1.)
    h_ztt.Add(h_ztt,h_zll,1.,1.)

    nbins = h_ztt.GetNbinsX()
    h_tot = h_ztt.Clone("total")
    print('Total (cross-check) : %7.0f'%(h_tot.GetSumOfWeights()))
    print('')

    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    
    styles.InitTotalHist(h_tot)

    if blind_data and var=='bdt_signal':
        for ib in range(1,nbins+1):
            bdt = h_data.GetXaxis().GetBinLowEdge(ib)
            if bdt>0.6:
                h_data.SetBinContent(ib,10000.)
                h_data.SetBinError(ib,0.)
    
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
    h_top.Draw('hsame')
    h_wjets.Draw('hsame')
    h_qcd.Draw('hsame')
    h_vv.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.65,0.4,0.85,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.042)
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_ztt,'Z#rightarrow#tau#tau','f')
    if chan=='mt': leg.AddEntry(h_zll,'Z#rightarrow#mu#mu','f')
    else: leg.AddEntry(h_zll,'Z#rightarrowee','f')
    leg.AddEntry(h_top,'t#bar{t}','f')
    leg.AddEntry(h_wjets,'W+jets','f')
    leg.AddEntry(h_qcd,'QCD','f')
    leg.AddEntry(h_vv,'electroweak','f')
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
    fig_folder = '/eos/home-r/rasp/php-plots/plots/%s_classic/'%(chan)
    outputGraphics = '%s/%s.png'%(fig_folder,var)
    canvas.Print(outputGraphics)
    return h_data, h_ztt

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3',choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-variable','--variable',dest='variable',default='m_vis')
    parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et'])
    parser.add_argument('-nbins','--nbins', dest='nbins',type=int,default=50)
    parser.add_argument('-xmin','--xmin',dest='xmin',type=float,default=0.0)
    parser.add_argument('-xmax','--xmax',dest='xmax',type=float,default=250.)
    parser.add_argument('-ymin','--ymin',dest='ymin',type=float,default=0.501)
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.499)
    parser.add_argument('-suffix','--suffix',dest='suffix',default='x_ipcut1_ff_noipcut')
    parser.add_argument('-region','--region',dest='region',default='lowmt')
    parser.add_argument('-qcdNorm','--qcdNorm',dest='qcdNorm',type=float,default=1.0)
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    suffix = args.suffix
    region = args.region
    qcdNorm = args.qcdNorm

    suffixIP = 'no_ipcut'
    if 'ipcut1' in suffix:
        suffixIP = 'ipcut'
    
    plotLegend = True
    if var in ['eta_1','eta_2','bdt_ditau','bdt_fakes']:
        plotLegend = False

    basedir = utils.outputFolder+'/selection'
    bins = utils.createBins(nbins,xmin,xmax)

    inputFileName = '%s/baseline/%s_%s_%s.root'%(basedir,chan,era,suffix)
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
         suffix=suffixIP,ymin=ymin,ymax=ymax,
         plotLegend=plotLegend,qcdNorm=qcdNorm)
    

