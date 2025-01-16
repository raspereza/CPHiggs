#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

XTitle = {
    'mm': {
        'mt_1'  : "m_{T} (GeV)",
        'pt_1'  : "leading #mu p_{T} (GeV)",
        'eta_1' : "leading #mu #eta",
        'pt_2'  : "trailing #mu p_{T} (GeV)",
        'eta_2' : "trailing #mu #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{#mu#mu} (GeV)",
        'ipsig_1': "leading #mu IP sig",
        'ipsig_2': "trailing #mu IP sig"
    },
    'ee': {
        'mt_1'  : "m_{T} (GeV)",
        'pt_1'  : "leading elec p_{T} (GeV)",
        'eta_1' : "leading elec #eta",
        'pt_2'  : "trailing elec p_{T} (GeV)",
        'eta_2' : "trailing elec #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{ee} (GeV)",
	'ipsig_1': "leading elec IP sig",
        'ipsig_2': "trailing elec IP sig"
    }
}

def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.5)
    ymax = kwargs.get('ymax',1.5)
    
    # histograms
    h_data = hists['data_'+var+'_os_iso_all'].Clone('h_data')
    h_zll = hists['dy_'+var+'_os_iso_all'].Clone('h_zll')
    h_top = hists['top_'+var+'_os_iso_all'].Clone('h_top')
    h_vv = hists['vv_'+var+'_os_iso_all'].Clone('h_vv')
    h_wjets = hists['wjets_'+var+'_os_iso_all'].Clone('h_wjets')

    styles.InitData(h_data)

    xtitle = XTitle[chan][var]
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_vv,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    styles.InitHist(h_wjets,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    
    x_data = h_data.GetSumOfWeights()
    
    x_zll   = h_zll.GetSumOfWeights() 
    x_top   = h_top.GetSumOfWeights()
    x_vv    = h_vv.GetSumOfWeights() 
    x_wjets = h_wjets.GetSumOfWeights()
    x_tot   = x_zll + x_top + x_vv + x_wjets

    print('')
    print('Yields ->')
    print('DY    : %7.0f'%(x_zll))
    print('Top   : %7.0f'%(x_top))
    print('VV    : %7.0f'%(x_vv))
    print('WJets : %7.0f'%(x_wjets))
        
    print('Total : %7.0f'%(x_tot))
    print('Data  : %7.0f'%(x_data))
    print('')

    h_vv.Add(h_vv,h_wjets,1.,1.)
    h_top.Add(h_top,h_vv,1.,1.)
    h_zll.Add(h_zll,h_top,1.,1.)

    h_tot = h_zll.Clone("total")
    styles.InitTotalHist(h_tot)

    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_zll)
    utils.zeroBinErrors(h_top)
    utils.zeroBinErrors(h_vv)
    utils.zeroBinErrors(h_wjets)
    
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
    h_zll.Draw('hsame')
    h_top.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.7,0.45,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_zll,'Drell-Yan','f')
    leg.AddEntry(h_top,'rest','f')
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
    outputGraphics = os.getenv('CMSSW_BASE') + '/src/CPHiggs/IP/figures/' + var + '_' + chan + '_' + era + '.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mm','ee'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='m_vis')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=40)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=50.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=130.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    parser.add_argument('-applyIP1','--applyIP1',dest='applyIP1',type=int,default=0)
    parser.add_argument('-applyIP2','--applyIP2',dest='applyIP2',type=int,default=0)
    parser.add_argument('-applySF', '--applySF', dest='applySF' ,type=int,default=0)

    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax

    plotLegend = True
    if var=='eta_1' or var=='eta_2':
        plotLegend = False

    
    applyIPSigLep1Cut = args.applyIP1
    applyIPSigLep2Cut = args.applyIP2
    applyIPSigPromptLepSF = args.applySF
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    basedir = os.getenv('CMSSW_BASE')+'/src/CPHiggs/IP/selection'

    suffix_ip1 = ''
    suffix_ip2 = ''
    suffix_prompt = ''
    if applyIPSigLep1Cut==1:
        suffix_ip1 = '_ipcut1'
    if applyIPSigLep2Cut==1:
        suffix_ip2 = '_ipcut2'
    if applyIPSigPromptLepSF==1:
        suffix_prompt = '_promptSF'

    suffix = suffix_ip1+suffix_ip2+suffix_prompt 
    
    inputFileName = '%s/%s_%s%s.root'%(basedir,chan,era,suffix)
    if os.path.isfile(inputFileName):        
        print('')
        print('Loading histograms from file %s'%(inputFileName))
        print('')
    else:
        print('')
        print('Input ROOT file %s not found'%(inputFileName))
        print('Quitting')
        print('')
        exit()
    inputFile = ROOT.TFile(inputFileName,'read')
    hists = utils.extractHistos(inputFile,var,bins)
    Plot(hists,era=era,var=var,channel=chan,ymin=ymin,ymax=ymax,plotLegend=plotLegend)

    
