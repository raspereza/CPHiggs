#! /usr/bin/env python3
# Author: Alexei Raspereza (October 2025)
# Selection of samples for measurement of the SF related to IPSig cut
import ROOT
import math
import yaml
from yaml.loader import SafeLoader
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

XTitle = {
    'mm': {
        'pt_1'  : "leading #mu p_{T} (GeV)",
        'eta_1' : "leading #mu #eta",
        'pt_2'  : "trailing #mu p_{T} (GeV)",
        'eta_2' : "trailing #mu #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{#mu#mu} (GeV)",
        'ipsig_1': "leading #mu IP sig",
        'ipsig_2': "trailing #mu IP sig",
        'jpt_1': "leading jet p_{T} (GeV)",
        'jpt_2': "trailing jet p_{T} (GeV)",
        'jeta_1': "leading jet #eta",
        'jeta_2': "trailing jet #eta",
        'n_jets': "Number of jets",
        'n_bjets': "Number of b-tagged jets"
    },
    'ee': {
        'pt_1'  : "leading elec p_{T} (GeV)",
        'eta_1' : "leading elec #eta",
        'pt_2'  : "trailing elec p_{T} (GeV)",
        'eta_2' : "trailing elec #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{ee} (GeV)",
	'ipsig_1': "leading elec IP sig",
        'ipsig_2': "trailing elec IP sig",
        'jpt_1': "leading jet p_{T} (GeV)",
        'jpt_2': "trailing jet p_{T} (GeV)",
        'jeta_1': "leading jet #eta",
        'jeta_2': "trailing jet #eta",
        'n_jets': "Number of jets",
        'n_bjets': "Number of b-tagged jets"
    }
}

    
def Plot(histData,histMC,**kwargs):

    era = kwargs.get('era','Run3_2022')
    var = kwargs.get('variable','m_vis')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.5)
    ymax = kwargs.get('ymax',1.5)
    
    # histograms
    h_data = histData.Clone('h_data')
    h_zll = histMC.Clone('h_zll')

    styles.InitData(h_data)

    xtitle = XTitle[chan][var]
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)

    x_data = h_data.GetSumOfWeights()
    x_zll  = h_zll.GetSumOfWeights() 
    
    print('')
    print('Era   : %s'%(era))
    print('Lumi  : %3.0f'%(utils.eraLumi[era]))
    print('')
    print('DY    : %7.0f'%(x_zll))
    print('Data  : %7.0f'%(x_data))
    print('')

    h_tot = h_zll.Clone("total")
    styles.InitTotalHist(h_tot)

    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_zll)
    
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
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.7,0.55,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.055)
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_zll,'DY MC','f')
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
    basedir = '%s/figures'%(utils.outputFolder)
    outputGraphics = '%s/%s_%s_%s_tree.png'%(basedir,var,chan,era)
    canvas.Print(outputGraphics)


if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3',choices=['Run3','Run3_2022','Run3_2023','Run3_2022preEE','Run3_2022postEE','Run3_2023preBPix','Run3_2023postBPix'])
    parser.add_argument('-channel','--channel',dest='channel',default='mm',choices=['ee','mm'])
    parser.add_argument('-var','--var',dest='var',default='m_vis')
    parser.add_argument('-nbins','--nbins',dest='nbins',type=int,default=50)
    parser.add_argument('-xmin','--xmin',dest='xmin',type=float,default=50.)
    parser.add_argument('-xmax','--xmax',dest='xmax',type=float,default=150.)
    parser.add_argument('-hideLegend','--hideLegend',dest='hideLegend',action='store_true')
    args = parser.parse_args()

    baseFolder = '%s/src/CPHiggs/Analysis'%(os.getenv('CMSSW_BASE'))
    
    periods = utils.periods[args.era]
    
    dy_samples = ['DYto2L_M_50_amcatnloFXFX','DYto2L_M_50_0J_amcatnloFXFX','DYto2L_M_50_1J_amcatnloFXFX','DYto2L_M_50_2J_amcatnloFXFX']
    if args.era in ['Run3_2022preEE','Run3_2022postEE']:
        dy_samples.append('DYto2L_M_50_amcatnloFXFX_ext1')
    
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    channel = args.channel
    var = args.var
    era = args.era

    plotLegend = True
    if args.hideLegend:
        plotLegend = False
    
    bins = utils.createBins(nbins,xmin,xmax)

    ROOT.gSystem.cd('')
    histData = ROOT.TH1D('histData','',nbins,array('d',list(bins)))
    histMC = ROOT.TH1D('histMC','',nbins,array('d',list(bins)))
    histsData = {}
    histsMC = {}
    
    cut = 'weight*(pt_1>26.&&fabs(eta_1)<2.4&&pt_2>25.&&fabs(eta_2)<2.4&&iso_1<0.15&&iso_2<0.15&&trg_singlemuon&&os)'
    if args.channel=='ee':
        cut = 'weight*(pt_1>32.&&fabs(eta_1)<2.1&&pt_2>25.&&fabs(eta_2)<2.1&&iso_1<0.15&&iso_2<0.15&&trg_singleelectron&&os)'

    print('')
    dummyCanv = ROOT.TCanvas('dummyCanv','',500,500)
    for period in periods:
        samples = utils.muonSamples[period]
        if args.channel=='ee':
            samples = utils.elecSamples[period]
        for sample in samples:
            name = '%s_%s'%(period,sample)
            filename = '%s/%s/%s/%s/nominal/merged.root'%(utils.tupleFolderV2,period,channel,sample)
            inputfile = ROOT.TFile(filename)
            tree = inputfile.Get('ntuple')
            ROOT.gSystem.cd('')
            hist = ROOT.TH1D(name,'',nbins,array('d',list(bins)))
            dump = '%s>>%s'%(var,name)
            nevents = tree.Draw(dump,cut)
            print('Data : %s -> %3.0f , %3.0f'%(filename,nevents,hist.GetSumOfWeights()))
            histData.Add(histData,hist,1.,1.)

    print('')
    for period in periods:
        yaml_file = baseFolder+'/params/'+period+'.yaml'
        metafile = open(yaml_file,'r')
        metadata = list(yaml.load_all(metafile,Loader=SafeLoader))
        for sample in dy_samples:
            name = '%s_%s'%(period,sample)
            norm = metadata[0]['lumi']*metadata[0][sample]['xs']*metadata[0][sample]['filter_efficiency']/metadata[0][sample]['eff']
            filename = '%s/%s/%s/%s/nominal/merged.root'%(utils.tupleFolderV2,period,channel,sample)
            inputfile = ROOT.TFile(filename)
            tree = inputfile.Get('ntuple')
            ROOT.gSystem.cd('')
            hist = ROOT.TH1D(name,'',nbins,array('d',list(bins)))
            dump = '%s>>%s'%(var,name)
            nevents = tree.Draw(dump,cut)
            print('MC   : %s -> %3.0f , %3.0f'%(filename,nevents,hist.GetSumOfWeights()))
            histMC.Add(histMC,hist,1.,norm)
            

    ROOT.gSystem.cd('')
    
    Plot(histData,histMC,
         variable=var,
         era=era,
         channel=channel,
         plotLegend=plotLegend)
