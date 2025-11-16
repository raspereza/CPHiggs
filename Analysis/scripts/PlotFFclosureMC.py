#! /usr/bin/env python3
# Author: Alexei Raspereza (November 2025)
# Closure of MC 
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def ExtractHistosFF(f,sample,var,region,ff,bins):

    hists = {}

    for typ in utils.type_labels:
        nameInput='%s_%s_%s_%s'%(sample,var,region,typ)
        name='%s_%s'%(sample,typ)
        hists[name]=utils.rebinHisto(f.Get(nameInput),bins,'rebinned')
        nameInput='%s_%s_%s_%s_%s'%(sample,var,region,ff,typ)
        print(nameInput)
        name='%s_%s_%s'%(sample,ff,typ)
        hists[name] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')

    return hists
    
def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    sample = kwargs.get('sample','wjets')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.501)
    ymax = kwargs.get('ymax',1.499)
    ff = kwargs.get('ff','wj')
    ff_version = kwargs.get('ff_version','v3')

    h_data = hists['%s_had'%(sample)]
    h_model = hists['%s_%s_had'%(sample,ff)]

    x_data = h_data.GetSumOfWeights()
    x_model = h_model.GetSumOfWeights()

    styles.InitData(x_data)
    xtitle = utils.XTitle[chan][var]
    styles.InitHist(h_model,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)
    
    print('')
    print('Yields ->')
    print('Model   : %7.0f'%(x_model))
    print('Data    : %7.0f'%(x_data))
    print('')
                
    h_tot = h_model.Clone("total")
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoRatio(h_data,h_model,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_model,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_model)

    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    h_data.GetYaxis().SetRangeUser(0.,1.2*YMax)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("events")
#    h_ratio.GetXaxis().SetTitleSize(0.05)
#    h_ratio.GetYaxis().SetTitleSize(0.05)
    h_ratio.GetYaxis().SetTitle("obs/exp")
    h_ratio.GetXaxis().SetTitle(xtitle)


    #######################
    #### Closure tests ####
    #######################
    
    # canvas and pads
    canvas = styles.MakeCanvas("canv_den","",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    
    
    h_data.Draw('e1')
    h_model.Draw('e2same')
    h_tot.Draw('e2same')
    h_data.Draw('e1same')

    leg = ROOT.TLegend(0.72,0.53,0.85,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.SetHeader(sample+' closure')
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_model,'model','f')
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
    outputFolder = '/eos/home-r/rasp/php-plots/plots/FFclosure_%s/%s/mc_%s'%(chan,ff_version,sample)
    outputGraphics = '%s/%s.png'%(outputFolder,var)
    canvas.Print(outputGraphics)

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
    parser.add_argument('-ymin','--ymin',dest='ymin',type=float,default=0.201)
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.899)
    parser.add_argument('-region','--region',dest='region',default='lowmt_os_iso')
    parser.add_argument('-ff_version','--ff_version',dest='ff_version',default='v4')
    parser.add_argument('-sample','--sample',dest='sample',default='wjets')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    sample = args.sample
    region = args.region
    ff_version = args.ff_version
    plotLegend = True

    if var in ['eta_1','eta_2','CMetQCD','CMetW','bdt_ditau','bdt_fakes','bdt_signal']:
        plotLegend = False

    ff = 'mc_wj'
    if sample=='top':
        ff = 'mc_top'
    
    basedir = utils.outputFolder+'/selection/baseline'

    bins = utils.createBins(nbins,xmin,xmax)

    inputFileName = '%s/%s_%s_x_ipcut1_ff_%s.root'%(basedir,chan,era,ff_version)
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
    hists = ExtractHistosFF(inputFile,sample,var,region,ff,bins)
    Plot(hists,era=era,var=var,channel=chan,
         sample=sample,ymin=ymin,ymax=ymax,
         plotLegend=plotLegend,ff=ff,ff_version=ff_version)

    

