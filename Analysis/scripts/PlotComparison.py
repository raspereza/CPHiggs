#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# 
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

varTitle = {
    'dPt': 'p_{T}^{reco}/p_{T}^{gen}',
    'dTheta': '#Delta#theta(reco,gen)',
    'dAlpha': '#Delta#alpha(reco,gen)',
    'mass': 'm_{#tau#tau} (GeV)',
    'mass_a1_a1': 'm_{#tau#tau} (GeV)',
}

def Plot(h1,h2,**kwargs):

    var = kwargs.get('var','dPt')
    chan = kwargs.get('chan','mt')

    if 'dPt' in var or 'mass' in var:
        leg1 = 'FastMTT'
    else:
        leg1 = 'SV'
    
    # histograms

    xtitle = varTitle[var]
    ytitle = 'Events'

    styles.InitModel(h1,xtitle,ytitle,2)
    styles.InitModel(h2,xtitle,ytitle,4)

    utils.zeroBinErrors(h1)
    utils.zeroBinErrors(h2)
    
    YMax = h1.GetMaximum()
    if h2.GetMaximum()>YMax: YMax = h2.GetMaximum()
    
    h1.GetYaxis().SetRangeUser(0.,1.2*YMax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h1.Draw('h')
    h2.Draw('hsame')
    
    leg = ROOT.TLegend(0.7,0.5,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetHeader('gg#rightarrowH')
    leg.SetTextSize(0.046)
    leg.AddEntry(h1,leg1,'l')
    leg.AddEntry(h2,'KinFit','l')
    leg.Draw()

    styles.CMS_label(canvas,era='Run3_all',extraText='Simulation')
    
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    outputGraphics = '/eos/home-r/rasp/php-plots/plots/phiCP/'+var+'_'+chan+'.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022postEE', choices=['Run3_2022postEE'])
    parser.add_argument('-channel','--channel', dest='channel', default='tt',choices=['mt','tt'])
    parser.add_argument('-sample', '--sample', dest='sample', default='ggH_sm')
    parser.add_argument('-variable' ,'--variable', dest='variable', default='dAlpha')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=40)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=0.02)
    
    args = parser.parse_args()

    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    era = args.era
    chan = args.channel
    sample = args.sample

    var1 = args.sample+'_'+args.variable+'_sm'
    var2 = args.sample+'_'+args.variable+'_KinFit_sm'
    
    bins = utils.createBins(nbins,xmin,xmax)
    basefolder = '/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection/phiCP'
    filename = '%s/%s_%s_%s.root'%(basefolder,sample,chan,era)
    inputfile = ROOT.TFile(filename,'')
    print(inputfile)
    
    hist1 = inputfile.Get(var1)
    hist2 = inputfile.Get(var2)
    lumi_scale = 10.7
    hist1.Scale(lumi_scale)
    hist2.Scale(lumi_scale)
    
    h1 = utils.rebinHisto(hist1,bins,'_rebinned')
    h2 = utils.rebinHisto(hist2,bins,'_rebinned')
    
    print('')
    print('Nominal : Mean = %6.4f  RMS = %6.4f'%(hist1.GetRMS(),hist1.GetMean()))
    print('KinFit  : Mean = %6.4f  RMS = %6.4f'%(hist2.GetRMS(),hist2.GetMean()))
    print('')
    
    Plot(h1,h2,var=args.variable,chan=chan)
    
