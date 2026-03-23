#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# 
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def Plot(h1,**kwargs):

    var = kwargs.get('var','chi2_pi_a1')
    title = kwargs.get('title','#chi2^2')
    
    # histograms

    xtitle = title
    ytitle = 'Events'

    styles.InitModel(h1,xtitle,ytitle,2)

    utils.zeroBinErrors(h1)
    
    YMax = h1.GetMaximum()
    
    h1.GetYaxis().SetRangeUser(0.,1.2*YMax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h1.Draw('h')
    
    styles.CMS_label(canvas,era='Run3_all',extraText='Simulation')

    text = ROOT.TLatex()
    text.SetTextSize(0.07)
    text.DrawLatexNDC(0.3,0.8,'a_{1}+a_{1}')
    
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
    parser.add_argument('-channel','--channel', dest='channel', default='tt',choices=['mt','et'])
    parser.add_argument('-sample', '--sample', dest='sample', default='ggH_sm')
    parser.add_argument('-variable' ,'--variable', dest='variable', default='chi2_a1_a1')
    parser.add_argument('-title','--title',dest='title',default='#chi^{2}')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=50)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=5.)
    
    args = parser.parse_args()

    var  = args.variable
    title = args.title
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    era = args.era
    chan = args.channel
    sample = args.sample

    var = args.sample+'_'+args.variable+'_sm'
    
    bins = utils.createBins(nbins,xmin,xmax)
    basefolder = '/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection/phiCP'
    filename = '%s/%s_%s_%s.root'%(basefolder,sample,chan,era)
    inputfile = ROOT.TFile(filename,'')
    print(inputfile)

    lumi_scale = 10.7
    hist = inputfile.Get(var)
    h = utils.rebinHisto(hist,bins,'_rebinned')
    h.Scale(lumi_scale)

    print('Mean = %6.4f   RMS = %6.4f'%(h.GetMean(),h.GetRMS()))
    
    Plot(h,title=title,chan=chan,var=var)
    
