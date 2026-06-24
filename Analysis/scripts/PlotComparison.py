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
    'chi2_kinfit': 'Fit #chi^{2}',
    'dPt': 'p_{T}^{reco}/p_{T}^{gen}',
    'dTheta': '#Delta#theta(reco,gen)',
    'dAlpha': '#Delta#alpha(reco,gen)',
    'mass': 'm_{#tau#tau} (GeV)',
    'mass_a1_a1': 'm_{#tau#tau} (GeV)',
}

def computeSB(h1,h2,h3):

    sb = h1.Clone('sb')
    nbins = sb.GetNbinsX()
    for i in range(1,nbins+1):
        bkgd = 0
        sig = 0
        for j in range(1,i+1):
            bkgd += h2.GetBinContent(j)+h3.GetBinContent(j)
            sig += h1.GetBinContent(j)
        x = sig/math.sqrt(bkgd)
        sb.SetBinContent(i,x)
        sb.SetBinError(i,0)

    sb.GetYaxis().SetRangeUser(1.79,4.)
    # canvas and pads
    canvas = styles.MakeCanvas("canv1","",800,700)
    sb.GetXaxis().SetTitle('#chi^{2} cut')
    sb.GetYaxis().SetTitle('S#sqrt{B}')
    
    sb.Draw('h')
    styles.CMS_label(canvas,era='Run3_all',extraText='')
    
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    outputGraphics = '/eos/home-r/rasp/php-plots/plots/phiCP_datacards/sb_chi2.png'    
    canvas.Print(outputGraphics)


def Plot(h1,h2,h3,**kwargs):

    var = kwargs.get('var','chi2_kinfit')
    chan = kwargs.get('chan','tt')

    if 'dPt' in var or 'mass' in var:
        leg1 = 'FastMTT'
    else:
        leg1 = 'SV'
    
    # histograms

    xtitle = varTitle[var]
    ytitle = 'normalized to unity'

    styles.InitModel(h1,xtitle,ytitle,2)
    styles.InitModel(h2,xtitle,ytitle,4)
    styles.InitModel(h3,xtitle,ytitle,1)

    utils.zeroBinErrors(h1)
    utils.zeroBinErrors(h2)
    utils.zeroBinErrors(h3)
    
    YMax = h1.GetMaximum()
    if h2.GetMaximum()>YMax: YMax = h2.GetMaximum()
    if h3.GetMaximum()>YMax: YMax = h3.GetMaximum()
    
    h1.GetYaxis().SetRangeUser(0.,1.2*YMax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h1.Draw('h')
    h2.Draw('hsame')
    h3.Draw('hsame')
    
    leg = ROOT.TLegend(0.3,0.5,0.6,0.7)
    styles.SetLegendStyle(leg)
    leg.SetHeader('all 3-prongs')
    leg.SetTextSize(0.046)
    leg.AddEntry(h1,'gg#rightarrow H','l')
    leg.AddEntry(h2,'Z#rightarrow#tau#tau','l')
    leg.AddEntry(h3,'j#rightarrow#tau fakes','l')
    leg.Draw()

#    styles.CMS_label(canvas,era='Run3_all',extraText='Simulation')
    
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    outputGraphics = '/eos/home-r/rasp/php-plots/plots/phiCP_datacards/'+var+'_'+chan+'.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022postEE', choices=['Run3','Run3_2022postEE'])
    parser.add_argument('-channel','--channel', dest='channel', default='tt',choices=['mt','tt'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='chi2_kinfit')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=20)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=10.)
    
    args = parser.parse_args()

    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    era = args.era
    chan = args.channel
    variable = args.variable
    
    var1 = 'ggH_sm_'+variable+'_os_iso_all_all_sm'
    var2 = 'ztt_'+variable+'_os_iso_all_all'
    var3 = 'data_'+variable+'_ss_iso_all_all'
    
    bins = utils.createBins(nbins,xmin,xmax)
    
    basefolder = '/eos/home-r/rasp/CPHiggs/Analysis/selection/datacardsPhiCP'
    filename = '%s/%s_%s_x_ipcut1.root'%(basefolder,chan,era)
    inputfile = ROOT.TFile(filename,'READ')
    print(inputfile)
    
    
    hist1 = inputfile.Get(var1)
    hist2 = inputfile.Get(var2)
    hist3 = inputfile.Get(var3)

    lumi_scale = 10.7
    hist1.Scale(lumi_scale)
    hist2.Scale(lumi_scale)
    hist3.Scale(lumi_scale)
    
    norm1 = hist1.GetSumOfWeights()
    norm2 = hist2.GetSumOfWeights()
    norm3 = hist3.GetSumOfWeights()
    
    print('ggH   : %3.0f'%(norm1))
    print('DY    : %3.0f'%(norm2))
    print('Fakes : %3.0f'%(norm3))

    computeSB(hist1,hist2,hist3)
    
    hist1.Scale(1.0/norm1)
    hist2.Scale(1.0/norm2)
    hist3.Scale(1.0/norm3)
    
    h1 = utils.rebinHisto(hist1,bins,'_rebinned')
    h2 = utils.rebinHisto(hist2,bins,'_rebinned')
    h3 = utils.rebinHisto(hist3,bins,'_rebinned')
    
    print('')
    print('ggH   : Mean = %6.4f  RMS = %6.4f'%(hist1.GetRMS(),hist1.GetMean()))
    print('DY    : Mean = %6.4f  RMS = %6.4f'%(hist2.GetRMS(),hist2.GetMean()))
    print('Fakes : Mean = %6.4f  RMS = %6.4f'%(hist3.GetRMS(),hist3.GetMean()))
    print('')
    
    Plot(h1,h2,h3,var=variable,chan=chan)
    
