#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

def header(name,channel):
    title = ''
    if 'aco_lep_pi' in name:
        title = '#mu+#pi'
    if 'aco_lep_rho' in name:
        title = '#mu+#rho'
    if 'aco_lep_a1' in name:
        title = '#mu+a_{1}'

    if channel=='et':
        if 'aco_lep_pi' in name:
            title = 'e+#pi'
        if 'aco_lep_rho' in name:
            title = 'e+#rho'
        if 'aco_lep_a1' in name:
            title = 'e+a_{1}'

    return title


def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    suffix = kwargs.get('suffix','')
    asym = kwargs.get('asym',0.1)
    
    # histograms
    h_even = hists['even'].Clone('h_even')
    h_odd = hists['odd'].Clone('h_odd')

    xtitle = '#phi_{CP} [deg]'
    ytitle = 'normalized'

    styles.InitModel(h_even,xtitle,ytitle,2)
    styles.InitModel(h_odd,xtitle,ytitle,4)

    utils.zeroBinErrors(h_even)
    utils.zeroBinErrors(h_odd)
    
    YMax = h_even.GetMaximum()
    if h_odd.GetMaximum()>YMax: YMax = h_odd.GetMaximum()
    
    h_even.GetYaxis().SetRangeUser(0.,1.4*YMax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h_even.Draw('h')
    h_odd.Draw('hsame')

    leg = ROOT.TLegend(0.65,0.25,0.85,0.5)
    styles.SetLegendStyle(leg)
    leg.SetHeader('%s (A=%5.3f)'%(header(var,chan),asym))
    leg.SetTextSize(0.046)
    leg.AddEntry(h_even,'CP-even','l')
    leg.AddEntry(h_odd,'CP-odd','l')
    leg.Draw()

    styles.CMS_label(canvas,era=era,extraText='Simulation')

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    outputGraphics = os.getenv('CMSSW_BASE')+'/src/CPHiggs/IP/figures/'+var+'_'+chan+'_'+era+suffix+'.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3v2', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix','Run3','Run3v0','Run3v1','Run3v2','Run3v3'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='aco_lep_rho')
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-applyIPSigCut','--applyIPSigCut', dest='applyIPSigCut',action='store_true')
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',action='store_true')
    parser.add_argument('-mtCut','--mtCut',dest='mtCut',action='store_true')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=8)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=360.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    
    args = parser.parse_args()

    era = args.era
    applyMTCut = args.mtCut
    useCrossTrigger = args.useCrossTrigger
    applyIPSigCut = args.applyIPSigCut
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)


    suffix = ''
    suffix_mt = ''
    suffix_xtrig = ''
    suffix_ip = ''
    if useCrossTrigger:
        suffix_xtrig = '_xtrig'
    if applyMTCut:
        suffix_mt = '_mtcut'
    if applyIPSigCut:
        suffix_ip = '_ipcut1'

    suffix = suffix_mt+suffix_xtrig+suffix_ip 
        
    filename = '%s/src/CPHiggs/IP/selection/signal_%s_%s%s.root'%(os.getenv('CMSSW_BASE'),chan,era,suffix)
    
    inputfile = ROOT.TFile(filename,'read')

    histEven = inputfile.Get('even_%s_os_iso_all'%(var))
    histOdd  = inputfile.Get('odd_%s_os_iso_all'%(var))

    histEven.Scale(2.12)
    histOdd.Scale(2.12)
    
#    utils.zeroBinErrors(histEven)
#    utils.zeroBinErrors(histOdd)
    
    normEven = histEven.GetSumOfWeights()
    normOdd  = histOdd.GetSumOfWeights()
    histOdd.Scale(normEven/normOdd)
    print('')
    print('Norm(CP-odd) %5.1f  :  Norm(CP-even) = %5.1f'%(normOdd,normEven))
    print('')
    print('   CP-Odd  |   CP-even')
    print('-----------+-------------')
    asym = 0
    for ib in range(1,nbins+1):
        xOdd = histOdd.GetBinContent(ib)
        xEven = histEven.GetBinContent(ib)
        print('  %6.2f   |    %6.2f'%(xOdd,xEven))
        asym += abs(xEven-xOdd)/(xEven+xOdd)
    print('-----------+-------------')
    asym /= float(nbins)

    print('CP Asymetry = %7.5f'%(asym))
    print('')
    
    hists = {}
    hists['even'] = histOdd
    hists['odd'] = histEven
#    hists['even'] = histEven
#    hists['odd'] = histOdd
    Plot(hists,era=era,var=var,channel=chan,suffix='',asym=asym)
    
    
