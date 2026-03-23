#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def header(name,channel):
    title = ''
    if 'aco_lep_pi' in name:
        title = '#mu+#pi'
    if 'aco_lep_rho' in name:
        title = '#mu+#rho'
    if 'aco_lep_a1' in name:
        title = '#mu+a_{1}'

    if channel=='tt':
        if 'aco_pi_a1' in name:
            title = '#pi+a_{1}'
        if 'aco_rho_a1' in name:
            title = '#rho+a_{1}'
        if 'aco_a1_a1' in name:
            title = 'a_{1}+a_{1}'
        
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
    var = kwargs.get('var','theta_GJ')
    chan = kwargs.get('channel','mt')
    asym = kwargs.get('asym',0.1)
    
    # histograms
    h_even = hists['even'].Clone('h_even')
    h_odd = hists['odd'].Clone('h_odd')

    xtitle = '#phi_{CP} [deg]'
    ytitle = 'Events'

    styles.InitModel(h_even,xtitle,ytitle,2)
    styles.InitModel(h_odd,xtitle,ytitle,4)

    utils.zeroBinErrors(h_even)
    utils.zeroBinErrors(h_odd)
    
    YMax = h_even.GetMaximum()
    if h_odd.GetMaximum()>YMax: YMax = h_odd.GetMaximum()
    
    h_even.GetYaxis().SetRangeUser(0.,1.2*YMax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h_even.Draw('h')
    h_odd.Draw('hsame')

    leg = ROOT.TLegend(0.65,0.17,0.85,0.37)
    styles.SetLegendStyle(leg)
    leg.SetHeader('%s (A=%5.3f)'%(header(var,chan),asym))
    leg.SetTextSize(0.047)
    leg.AddEntry(h_even,'CP-even','l')
    leg.AddEntry(h_odd,'CP-odd','l')
    leg.Draw()

    styles.CMS_label(canvas,era='Run3_all',extraText='Simulation')

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    graphicFolder = '/eos/home-r/rasp/php-plots/plots/phiCP'
    outputGraphics = graphicFolder+'/'+var+'_'+chan+'_'+era+'.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022postEE', choices=['Run3_2022postEE','Run3_2022preEE'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='aco_rho_a1')
    parser.add_argument('-channel','--channel', dest='channel', default='tt',choices=['mt','tt','et'])
    parser.add_argument('-sample','--sample',dest='sample',default='ggH_sm')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=8)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=360.)
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    sample = args.sample
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    basefolder = '/afs/cern.ch/work/r/rasp/CPHiggs/Analysis/selection/phiCP'
    filename = '%s/%s_%s_%s.root'%(basefolder,sample,chan,era)
    if os.path.isfile(filename):
        print('')
        print('Opening file %s'%(filename))
        print('')
    else:
        print('')
        print('File %s does not exist'%(filename))
        print('')
        exit()
    
    inputfile = ROOT.TFile(filename,'read')

    bins = utils.createBins(nbins,xmin,xmax)
    
    name_sm = '%s_%s_sm'%(sample,var)
    name_ps = '%s_%s_ps'%(sample,var)
    histEven = utils.rebinHisto(inputfile.Get(name_sm),bins,'_rebinned')
    histOdd  = utils.rebinHisto(inputfile.Get(name_ps),bins,'_rebinned')

    lumi_scale = 10.7
    histEven.Scale(lumi_scale)
    histOdd.Scale(lumi_scale)
    
    normEven = histEven.GetSumOfWeights()
    normOdd  = histOdd.GetSumOfWeights()
    
    if 'aco' in var:
        utils.symmetrize(histEven)
        utils.symmetrize(histOdd)
    
    asym = 0
    if 'aco' in var:
        print('')
        print('Norm(CP-odd) %5.1f  :  Norm(CP-even) = %5.1f'%(normOdd,normEven))
        print('')
        print('   CP-Odd  |   CP-even')
        print('-----------+-------------')
        LogL_Even = 0.
        LogL_Odd = 0.
        for ib in range(1,nbins+1):
            xOdd = histOdd.GetBinContent(ib)
            xEven = histEven.GetBinContent(ib)
            print('  %6.2f   |    %6.2f'%(xOdd,xEven))
            asym += abs(xEven-xOdd)/(xEven+xOdd)
            term = (xEven-xOdd)*math.log(xEven/xOdd)
            LogL_Even += term
            term = (xOdd-xEven)*math.log(xOdd/xEven)
            LogL_Odd += term
        print('-----------+-------------')
        asym /= float(nbins)
        print('CP Asymetry   = %7.5f'%(asym))
        print('LogL(CP-even) = %7.5f'%(LogL_Even))
        print('LogL(CP-odd)  = %7.5f'%(LogL_Odd))
        print('')
    
    hists = {}
    hists['even'] = histOdd
    hists['odd'] = histEven
    Plot(hists,era=era,var=var,channel=chan,asym=asym)
    
    
