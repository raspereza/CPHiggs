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
    if 'theta_GJ' in name:
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
        if 'theta_GJ' in name:
            title = 'e+a_{1}'

    return title


def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','theta_GJ')
    chan = kwargs.get('channel','mt')
    asym = kwargs.get('asym',0.1)
    asym2 = kwargs.get('asym2',0.1)
    ymin = kwargs.get('ymin',0.)
    ymax = kwargs.get('ymax',60.)
    
    # histograms
    h_even = hists['sm'].Clone('h_sm')
    h_odd = hists['ps'].Clone('h_ps')
    h2_even = hists['sm_2'].Clone('h_sm_2')
    h2_odd = hists['ps_2'].Clone('h_ps_2')    
    
    xtitle = '#phi_{CP} [deg]'
    ytitle = 'Events'
    
    styles.InitModel(h_even,xtitle,ytitle,2)
    styles.InitModel(h_odd,xtitle,ytitle,4)
    styles.InitModel(h2_even,xtitle,ytitle,2)
    styles.InitModel(h2_odd,xtitle,ytitle,4)

    h_even.SetLineStyle(2)
    h_odd.SetLineStyle(2)
    h2_even.SetLineStyle(1)
    h2_odd.SetLineStyle(1)
    
    h_even.SetLineWidth(3)
    h_odd.SetLineWidth(3)
    h2_even.SetLineWidth(3)
    h2_odd.SetLineWidth(3)
    
    #    utils.zeroBinErrors(h_even)
    #    utils.zeroBinErrors(h_odd)
    utils.zeroBinErrors(h2_even)
    utils.zeroBinErrors(h2_odd)
    
    YMax = h_even.GetMaximum()
    h_even.GetYaxis().SetRangeUser(ymin,ymax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h_even.Draw('h')
    h_odd.Draw('hsame')
    h2_even.Draw('hsame')
    h2_odd.Draw('hsame')

    leg = ROOT.TLegend(0.25,0.17,0.45,0.33)
    styles.SetLegendStyle(leg)
    leg.SetHeader('Nominal A=%5.3f'%(asym))
    leg.SetTextSize(0.035)
    leg.AddEntry(h_even,'CP-even','l')
    leg.AddEntry(h_odd,'CP-odd','l')
    leg.Draw()

    leg1 = ROOT.TLegend(0.65,0.17,0.85,0.33)
    styles.SetLegendStyle(leg1)
    leg1.SetHeader('KinFit A=%5.3f'%(asym2))
    leg1.SetTextSize(0.035)
    leg1.AddEntry(h2_even,'CP-even','l')
    leg1.AddEntry(h2_odd,'CP-odd','l')
    leg1.Draw()

    styles.CMS_label(canvas,era='Run3_all',extraText='Simulation')

    text = ROOT.TLatex()
    text.SetTextSize(0.07)
    text.DrawLatexNDC(0.25,0.86,header(var,chan))
    
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    graphicFolder = '/eos/home-r/rasp/php-plots/plots/phiCP'
    outputGraphics = graphicFolder+'/'+var+'_'+chan+'_'+era+'_comp.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022postEE', choices=['Run3_2022postEE','Run3_2022preEE'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='aco_rho_a1')
    parser.add_argument('-suffix','--suffix',dest='suffix',default='KinFit')
    parser.add_argument('-channel','--channel', dest='channel', default='tt',choices=['mt','tt'])
    parser.add_argument('-sample','--sample',dest='sample',default='ggH_sm')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=8)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=360.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=60)
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    sample = args.sample
    suffix = args.suffix
    ymin = args.ymin
    ymax = args.ymax
    
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
    name_2_sm = '%s_%s_%s_sm'%(sample,var,suffix)
    name_2_ps = '%s_%s_%s_ps'%(sample,var,suffix)
    hist_ps = utils.rebinHisto(inputfile.Get(name_sm),bins,'_rebinned')
    hist_sm  = utils.rebinHisto(inputfile.Get(name_ps),bins,'_rebinned')
    hist_2_ps = utils.rebinHisto(inputfile.Get(name_2_sm),bins,'_rebinned')
    hist_2_sm  = utils.rebinHisto(inputfile.Get(name_2_ps),bins,'_rebinned')

    lumi_scale = 10.7
    hist_sm.Scale(lumi_scale)
    hist_ps.Scale(lumi_scale)
    hist_2_sm.Scale(lumi_scale)
    hist_2_ps.Scale(lumi_scale)
    
    
    norm_sm = hist_sm.GetSumOfWeights()
    norm_ps  = hist_ps.GetSumOfWeights()
    norm_2_sm = hist_2_sm.GetSumOfWeights()
    norm_2_ps  = hist_2_ps.GetSumOfWeights()
    
    utils.symmetrize(hist_sm)
    utils.symmetrize(hist_ps)
    utils.symmetrize(hist_2_sm)
    utils.symmetrize(hist_2_ps)

    
    asym = 0
    logL = 0.
    print('')
    print('Default option ->')
    print('Norm(CP-even) %5.1f  :  Norm(CP-odd) = %5.1f'%(norm_2_ps,norm_2_sm))
    print('')
    print('   CP-Odd  |   CP-even')
    print('-----------+-------------')
    for ib in range(1,nbins+1):
        xOdd = hist_sm.GetBinContent(ib)
        xEven = hist_ps.GetBinContent(ib)
        print('  %6.2f   |    %6.2f'%(xOdd,xEven))
        asym += abs(xEven-xOdd)/(xEven+xOdd)
        term = xOdd - xEven + xEven*math.log(xEven/xOdd)
        logL += term
    print('-----------+-------------')
    asym /= float(nbins)
    print('CP Asymetry = %5.3f'%(asym))
    print('LogL = %5.2f'%(logL))
    print('')
    
    asym2 = 0
    logL2 = 0.
    print('')
    print('Kinematic fit ->')
    print('Norm(CP-even) %5.1f  :  Norm(CP-odd) = %5.1f'%(norm_ps,norm_sm))
    print('')
    print('   CP-Odd  |   CP-even')
    print('-----------+-------------')
    for ib in range(1,nbins+1):
        xOdd = hist_2_sm.GetBinContent(ib)
        xEven = hist_2_ps.GetBinContent(ib)
        print('  %6.2f   |    %6.2f'%(xOdd,xEven))
        asym2 += abs(xEven-xOdd)/(xEven+xOdd)
        term = xOdd - xEven + xEven*math.log(xEven/xOdd)
        logL2 += term
    print('-----------+-------------')
    asym2 /= float(nbins)
    print('CP Asymetry = %5.3f'%(asym2))
    print('LogL = %5.2f'%(logL2))
    print('')
    
    hists = {}
    hists['sm'] = hist_sm
    hists['ps'] = hist_ps
    hists['sm_2'] = hist_2_sm
    hists['ps_2'] = hist_2_ps
    if 'aco_a1_a1' in var:
        for ib in [1,8]:
            hists['sm_2'].SetBinContent(ib,0.98*hists['sm_2'].GetBinContent(ib))
            hists['ps_2'].SetBinContent(ib,0.98*hists['ps_2'].GetBinContent(ib))
        for ib in [4,5]:
            hists['sm_2'].SetBinContent(ib,1.04*hists['sm_2'].GetBinContent(ib))
            hists['ps_2'].SetBinContent(ib,1.04*hists['ps_2'].GetBinContent(ib))

    Plot(hists,era=era,var=var,channel=chan,asym=asym,asym2=asym2,ymin=ymin,ymax=ymax)
    
    
