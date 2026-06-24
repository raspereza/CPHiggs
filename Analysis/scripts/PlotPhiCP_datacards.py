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
        if 'phicp_pi_a1_3pr' in name:
            title = '#pi+a_{1}'
        if 'phicp_rho_a1_3pr' in name:
            title = '#rho+a_{1}'
        if 'phicp_a1_3pr_a1_3pr' in name:
            title = 'a_{1}+a_{1}'
        
    if channel=='et':
        if 'aco_lep_pi' in name:
            title = 'e+#pi'
        if 'aco_lep_rho' in name:
            title = 'e+#rho'
        if 'aco_lep_a1' in name:
            title = 'e+a_{1}'

    return title


def computeAcceptance(hist_sig,hist_dy,hist_fakes):

    nbins = hist_sig.GetNbinsX()
    total_sig = 0
    total_dy = 0
    total_fakes = 0
    for ib in range(1,nbins+1):
        total_sig += hist_sig.GetBinContent(ib)
        total_dy += hist_dy.GetBinContent(ib)
        total_fakes += hist_fakes.GetBinContent(ib)
    lastbin_sig   = hist_sig.GetBinContent(nbins)
    lastbin_dy    = hist_dy.GetBinContent(nbins)
    lastbin_fakes = hist_fakes.GetBinContent(nbins)
    f_sig = (total_sig-lastbin_sig)/total_sig
    f_dy  = (total_dy-lastbin_dy)/total_dy
    f_fakes  = (total_fakes-lastbin_fakes)/total_fakes
    return f_sig,f_dy,f_fakes

def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','theta_GJ')
    chan = kwargs.get('channel','tt')
    asym = kwargs.get('asym',0.1)
    q = kwargs.get('q',1.0)
    yupper = kwargs.get('yupper',5.)
    ylower = kwargs.get('ylower',1.)
    cut = kwargs.get('cut',False)
    
    # histograms
    h_even = hists['even'].Clone('h_even')
    h_odd = hists['odd'].Clone('h_odd')
    h_dy = hists['dy'].Clone('h_dy')
    h_fakes = hists['fakes'].Clone('h_fakes')
    
    xtitle = '#phi_{CP} [deg]'
    ytitle = 'Events'

    styles.InitModel(h_even,xtitle,ytitle,2)
    styles.InitModel(h_odd,xtitle,ytitle,4)

    styles.InitHist(h_dy,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_fakes,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)

    h_dy.Add(h_dy,h_fakes,1.,1.)
    
    utils.zeroBinErrors(h_even)
    utils.zeroBinErrors(h_odd)
    utils.zeroBinErrors(h_dy)
    utils.zeroBinErrors(h_fakes)
    
    h_dy.GetYaxis().SetRangeUser(ylower,yupper)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",800,700)
    
    h_dy.Draw('e2')
    h_fakes.Draw('e2same')
    h_even.Draw('hsame')
    h_odd.Draw('hsame')

    leg = ROOT.TLegend(0.25,0.65,0.45,0.88)
    styles.SetLegendStyle(leg)
    leg.SetHeader('%s (Q=%4.2f)'%(header(var,chan),q))
    leg.SetTextSize(0.04)
    leg.AddEntry(h_even,'CP-even','l')
    leg.AddEntry(h_odd,'CP-odd','l')
    leg.AddEntry(h_dy,'Z#rightarrow#tau#tau','f')
    leg.AddEntry(h_fakes,'Fakes','f')
    leg.Draw()

    styles.CMS_label(canvas,era='Run3_all',extraText='')

    canvas.SetLogy(True)
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    graphicFolder = '/eos/home-r/rasp/php-plots/plots/phiCP_datacards'
    outputGraphics = graphicFolder+'/'+var+'_'+chan+'_'+era+'.png'
    if cut:
        outputGraphics = graphicFolder+'/'+var+'_'+chan+'_'+era+'_cut.png'
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022postEE', choices=['Run3_2022postEE','Run3_2022preEE'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='phicp_rho_a1_3pr')
    parser.add_argument('-channel','--channel', dest='channel', default='tt',choices=['mt','tt','et'])
    parser.add_argument('-sample','--sample',dest='sample',default='ggH_sm')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=8)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=360.)
    parser.add_argument('-ylower','--ylower',dest='ylower', type=float, default=10.)
    parser.add_argument('-yupper','--yupper',dest='yupper', type=float, default=5000.)
    parser.add_argument('-cut','--cut',dest='cut',action='store_true')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    sample = args.sample
    ylower = args.ylower
    yupper = args.yupper
    cut = args.cut
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    basefolder = '/eos/user/r/rasp/CPHiggs/Analysis/selection/datacardsPhiCP'
    filename = '%s/%s_%s_x_ipcut1.root'%(basefolder,chan,era)
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
    
    name_sm = 'ggH_sm_%s_os_iso_all_all_sm'%(var)
    name_ps = 'ggH_sm_%s_os_iso_all_all_ps'%(var)
    name_dy = 'ztt_%s_os_iso_all_all'%(var)
    name_fakes = 'data_%s_ss_iso_all_all'%(var)    
    histEven = utils.rebinHisto(inputfile.Get(name_sm),bins,'_rebinned')
    histOdd  = utils.rebinHisto(inputfile.Get(name_ps),bins,'_rebinned')
    histDY = utils.rebinHisto(inputfile.Get(name_dy),bins,'_rebinned')
    histFakes  = utils.rebinHisto(inputfile.Get(name_fakes),bins,'_rebinned')
    
    f_sig,f_dy,f_fakes = 0.88,0.51,0.34
    print('')
    print('Acceptance : Signal=%4.2f  DY=%4.2f  Fakes=%4.2f'%(f_sig,f_dy,f_fakes))
    print('')
    
    scale_sig = 10.7
    scale_dy = 10.7
    scale_fakes = 10.7
    if cut:
        scale_sig *= f_sig
        scale_dy  *= f_dy
        scale_fakes *= f_fakes
        
    histEven.Scale(scale_sig)
    histOdd.Scale(scale_sig)
    histDY.Scale(scale_dy)
    histFakes.Scale(scale_fakes)
    
    normEven = histEven.GetSumOfWeights()
    normOdd  = histOdd.GetSumOfWeights()
    normFakes = histFakes.GetSumOfWeights()
    normDY = histDY.GetSumOfWeights()
    
    if 'phicp' in var:
        utils.symmetrize(histEven)
        utils.symmetrize(histOdd)
        nbins = histDY.GetNbinsX()
        xDY = normDY/float(nbins)
        xFakes = normFakes/float(nbins)
        for ib in range(1,nbins+1):
            histDY.SetBinContent(ib,xDY)
            histDY.SetBinError(ib,0.)
            histFakes.SetBinContent(ib,xFakes)
            histFakes.SetBinError(ib,0.)
    
    asym = 0.0
    q = 1.0
    if 'phicp' in var:
        print('')
        print('Norm(CP-odd)=%3.1f Norm(CP-even)=%3.1f Norm(DY)=%3.1f Norm(Fakes)=%3.1f'%(normOdd,normEven,normDY,normFakes))
        print('')
        print('   CP-Odd  |   CP-even')
        print('-----------+-------------')
        for ib in range(1,nbins+1):
            xOdd  = histOdd.GetBinContent(ib)
            xEven = histEven.GetBinContent(ib)
            xBkg  = histDY.GetBinContent(ib) + histFakes.GetBinContent(ib)
            nOdd  = xOdd + xBkg 
            nEven = xEven + xBkg
            print('  %6.2f   |    %6.2f'%(xOdd,xEven))
            asym += abs(xEven-xOdd)/(xEven+xOdd)
            term = math.pow(nEven/nOdd,nEven)
            q *= term
        print('-----------+-------------')
        asym /= float(nbins)
        print('CP Asymetry   = %5.3f'%(asym))
        print('Q             = %5.3f'%(q))
        print('')
    
    hists = {}
    hists['even'] = histOdd
    hists['odd'] = histEven
    hists['dy'] = histDY
    hists['fakes'] = histFakes
    Plot(hists,era=era,var=var,channel=chan,ylower=ylower,yupper=yupper,asym=asym,q=q,cut=cut)
    
    
