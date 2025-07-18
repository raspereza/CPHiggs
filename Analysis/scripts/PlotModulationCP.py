#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

header = {
    'aco_mu_pi': '#mu+#pi',
    'aco_mu_rho': '#mu+#rho',
    'aco_mu_a1': '#mu+a_{1}',
    'aco_mu_a1_FASTMTT_NoMassConstraint': '#mu+a_{1}',
    'aco_mu_a1_FASTMTT_MassConstraint': '#mu+a_{1}',
    'aco_e_pi': 'e+#pi',
    'aco_e_rho': 'e+#rho'    
}


def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    suffix = kwargs.get('suffix','')
    
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
    leg.SetHeader(header[var])
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
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022EE', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='aco_mu_rho')
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',type=int,default=0)
    parser.add_argument('-mtCut','--mtCut',dest='mtCut',type=int,default=0)
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=8)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=360.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    
    args = parser.parse_args()

    era = args.era
    applyMTCut = args.mtCut
    useCrossTrigger = args.useCrossTrigger
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


    folder = '/eos/cms/store/group/phys_tau/ksavva/For_Aliaksei/files/testingzpt/%s/%s/'%(era,chan)
    fileNameOdd  = '%s/GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay/nominal/merged.root'%(folder)
    fileNameEven = '%s/GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay/nominal/merged.root'%(folder)
    
    fileOdd = ROOT.TFile(fileNameOdd,'read')
    fileEven = ROOT.TFile(fileNameEven,'read')

    treeOdd = fileOdd.Get('ntuple')
    treeEven = fileEven.Get('ntuple')
    
    cutDeepTau1 = 'idDeepTau2018v2p5VSjet_1>=4&&idDeepTau2018v2p5VSe_1>=2&&idDeepTau2018v2p5VSmu_1>=1'
    cutDeepTau2 = 'idDeepTau2018v2p5VSjet_2>=4&&idDeepTau2018v2p5VSe_2>=2&&idDeepTau2018v2p5VSmu_2>=1'
    
    cutTau1 = 'pt_1>25&&fabs(eta_1)<2.4&&iso_1<0.15&&os>0.5'
    cutTau2 = 'pt_2>20&&fabs(eta_2)<2.3'

    if chan=='et':
        cutTau1 = 'pt_1>31&&fabs(eta_1)<2.1&&iso_1<0.15&&os>0.5'

    
    weightOdd = '(weight*wt_cp_ps)'
    weightEven = '(weight*wt_cp_sm)'
    weightOddCut = '((fabs%s)<1000.)'%(weightOdd)
    weightEvenCut = '((fabs%s)<1000.)'%(weightEven)
    cutDM = 'decayModePNet_2==1'


    cutIP = 'fabs(ip_LengthSig_1)>1.0'
    if var=='aco_e_pi' or var=='aco_mu_pi':
        cutDM = 'decayModePNet_2==0'
        cutIP = 'fabs(ip_LengthSig_1)>1.0&&fabs(ip_LengthSig_2)>1.0'


    if var=='aco_mu_a1' or var=='aco_mu_a1_FASTMTT_NoMassConstraint' or var=='aco_mu_a1_FASTMTT_MassConstraint':
        cutDM = 'decayModePNet_2==10'

    if var=='aco_e_a1' or var=='aco_e_a1_FASTMTT_NoMassConstraint' or var=='aco_e_a1_FASTMTT_MassConstraint':
        cutDM = 'decayModePNet_2==10'
        
    cutEven = '(%s&&%s&&%s&&%s&&%s)*%s'%(cutTau1,cutTau2,cutDeepTau2,cutDM,weightEvenCut,weightEven)
    cutOdd  = '(%s&&%s&&%s&&%s&&%s)*%s'%(cutTau1,cutTau2,cutDeepTau2,cutDM,weightOddCut,weightOdd)
    cutEvenIP = '(%s&&%s&&%s&&%s&&%s&&%s)*%s'%(cutTau1,cutTau2,cutDeepTau2,cutDM,cutIP,weightEvenCut,weightEven)
    cutOddIP  = '(%s&&%s&&%s&&%s&&%s&&%s)*%s'%(cutTau1,cutTau2,cutDeepTau2,cutDM,cutIP,weightOddCut,weightOdd)

    histEven = ROOT.TH1D('histEven','',nbins,array('d',list(bins)))
    histOdd  = ROOT.TH1D('histOdd', '',nbins,array('d',list(bins)))
    histEvenIP = ROOT.TH1D('histEvenIP','',nbins,array('d',list(bins)))
    histOddIP  = ROOT.TH1D('histOddIP', '',nbins,array('d',list(bins)))

    varToPlot = '(180.*(%s/TMath::Pi()))'%(var)


    dummy = ROOT.TCanvas('dummy','dummy',500,500)
    
    treeOdd.Draw(varToPlot+'>>histOdd',cutOdd)
    treeEven.Draw(varToPlot+'>>histEven',cutEven)

    treeOdd.Draw(varToPlot+'>>histOddIP',cutOddIP)
    treeEven.Draw(varToPlot+'>>histEvenIP',cutEvenIP)


    normEven = histEven.GetSumOfWeights()
    normOdd  = histOdd.GetSumOfWeights()
    histEven.Scale(1.0/normEven)
    histOdd.Scale(1.0/normOdd)
    histEvenIP.Scale(1.0/normEven)
    histOddIP.Scale(1.0/normOdd)
    print('')
    print('Norm(CP-odd) %4.2f  :  Norm(CP-even) = %4.2f'%(normOdd,normEven))
    print('')
    print('  CP-Odd  |   CP-even')
    print('----------+------------')
    asym = 0
    for ib in range(1,nbins+1):
        xOdd = histOdd.GetBinContent(ib)
        xEven = histEven.GetBinContent(ib)
        print('  %5.3f   |    %5.3f'%(xOdd,xEven))
        asym += abs(xEven-xOdd)/(xEven+xOdd)
    print('----------+------------')
    asymIP = 0
    for ib in range(1,nbins+1):
        xOdd = histOddIP.GetBinContent(ib)
        xEven = histEvenIP.GetBinContent(ib)
        print('  %5.3f   |    %5.3f'%(xOdd,xEven))
        asymIP += abs(xEven-xOdd)/(xEven+xOdd)
    print('')

    asym /= float(nbins)
    asymIP /= float(nbins)

    print('CP Asymetry : no IP cut = %5.3f  ---  IP cut = %5.3f'%(asym,asymIP))
    print('')
    
    hists = {}
    hists['even'] = histEven
    hists['odd'] = histOdd
    Plot(hists,era=era,var=var,channel=chan,suffix='')
    
    histsIP = {}
    histsIP['even'] = histEvenIP
    histsIP['odd'] = histOddIP
    Plot(histsIP,era=era,var=var,channel=chan,suffix='_ip')
    
