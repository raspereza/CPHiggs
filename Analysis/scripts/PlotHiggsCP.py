#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting phi(CP) variables
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

import yaml
from yaml.loader import SafeLoader

header = {
    'aco_mu_pi': '#mu+#pi',
    'aco_mu_rho': '#mu+#rho',
    'aco_mu_a1_FASTMTT_MassConstraint': '#mu+a_{1}',
    'aco_e_pi': 'e+#pi',
    'aco_e_rho': 'e+#rho',
    'aco_e_a1_FASTMTT_MassConstraint' : 'e+a_{1}'    
}

varToChan = {
    'aco_mu_pi' : 'mt',
    'aco_mu_rho' : 'mt',
    'aco_mu_a1_FASTMTT_MassConstraint' : 'mt',
    'aco_e_pi' : 'et',
    'aco_e_rho' : 'et',
    'aco_e_a1_FASTMTT_MassConstraint' : 'et',
}


def Plot(hists,**kwargs):

    var = kwargs.get('var','m_vis')
    suffix = kwargs.get('suffix','')
    asym = kwargs.get('asym',1.0)
    
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
    
    h_even.GetYaxis().SetRangeUser(0.,1.4*YMax)

    # canvas and pads
    canvas = styles.MakeCanvas("canv%s"%(suffix),"",800,700)
    
    h_even.Draw('h')
    h_odd.Draw('hsame')
    legTitle = '%s  A = %5.3f'%(header[var],asym)
    
    leg = ROOT.TLegend(0.45,0.25,0.75,0.5)
    styles.SetLegendStyle(leg)
    leg.SetHeader(legTitle)
    leg.SetTextSize(0.046)
    leg.AddEntry(h_even,'CP-even','l')
    leg.AddEntry(h_odd,'CP-odd','l')
    leg.Draw()

    txt = 'no IPsig cut'
    if suffix=='_ip':
        txt = 'IPSig(lep)>1.0'
    text = ROOT.TText(0.3,0.8,txt)
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.05)
    text.Draw()
    styles.CMS_label(canvas,era='Run3_simulation',extraText='Simulation')
    
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    outputGraphics = '/eos/home-r/rasp/php-plots/plots/phiCP/'+var+suffix+'.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-variable' ,'--variable', dest='variable', default='aco_mu_rho',choices=['aco_mu_pi','aco_mu_rho','aco_mu_a1_FASTMTT_MassConstraint','aco_e_pi','aco_e_rho','aco_e_a1_FASTMTT_MassConstraint'])
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=8)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=360.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    parser.add_argument('-noIpCutTau','--noIpCutTau',dest='noIpCutTau',action='store_true')
    
    args = parser.parse_args()

    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    noIpCutTau = args.noIpCutTau

    sample = 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'
    
    chan = varToChan[var]
    bins = utils.createBins(nbins,xmin,xmax)
    
    cutDeepTau = 'idDeepTau2018v2p5VSjet_2>=7&&idDeepTau2018v2p5VSe_2>=2&&idDeepTau2018v2p5VSmu_2>=4'
    if chan=='et':
        cutDeepTau = 'idDeepTau2018v2p5VSjet_2>=7&&idDeepTau2018v2p5VSe_2>=6&&idDeepTau2018v2p5VSmu_2>=4'
        
    
    cutLep = 'pt_1>26&&fabs(eta_1)<2.4&&iso_1<0.15&&os>0.5&&trg_singlemuon>0.5&&mt_1<65.'
    cutTau = 'pt_2>20&&fabs(eta_2)<2.5&&os>0.5'

    if chan=='et':
        cutLep = 'pt_2>32&&fabs(eta_1)<2.1&&iso_1<0.15&&os>0.5&&trg_singleelectron>0.5&&mt_1<65.'

    
    weightOdd = '(weight*wt_cp_ps)'
    weightEven = '(weight*wt_cp_sm)'
    weightOddCut = '((fabs%s)<1000.)'%(weightOdd)
    weightEvenCut = '((fabs%s)<1000.)'%(weightEven)

    cutIP = 'fabs(ip_LengthSig_1)>1.0'
    cutDM = 'decayModePNet_2==1&&decayMode_2==1&&pion_E_split_2>0.2'

    if var=='aco_e_pi' or var=='aco_mu_pi':
        if noIpCutTau:
            cutDM = 'decayModePNet_2==0'
        else:
            cutDM = 'decayModePNet_2==0&&fabs(ip_LengthSig_2)>1.25'
            
    if var=='aco_e_a1_FASTMTT_MassConstraint' or var=='aco_mu_a1_FASTMTT_MassConstraint':
        cutDM = 'decayModePNet_2==10&&hasRefitSV_2'

        
    cutEven = '(%s&&%s&&%s&&%s&&%s)*%s'%(cutLep,cutTau,cutDeepTau,cutDM,weightEvenCut,weightEven)
    cutOdd  = '(%s&&%s&&%s&&%s&&%s)*%s'%(cutLep,cutTau,cutDeepTau,cutDM,weightOddCut,weightOdd)
    cutEvenIP = '(%s&&%s&&%s&&%s&&%s&&%s)*%s'%(cutLep,cutTau,cutDeepTau,cutDM,cutIP,weightEvenCut,weightEven)
    cutOddIP  = '(%s&&%s&&%s&&%s&&%s&&%s)*%s'%(cutLep,cutTau,cutDeepTau,cutDM,cutIP,weightOddCut,weightOdd)
    dummy = ROOT.TCanvas('dummy','dummy',500,500)
    
    varToPlot = '(180.*(%s/TMath::Pi()))'%(var)

    baseFolder = '%s/src/CPHiggs/Analysis/'%(os.getenv('CMSSW_BASE'))
    eras = ['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix']

    histEvenEra = {}
    histOddEra = {}
    histEvenIPEra = {}
    histOddIPEra = {}
    
    for era in eras:

        print('processing %s'%(era)) 
        yaml_file = baseFolder+'/params/'+era+'.yaml'
        metafile = open(yaml_file,'r')
        metadata = list(yaml.load_all(metafile,Loader=SafeLoader))
        norm = metadata[0]['lumi']*metadata[0][sample]['filter_efficiency']*metadata[0][sample]['xs']/metadata[0][sample]['eff']

        ROOT.gROOT.cd('')
        histEvenEra[era] = ROOT.TH1D('histEven_%s'%(era),'',nbins,array('d',list(bins)))
        histOddEra[era]  = ROOT.TH1D('histOdd_%s'%(era), '',nbins,array('d',list(bins)))
        histEvenIPEra[era] = ROOT.TH1D('histEvenIP_%s'%(era),'',nbins,array('d',list(bins)))
        histOddIPEra[era]  = ROOT.TH1D('histOddIP_%s'%(era), '',nbins,array('d',list(bins)))

        folder = '%s/%s/%s/%s'%(utils.tupleFolderV2,era,chan,sample)
        fileName  = '%s/nominal/merged.root'%(folder)
        inputFile = ROOT.TFile(fileName,'read')
        tree = inputFile.Get('ntuple')

        ROOT.gROOT.cd('')
        tree.Draw(varToPlot+'>>histOdd_%s'%(era),cutOdd)
        tree.Draw(varToPlot+'>>histEven_%s'%(era),cutEven)
        tree.Draw(varToPlot+'>>histOddIP_%s'%(era),cutOddIP)
        tree.Draw(varToPlot+'>>histEvenIP_%s'%(era),cutEvenIP)

        ROOT.gROOT.cd('')        
        histOddEra[era].Scale(norm)
        histEvenIPEra[era].Scale(norm)
        histOddIPEra[era].Scale(norm)
        histEvenEra[era].Scale(norm)

    ROOT.gROOT.cd('')
    histOdd = histOddEra['Run3_2022'].Clone('histOdd')
    histEven = histEvenEra['Run3_2022'].Clone('histEven')
    histOddIP = histOddIPEra['Run3_2022'].Clone('histOddIP')
    histEvenIP = histEvenIPEra['Run3_2022'].Clone('histEvenIP')

    for era in ['Run3_2022EE','Run3_2023','Run3_2023BPix']:
        histOdd.Add(histOdd,histOddEra[era])
        histEven.Add(histEven,histEvenEra[era])
        histOddIP.Add(histOddIP,histOddIPEra[era])
        histEvenIP.Add(histEvenIP,histEvenIPEra[era])

    
    normEven = histEven.GetSumOfWeights()
    normOdd  = histOdd.GetSumOfWeights()
    normEvenIP = histEvenIP.GetSumOfWeights()
    normOddIP  = histOddIP.GetSumOfWeights()

    utils.symmetrize(histEven)
    utils.symmetrize(histOdd)
    utils.symmetrize(histEvenIP)
    utils.symmetrize(histOddIP)
    
    print('')
    print('no IP  : Norm(CP-odd) %4.2f  :  Norm(CP-even) = %4.2f'%(normOdd,normEven))
    print('IP cut : Norm(CP-odd) %4.2f  :  Norm(CP-even) = %4.2f'%(normOddIP,normEvenIP))
    print('')
    print('  CP-Odd   |   CP-even')
    print('-----------+------------')
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

    suffix = ''
    if var=='aco_mu_pi' or var=='aco_e_pi':
        if noIpCutTau: suffix = '_no_ip2'
    
    hists = {}
    hists['even'] = histEven
    hists['odd'] = histOdd
    Plot(hists,var=var,channel=chan,asym=asym,suffix=suffix)

    suffix2 = suffix+'_ip'
    histsIP = {}
    histsIP['even'] = histEvenIP
    histsIP['odd'] = histOddIP
    Plot(histsIP,var=var,channel=chan,asym=asymIP,suffix=suffix2)
    
