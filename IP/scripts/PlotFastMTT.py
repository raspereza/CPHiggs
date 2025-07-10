#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro to test FastMTT
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

def Plot(hist1,hist2,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    chan = kwargs.get('chan','mt')
    isDY = kwargs.get('isDY',True) 
    isMass = kwargs.get('isMass',True)
    
    # histograms
    h_even = hist1
    h_odd = hist2

    xtitle = 'p_{T}^{rec}/p_{T}^{gen}'
    leg_even = 'w/o m_{H}'
    leg_odd = 'with m_{H}'
    if isMass:
        xtitle = 'mass (GeV)'
        leg_even = 'm_{vis}'
        leg_odd  = 'm_{#tau#tau}'

    header = 'H#rightarrow #tau_{#mu}#tau_{h}'
    if chan=='tt':
        header = 'H#rightarrow#tau_{h}#tau_{h}'
    if isDY:
        header = 'Z#rightarrow#tau_{#mu}#tau_{h}'
        if chan=='tt':
            header = 'Z#rightarrow#tau_{h}#tau_{h}'
        
    ytitle = 'normalised'

    styles.InitModel(h_even,xtitle,ytitle,2)
    styles.InitModel(h_odd,xtitle,ytitle,4)

    utils.zeroBinErrors(h_even)
    utils.zeroBinErrors(h_odd)
    
    YMax = h_even.GetMaximum()
    if h_odd.GetMaximum()>YMax: YMax = h_odd.GetMaximum()
    
    h_even.GetYaxis().SetRangeUser(0.,1.1*YMax)

    # canvas and pads
    canv_name = 'canv_pt'
    if isMass: canv_name = 'canv_mass'
    canvas = styles.MakeCanvas(canv_name,"",800,700)
    
    h_even.Draw('h')
    h_odd.Draw('hsame')

    leg = ROOT.TLegend(0.7,0.5,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetHeader(header)
    leg.SetTextSize(0.05)
    leg.AddEntry(h_even,leg_even,'l')
    leg.AddEntry(h_odd,leg_odd,'l')
    leg.Draw()

    styles.CMS_label(canvas,era=era,extraText='Simulation')

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    name_var = 'pt'
    if isMass:
        name_var = 'mass'
    name_sample = 'Higgs'
    if isDY:
        name_sample = 'DY'
    outputGraphics = os.getenv('CMSSW_BASE')+'/src/CPHiggs/IP/figures/FastMTT_'+name_var+'_'+name_sample+'_'+chan+'_'+era+'.png'
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022EE', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et','tt'])
    parser.add_argument('-isDY','--isDY',dest='isDY',action='store_true')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    isDY = args.isDY
    
    sample = 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'
    if isDY:
        sample = 'DYto2L_M_50_madgraphMLM'
    
    # mass plot binning
    nbins = 50
    xmin = 0.
    xmax = 250.
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    # pT resolution plot binning 
    xmin_pt = 0.
    xmax_pt = 2.5
    nbins_pt = 50
    bins_pt = []
    width_pt = (xmax_pt-xmin_pt)/float(nbins_pt)
    for i in range(0,nbins_pt+1):
        xb = xmin + width_pt*float(i)
        bins_pt.append(xb)
    
    folder = '/eos/cms/store/group/phys_tau/ksavva/For_Aliaksei/files/testingzpt' 
        
    filename = '%s/%s/%s/%s/nominal/merged.root'%(folder,era,chan,sample)
    inputfile = ROOT.TFile(filename,'read')
    tree = inputfile.Get('ntuple')
    
    var_mvis = 'm_vis'
    var_mtt  = 'FastMTT_mass'
    var_pt1 = '(FastMTT_pt_1/genPart_pt_1)'
    var_pt2 = '(FastMTT_pt_2/genPart_pt_2)'
    var_pt1_cons = '(FastMTT_pt_1_constraint/genPart_pt_1)'
    var_pt2_cons = '(FastMTT_pt_2_constraint/genPart_pt_2)'
    
    cuts_tau2 = 'idDeepTau2018v2p5VSe_2>=6&&idDeepTau2018v2p5VSmu_2>=4&&idDeepTau2018v2p5VSjet_2>=5'
    cuts_tau2 += '&&pt_2>20.0&&fabs(eta_2)<2.3'
    cuts_tau1 = ''
    cuts_id = ''
    if chan=='mt':
        cuts_id = '(genPart_pdgId_1==13||genPart_pdgId_1==-13)&&(genPart_pdgId_2==15||genPart_pdgId_2==-15)'
        cuts_tau1 = 'iso_1<0.2'
        cuts_tau1 += '&&pt_1>25.0&&fabs(eta_1)<2.4'
    else:
        cuts_id = '(genPart_pdgId_1==15||genPart_pdgId_1==-15)&&(genPart_pdgId_2==15||genPart_pdgId_2==-15)'
        cuts_tau1 = 'idDeepTau2018v2p5VSe_1>=6&&idDeepTau2018v2p5VSmu_1>=4&&idDeepTau2018v2p5VSjet_1>=5' 
        cuts_tau1 += '&&pt_2>20.0&&fabs(eta_2)<2.3'

    cuts = cuts_id+'&&'+cuts_tau1+'&&'+cuts_tau2

    inputfile = ROOT.TFile(filename,'read')
    tree = inputfile.Get('ntuple')

    # declaring histograms
    
    h_mvis = ROOT.TH1D('h_mvis','',nbins,array('d',list(bins)))
    h_mtt = ROOT.TH1D('h_mtt','',nbins,array('d',list(bins)))

    h_pt1 = ROOT.TH1D('h_pt1','',nbins_pt,array('d',list(bins_pt)))
    h_pt1_cons = ROOT.TH1D('h_pt1_cons','',nbins_pt,array('d',list(bins_pt)))
    
    h_pt2 = ROOT.TH1D('h_pt2','',nbins_pt,array('d',list(bins_pt)))
    h_pt2_cons = ROOT.TH1D('h_pt2_cons','',nbins_pt,array('d',list(bins_pt)))
    
    tree.Draw(var_mvis+'>>h_mvis',cuts)
    tree.Draw(var_mtt+'>>h_mtt',cuts)
    tree.Draw(var_pt2+'>>h_pt2',cuts)
    tree.Draw(var_pt2_cons+'>>h_pt2_cons',cuts)
    if chan=='tt':
        tree.Draw(var_pt1+'>>h_pt1',cuts)
        tree.Draw(var_pt1_cons+'>>h_pt1_cons',cuts)
        h_pt2.Add(h_pt2,h_pt1)
        h_pt2_cons.Add(h_pt2_cons,h_pt1_cons)


    # normlization
    x_mvis = 0.
    x_mtt = 0.
    for i in range(1,nbins+1):
        x_mvis += h_mvis.GetBinContent(i)
        x_mtt += h_mtt.GetBinContent(i)

    x_pt2 = 0.
    x_pt2_cons = 0.
    for i in range(1,nbins_pt+1):
        x_pt2 += h_pt2.GetBinContent(i)
        x_pt2_cons += h_pt2_cons.GetBinContent(i)
    
    h_mvis.Scale(1./x_mvis)
    h_mtt.Scale(1./x_mtt)
    h_pt2.Scale(1./x_pt2)
    h_pt2_cons.Scale(1./x_pt2_cons)

    mean_pt = h_pt2.GetMean()
    rms_pt = h_pt2.GetRMS()

    mean_pt_cons = h_pt2_cons.GetMean()
    rms_pt_cons = h_pt2_cons.GetRMS()

    print('     :   w/o mH : with mH  ')
    print('Mean :   %5.3f  :  %5.3f'%(mean_pt,mean_pt_cons))
    print('RMS  :   %5.3f  :  %5.3f'%(rms_pt,rms_pt_cons))
    print('')

    Plot(h_mvis,h_mtt,era=era,isDY=isDY,chan=chan,isMass=True)
    Plot(h_pt2,h_pt2_cons,era=era,isDY=isDY,chan=chan,isMass=False)
    
    
