#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# 
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

import yaml
from yaml.loader import SafeLoader



def Plot(h1,h2,**kwargs):

    var = kwargs.get('var','ip_x')
    chan = kwargs.get('chan','mt')
    
    # histograms

    xtitle = var+' [cm]'
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
    leg.AddEntry(h1,'electron','l')
    leg.AddEntry(h2,'pion','l')
    leg.Draw()

    styles.CMS_label(canvas,era='Run3_simulation',extraText='Simulation')
    
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
    parser.add_argument('-variable' ,'--variable', dest='variable', default='ip_x')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=40)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=-0.02)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=0.02)
    parser.add_argument('-dm','--dm',dest='dm',type=int, default=0)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    
    args = parser.parse_args()

    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    dm = args.dm

    sample = 'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'
    
    chan = 'et'
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

    if dm==0:
        cutDM = 'decayModePNet_2==0'
            
    
    cut = '%s&&%s&&%s&&%s'%(cutLep,cutTau,cutDeepTau,cutDM)
    dummy = ROOT.TCanvas('dummy','dummy',500,500)
    
    varToPlot1 = '1.25*('+var+'_1)'
    varToPlot2 = '0.75*('+var+'_2)'

    baseFolder = '%s/src/CPHiggs/Analysis/'%(os.getenv('CMSSW_BASE'))
    eras = ['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix']

    hist1 = {}
    hist2 = {}
    
    for era in eras:

        print('processing %s'%(era)) 
        yaml_file = baseFolder+'/params/'+era+'.yaml'
        metafile = open(yaml_file,'r')
        metadata = list(yaml.load_all(metafile,Loader=SafeLoader))
        norm = metadata[0]['lumi']*metadata[0][sample]['filter_efficiency']*metadata[0][sample]['xs']/metadata[0][sample]['eff']

        ROOT.gROOT.cd('')
        hist1[era] = ROOT.TH1D('hist1_%s'%(era),'',nbins,array('d',list(bins)))
        hist2[era] = ROOT.TH1D('hist2_%s'%(era), '',nbins,array('d',list(bins)))

        folder = '%s/%s/%s/%s'%(utils.tupleFolderV2,era,chan,sample)
        fileName  = '%s/nominal/merged.root'%(folder)
        inputFile = ROOT.TFile(fileName,'read')
        tree = inputFile.Get('ntuple')

        ROOT.gROOT.cd('')
        tree.Draw(varToPlot1+'>>hist1_%s'%(era),cut)
        tree.Draw(varToPlot2+'>>hist2_%s'%(era),cut)

        ROOT.gROOT.cd('')        
        hist1[era].Scale(norm)
        hist2[era].Scale(norm)

    ROOT.gROOT.cd('')
    h1 = hist1['Run3_2022'].Clone('h1')
    h2 = hist2['Run3_2022'].Clone('h2')

    for era in ['Run3_2022EE','Run3_2023','Run3_2023BPix']:
        h1.Add(h1,hist1[era])
        h2.Add(h2,hist2[era])
        
    Plot(h1,h2,var=var,chan=chan)
    
