#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->ll  selection
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

# extracting histos from ROOT file created by RunSelection.py 
def extractHistos(f,var,bins):

    hists = {}
    
    for sample in utils.samples:
        nameInput = '%s_%s_incl_os_iso_all'%(sample,var)
        hists[sample] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')
        
    name = 'zll_incl'
    nameout = 'zll'
    hists[nameout] = hists[name].Clone(hists[name].GetName()+'_zll')
    for s in ['zll_ext','zll_0j','zll_1j','zll_2j']:
        hists[nameout].Add(hists[nameout],hists[s],1.,1.)
        
    name = 'ztt_0j'
    nameout = 'ztt'
    hists[nameout] = hists[name].Clone(hists[name].GetName()+'_ztt')
    for s in ['ztt_1j','ztt_2j']:
        hists[nameout].Add(hists[nameout],hists[s],1.,1.)
        
    return hists
    
def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.5)
    ymax = kwargs.get('ymax',1.5)
    suffix = kwargs.get('suffix','_x')
    calibrDY = kwargs.get('calibrDY',1.0)
    
    # histograms
    h_data = hists['data'].Clone('h_data')
    h_zll = hists['zll'].Clone('h_zll')
    h_zll.Add(h_zll,hists['ztt'],1.,1.)
    h_top = hists['top'].Clone('h_top')
    h_vv = hists['vv'].Clone('h_vv')
    h_wjets = hists['wjets'].Clone('h_wjets')

    styles.InitData(h_data)

    xtitle = utils.XTitle[chan][var]
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_vv,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    styles.InitHist(h_wjets,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    
    x_data = h_data.GetSumOfWeights()
    
    x_zll   = h_zll.GetSumOfWeights() 
    x_top   = h_top.GetSumOfWeights()
    x_vv    = h_vv.GetSumOfWeights() 
    x_wjets = h_wjets.GetSumOfWeights()
    x_tot   = x_zll + x_top + x_vv + x_wjets
    
    scale = 1.0 # (x_data-x_top-x_wjets-x_vv)/x_zll
    
    h_zll.Scale(calibrDY)
    x_zll   = h_zll.GetSumOfWeights() 
    x_tot   = x_zll + x_top + x_vv + x_wjets
    
    print('')
    print('Yields ->')
    print('DY    : %7.0f'%(x_zll))
    print('Top   : %7.0f'%(x_top))
    print('VV    : %7.0f'%(x_vv))
    print('WJets : %7.0f'%(x_wjets))

    x_tot_lumi = x_tot/utils.eraLumi[era]
    x_data_lumi = x_data/utils.eraLumi[era]
    print('')
    print('Lumi = %3.0f'%(utils.eraLumi[era]))
    print('Total : %7.0f  ->  per pb-1 : %5.3f'%(x_tot,x_tot_lumi))
    print('Data  : %7.0f  ->  per pb-1 : %5.3f'%(x_data,x_data_lumi))
    print('')

    h_vv.Add(h_vv,h_wjets,1.,1.)
    h_top.Add(h_top,h_vv,1.,1.)
    h_zll.Add(h_zll,h_top,1.,1.)

    h_tot = h_zll.Clone("total")
    styles.InitTotalHist(h_tot)

    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_zll)
    utils.zeroBinErrors(h_top)
    utils.zeroBinErrors(h_vv)
    utils.zeroBinErrors(h_wjets)
    
    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    
    h_data.GetYaxis().SetRangeUser(0.,1.2*YMax)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("events")
    h_ratio.GetYaxis().SetTitle("obs/exp")
    h_ratio.GetXaxis().SetTitle(xtitle)

    # canvas and pads
    canvas = styles.MakeCanvas("canv","",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    
    
    h_data.Draw('e1')
    h_zll.Draw('hsame')
    h_top.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.7,0.45,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_zll,'Drell-Yan','f')
    leg.AddEntry(h_top,'rest','f')
    if plotLegend: leg.Draw()

    styles.CMS_label(upper,era=era)

    upper.Draw("SAME")
    upper.RedrawAxis()
    upper.Modified()
    upper.Update()
    canvas.cd()

    # lower pad
    lower = ROOT.TPad("lower", "pad",0,0,1,0.30)
    lower.Draw()
    lower.cd()
    styles.InitLowerPad(lower)

    h_ratio.Draw('e1')
    h_tot_ratio.Draw('e2same')
    h_ratio.Draw('e1same')

    nbins = h_ratio.GetNbinsX()
    xmin = h_ratio.GetXaxis().GetBinLowEdge(1)    
    xmax = h_ratio.GetXaxis().GetBinLowEdge(nbins+1)
    line = ROOT.TLine(xmin,1.,xmax,1.)
    line.SetLineStyle(1)
    line.SetLineWidth(2)
    line.SetLineColor(4)
    line.Draw()
    lower.Modified()
    lower.RedrawAxis()

    canvas.cd()
    canvas.Modified()
    canvas.cd()
    canvas.SetSelected(canvas)
    canvas.Update()
    print('variable %s is plotted,   DY scale factor = %6.4f'%(var,scale))
    basedir = '%s/figures'%(utils.outputFolder)
    outputGraphics = '%s/%s_%s_%s%s.png'%(basedir,var,chan,era,suffix)
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3_2022',choices=['Run3_2022','Run3_2023','Run3_2022preEE','Run3_2022postEE','Run3_2023preBPix','Run3_2023postBPix','Run3'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mm','ee'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='m_vis')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=40)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=50.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=130.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    parser.add_argument('-applyIP1','--applyIP1',dest='applyIP1',type=int,default=0)
    parser.add_argument('-applyIP2','--applyIP2',dest='applyIP2',type=int,default=0)
    parser.add_argument('-applySF', '--applySF', dest='applySF' ,type=int,default=0)
    parser.add_argument('-applyJsonSF', '--applyJsonSF', dest='applyJsonSF',type=int,default=1)
    parser.add_argument('-generator', '--generator', dest='generator', default='amcatnlo',choices=['amcatnlo','MG','powheg'])
    parser.add_argument('-analysisType', '--analysisType', dest='analysisType', default='baseline',choices=['baseline','ipSig','datacardsPhiCP'])
    parser.add_argument('-calibrDY','--calibrDY',dest='calibrDY',type=float,default=1.0) # calibrate DY

    args = parser.parse_args()

    calibrDY = args.calibrDY
    
    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    generator = args.generator
    analysisType = args.analysisType
    
    plotLegend = True
    if var=='eta_1' or var=='eta_2':
        plotLegend = False

    
    applyIPSigLep1Cut = args.applyIP1
    applyIPSigLep2Cut = args.applyIP2
    applyIPSigPromptLepSF = args.applySF
    applyIPSigJsonSF = args.applyJsonSF
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)
        
    basedir = utils.outputFolder+'/selection/'+analysisType

    suffix = '_x'
    if applyIPSigLep1Cut==1:
        suffix += '_ipcut1'
    if applyIPSigLep2Cut==1:
        suffix += '_ipcut2'
    if applyIPSigPromptLepSF==1:
        suffix += '_promptSF'
    if applyIPSigJsonSF==1:
        suffix += '_json'

    inputFileName = '%s/%s_%s%s.root'%(basedir,chan,era,suffix)
    if os.path.isfile(inputFileName):        
        print('')
        print('Loading histograms from file %s'%(inputFileName))
        print('')
    else:
        print('')
        print('Input ROOT file %s not found'%(inputFileName))
        print('Quitting')
        print('')
        exit()
    inputFile = ROOT.TFile(inputFileName,'read')
    hists = utils.extractHistos(inputFile,var,bins,generator,era,chan)
    suffixOut = suffix + '_' + generator
    Plot(hists,era=era,var=var,channel=chan,ymin=ymin,ymax=ymax,plotLegend=plotLegend,suffix=suffixOut,calibrDY=calibrDY)

    
