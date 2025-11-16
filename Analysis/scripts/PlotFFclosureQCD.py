#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2025)
# QCD closure
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def ExtractHistosFF(f,var,region,ff,bins):
    hists = {}
    for sample in utils.samples:
        for typ in utils.type_labels:
            nameInput='%s_%s_%s_%s'%(sample,var,region,typ)
            name='%s_%s'%(sample,typ)
            hists[name]=utils.rebinHisto(f.Get(nameInput),bins,'rebinned')
            nameInput='%s_%s_%s_%s_%s'%(sample,var,region,ff,typ)
            name='%s_model_%s'%(sample,typ)
            hists[name] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')

    return hists
    
def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    region = kwargs.get('region','lowmt_os_antiiso')
    ff = kwargs.get('ff','ss_antiiso')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.501)
    ymax = kwargs.get('ymax',1.499)
    ff_version = kwargs.get('ff_version','v4')
    
    h_data = hists['data_all']
    h_model = hists['data_model_all']
    h_mc_lep = utils.ReplicaHist(h_data,'h_mc_lep')
    h_mc_had = utils.ReplicaHist(h_data,'h_mc_had')
    x_mc_direct = 0
    x_mc_invert = 0
    for s in utils.bkg_samples:
        nameMC = '%s_had'%(s)
        h_mc_had.Add(h_mc_had,hists[nameMC],1.,1.)
        for typ in ['lep','tau']:
            nameMC = '%s_%s'%(s,typ)
            h_data.Add(h_data,hists[nameMC],1.,-1.)
            h_mc_lep.Add(h_mc_lep,hists[nameMC],1.,1.)
            nameMC = '%s_model_%s'%(s,typ)
            h_model.Add(h_model,hists[nameMC],1.,-1.)
    

    styles.InitData(h_data)
    xtitle = utils.XTitle[chan][var]
    styles.InitHist(h_model,xtitle,"Events",ROOT.TColor.GetColor("#FFCCFF"),1001)

    x_data = h_data.GetSumOfWeights()
    x_model = h_model.GetSumOfWeights()
    x_mc_lep = h_mc_lep.GetSumOfWeights()
    x_mc_had = h_mc_had.GetSumOfWeights()
    
    print('')
    print('Yields ->')
    print('Model   : %7.0f'%(x_model))
    print('MC lep  : %7.0f'%(x_mc_lep))
    print('MC had  : %7.0f'%(x_mc_had))
    print('Data    : %7.0f'%(x_data))
    print('')

    h_tot = h_model.Clone("total")
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_model)

    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    h_data.GetYaxis().SetRangeUser(0.,1.5*YMax)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("events")
#    h_ratio.GetXaxis().SetTitleSize(0.05)
#    h_ratio.GetYaxis().SetTitleSize(0.05)
    h_ratio.GetYaxis().SetTitle("obs/exp")
    h_ratio.GetXaxis().SetTitle(xtitle)


    #######################
    #### Closure tests ####
    #######################
    
    # canvas and pads
    canvas = styles.MakeCanvas("canv_den","",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    
    
    h_data.Draw('e1')
    h_model.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.72,0.53,0.85,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.SetHeader('QCD closure')
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_model,'FF model','f')
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
    print('')
    outputFolder = '/eos/home-r/rasp/php-plots/plots/FFclosure_%s/%s/%s'%(chan,ff_version,region)
    outputGraphics = '%s/FF_%s/%s.png'%(outputFolder,ff,var)
    canvas.Print(outputGraphics)
    return h_data,h_model
    
if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3',choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-variable','--variable',dest='variable',default='m_vis')
    parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et'])
    parser.add_argument('-nbins','--nbins', dest='nbins',type=int,default=50)
    parser.add_argument('-xmin','--xmin',dest='xmin',type=float,default=0.0)
    parser.add_argument('-xmax','--xmax',dest='xmax',type=float,default=250.)
    parser.add_argument('-ymin','--ymin',dest='ymin',type=float,default=0.201)
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.799)
    parser.add_argument('-ff_version','--ff_version',dest='ff_version',default='v4')
    parser.add_argument('-region','--region',dest='region',default='lowmt_os_antiiso')
    parser.add_argument('-ff','--ff',dest='ff',default='ss_antiiso')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax

    plotLegend = True
    if var in ['eta_1','eta_2','CMetQCD','CMetW','bdt_ditau','bdt_fakes','bdt_signal']:
        plotLegend = False

    ff = args.ff
    region = args.region
    ff_version = args.ff_version
    
    basedir = utils.outputFolder+'/selection/baseline'

    bins = utils.createBins(nbins,xmin,xmax)

    inputFileName = '%s/%s_%s_x_ipcut1_ff_%s.root'%(basedir,chan,era,ff_version)
    if os.path.isfile(inputFileName):
        print('')
        print('Loading ROOT file %s'%(inputFileName))
        print('')
    else:
        print('')
        print('ROOT file %s not found'%(inputFileName))
        print('quitting')
        print('')
    inputFile = ROOT.TFile(inputFileName,'read')
    hists = ExtractHistosFF(inputFile,var,region,ff,bins)
    h_data,h_model = Plot(hists,era=era,var=var,channel=chan,
                          ymin=ymin,ymax=ymax,
                          plotLegend=plotLegend,
                          region=region,
                          ff=ff,
                          ff_version=ff_version)
    
    

