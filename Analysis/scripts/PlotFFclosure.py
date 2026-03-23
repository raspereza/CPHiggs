#! /usr/bin/env python3
# Author: Alexei Raspereza (November 2025)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def ExtractHistosFF(f,var,region,bins):
    hists = {}

    for sample in utils.samples:
        for typ in utils.type_labels:
            nameInput='%s_%s_%s_%s'%(sample,var,region,typ)
            name='%s_%s'%(sample,typ)
            hists[name]=utils.rebinHisto(f.Get(nameInput),bins,'rebinned')
            for ff in utils.ff_labels:
                nameInput='%s_%s_%s_%s_%s'%(sample,var,region,ff,typ)
                name='%s_%s_%s'%(sample,ff,typ)
                hists[name] = utils.rebinHisto(f.Get(nameInput),bins,'rebinned')

    return hists
    
def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.501)
    ymax = kwargs.get('ymax',1.499)
    ff_label = kwargs.get('ff','qcd')
    region = kwargs.get('region','qcd_ff')
    
    mc_reduced_samples = ['zll','vv','wjets']
    h_mc = {} # sample, {sig,ar,qcd,wj}, {jet, lepton} 
    # histograms
    
    h_data = hists['data_all'].Clone('h_data')
    h_data_ar = hists['data_ar_all'].Clone('h_data_ar')
    h_data_ff = hists['data_%s_all'%(ff_label)].Clone('h_data_ff')
    
    h_mc_lep = hists['ztt_lep'].Clone('h_mc_lep')
    h_mc_lep.Add(h_mc_lep,hists['ztt_tau'])

    h_mc_ff_lep = hists['ztt_%s_lep'%(ff_label)].Clone('h_mc_ff_lep')
    h_mc_ff_lep.Add(h_mc_ff_lep,hists['ztt_%s_tau'%(ff_label)])

    h_mc_ar_lep = hists['ztt_ar_lep'].Clone('h_mc_ar_lep')
    h_mc_ar_lep.Add(h_mc_ar_lep,hists['ztt_ar_tau'])
    
    h_mc_ar_had = hists['ztt_ar_had'].Clone('h_mc_ar_had')
    h_mc_had = hists['ztt_had'].Clone('h_mc_had')

    h_top_lep = hists['top_lep'].Clone('h_top_lep')
    h_top_lep.Add(h_top_lep,hists['top_tau'])

    h_top_ff_lep = hists['top_%s_lep'%(ff_label)].Clone('h_top_ff_lep')
    h_top_ff_lep.Add(h_top_ff_lep,hists['top_%s_tau'%(ff_label)])

    h_top_ar_lep = hists['top_ar_lep'].Clone('h_top_ar_lep')
    h_top_ar_lep.Add(h_top_ar_lep,hists['top_ar_tau'])

    h_top_ar_had = hists['top_ar_had'].Clone('h_top_ar_had')
    h_top_had = hists['top_had'].Clone('h_top_had')
    
    for mc in mc_reduced_samples:
        for leptype in ['lep','tau']:
            h_mc_lep.Add(h_mc_lep,hists['%s_%s'%(mc,leptype)])
            h_mc_ar_lep.Add(h_mc_ar_lep,hists['%s_ar_%s'%(mc,leptype)])
            h_mc_ff_lep.Add(h_mc_ff_lep,hists['%s_%s_%s'%(mc,ff_label,leptype)])
        h_mc_ar_had.Add(h_mc_ar_had,hists['%s_ar_had'%(mc)])
        h_mc_had.Add(h_mc_had,hists['%s_had'%(mc)])
            
    
    h_data_ff.Add(h_data_ff,h_mc_ff_lep,1.,-1.)
    h_data_ff.Add(h_data_ff,h_top_ff_lep,1.,-1.)
    
    h_mc_ar_lep.Add(h_mc_ar_lep,h_top_ar_lep)
    h_qcd_ar = h_data_ar.Clone('h_qcd_ar')
    h_qcd_ar.Add(h_qcd_ar,h_mc_ar_lep,1.,-1.)
    h_qcd_ar.Add(h_qcd_ar,h_mc_ar_had,1.,-1.)
    h_qcd_ar.Add(h_qcd_ar,h_top_ar_had,1.,-1.)
    utils.removeNegativeBins(h_qcd_ar)

    h_mc_lep.Add(h_mc_lep,h_top_lep)
    h_qcd = h_data.Clone('h_qcd')
    h_qcd.Add(h_qcd,h_mc_lep,1.,-1.)
    h_qcd.Add(h_qcd,h_mc_had,1.,-1.)
    h_qcd.Add(h_qcd,h_top_had,1.,-1.)
    utils.removeNegativeBins(h_qcd)
    
    xtitle = utils.XTitle[chan][var]
    styles.InitHist(h_mc_lep,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_data_ff,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)

    x_data = h_data.GetSumOfWeights()
    x_lep = h_mc_lep.GetSumOfWeights()
    x_fakes = h_data_ff.GetSumOfWeights()
    x_tot = x_lep+x_fakes
    
    print('')
    print('Yields ->')
    print('Leptons : %7.0f'%(x_lep))
    print('Total   : %7.0f'%(x_tot))
    print('Data    : %7.0f'%(x_data))
    print('')

    h_data_ff.Add(h_data_ff,h_mc_lep)
                
    h_tot = h_data_ff.Clone("total")
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_data_ff)
    utils.zeroBinErrors(h_mc_lep)

    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    h_data.GetYaxis().SetRangeUser(0.,1.2*YMax)
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
    h_data_ff.Draw('hsame')
    h_mc_lep.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.72,0.53,0.85,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.AddEntry(h_data,'data','lp')
    leg.AddEntry(h_data_ff,'jet#rightarrow#tau','f')
    leg.AddEntry(h_mc_lep,'leptons','f')
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
    outputFolder = '/eos/home-r/rasp/php-plots/plots/FFclosure_%s/%s'%(chan,region)
    outputGraphics = '%s/%s.png'%(outputFolder,var)
    canvas.Print(outputGraphics)

    ######################
    ###### denominator ###
    ######################
    utils.zeroBinErrors(h_top_ar_had)
    utils.zeroBinErrors(h_qcd_ar)
    utils.zeroBinErrors(h_mc_ar_lep)
    utils.zeroBinErrors(h_mc_ar_had)

    h_top_ar_had.Add(h_top_ar_had,h_mc_ar_lep)
    h_mc_ar_had.Add(h_mc_ar_had,h_top_ar_had)
    h_qcd_ar.Add(h_qcd_ar,h_mc_ar_had)

    styles.InitHist(h_mc_ar_lep,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_mc_ar_had,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)    
    styles.InitHist(h_top_ar_had,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_qcd_ar,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)
    h_qcd_ar.GetXaxis().SetTitle(xtitle)
    h_qcd_ar.GetYaxis().SetTitle('Events')
    h_qcd_ar.GetYaxis().SetTitleOffset(1.4)
    
    canvas = styles.MakeCanvas("canv_num","",600,600)
    
    h_qcd_ar.Draw('hsame')
    h_mc_ar_had.Draw('hsame')
    h_top_ar_had.Draw('hsame')
    h_mc_ar_lep.Draw('hsame')

    leg = ROOT.TLegend(0.7,0.55,0.85,0.85)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.SetHeader('den %s'%ff_label)
    leg.AddEntry(h_qcd_ar,'QCD','f')
    leg.AddEntry(h_mc_ar_had,'WJ j#rightarrow#tau','f')
    leg.AddEntry(h_top_ar_had,'Top j#rightarrow#tau','f')
    leg.AddEntry(h_mc_ar_lep,'leptons','f')
    if plotLegend: leg.Draw()

    canvas.Update()
    print('')
    outputFolder = '/eos/home-r/rasp/php-plots/plots/FFcomposition_%s/%s'%(chan,region)
    outputGraphics = '%s/%s_den.png'%(outputFolder,var)
    canvas.Print(outputGraphics)

    ######################
    ###### numerator #####
    ######################
    utils.zeroBinErrors(h_top_had)
    utils.zeroBinErrors(h_qcd)
    utils.zeroBinErrors(h_mc_lep)
    utils.zeroBinErrors(h_mc_had)

    h_top_had.Add(h_top_had,h_mc_lep)
    h_mc_had.Add(h_mc_had,h_top_had)
    h_qcd.Add(h_qcd,h_mc_had)

    styles.InitHist(h_mc_lep,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_mc_had,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)    
    styles.InitHist(h_top_had,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_qcd,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)
    h_qcd.GetXaxis().SetTitle(xtitle)
    h_qcd.GetYaxis().SetTitle('Events')
    h_qcd.GetYaxis().SetTitleOffset(1.4)
    
    canvas = styles.MakeCanvas("canv","",600,600)
    
    h_qcd.Draw('hsame')
    h_mc_had.Draw('hsame')
    h_top_had.Draw('hsame')
    h_mc_lep.Draw('hsame')

    leg = ROOT.TLegend(0.7,0.55,0.85,0.85)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.045)
    leg.SetHeader('num %s'%ff_label)
    leg.AddEntry(h_qcd,'QCD','f')
    leg.AddEntry(h_mc_had,'WJ j#rightarrow#tau','f')
    leg.AddEntry(h_top_had,'Top j#rightarrow#tau','f')
    leg.AddEntry(h_mc_lep,'leptons','f')
    if plotLegend: leg.Draw()

    canvas.Update()
    print('')
    outputFolder = '/eos/home-r/rasp/php-plots/plots/FFcomposition_%s/%s'%(chan,region)
    outputGraphics = '%s/%s_num.png'%(outputFolder,var)
    canvas.Print(outputGraphics)

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
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.899)
    parser.add_argument('-suffix','--suffix',dest='suffix',default='x_ipcut1_ff_ipcut')
    #    parser.add_argument('-suffix','--suffix',dest='suffix',default='x_ff')
    parser.add_argument('-ff','--ff',dest='ff',default='wj')
    parser.add_argument('-region','--region',dest='region',default='wj_ff')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    var = args.variable
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    suffix = args.suffix
    plotLegend = True
    ff = args.ff
    region = args.region
    no_legend_vars = ['eta_1','eta_2','CMetQCD','CMetW','bdt_ditau','bdt_fakes','bdt_signal','dR','mt_1','aco_lep_pi','aco_lep_rho','aco_lep_a1_1pr','aco_lep_a1_3pr'] 
    if var in no_legend_vars:
        plotLegend = False

    basedir = utils.outputFolder+'/selection'

    bins = utils.createBins(nbins,xmin,xmax)

    inputFileName = '%s/baseline/%s_%s_%s.root'%(basedir,chan,era,suffix)
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
    hists = ExtractHistosFF(inputFile,var,region,bins)
    Plot(hists,era=era,var=var,channel=chan,
         ymin=ymin,ymax=ymax,region=region,
         plotLegend=plotLegend,ff=ff)

    

