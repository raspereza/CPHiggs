#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

def Plot(f,d,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    chan = kwargs.get('channel','mt')
    binPt = kwargs.get('binPt','1')
    binEta = kwargs.get('binEta','1')
    postFit = kwargs.get('postFit',False)
    region = kwargs.get('region','pass')
    per_process = kwargs.get('per_process',True)
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.501)
    ymax = kwargs.get('ymax',1.499)
    
    typ_suffix = 'typ'
    if per_process: typ_suffix = 'proc'
    fit_suffix = 'prefit'
    if postFit: fit_suffix = 'postfit' 

    folder_fit = 'shapes_prefit'
    if postFit: folder_fit='shapes_fit_s'

    folder_region = 'ch2'
    if region=='pass':
        folder_region='ch1'

    folder = '%s/%s'%(folder_fit,folder_region)
    # histograms
    h_data = d.Get(region+'/data_obs').Clone('h_data')

    h_ztt = h_data.Clone('h_ztt')
    h_zll = h_data.Clone('h_zll')
    h_ttt = h_data.Clone('h_ttt')
    h_ttl = h_data.Clone('h_ttl')
    h_vvt = h_data.Clone('h_vvt')
    h_vvl = h_data.Clone('h_vvl')
    h_wjets = h_data.Clone('h_wjets')
    h_qcd = h_data.Clone('h_qcd')
    
    ZTT = f.Get(folder+'/ZTT_'+region)
    TTT = f.Get(folder+'/TTT_'+region)
    VVT = f.Get(folder+'/VVT_'+region)

    ZLL = f.Get(folder+'/ZLL').Clone('h_zll')
    TTL = f.Get(folder+'/TTL').Clone('h_ttl')
    VVL = f.Get(folder+'/VVL').Clone('h_vvl')

    WJ = f.Get(folder+'/WJ').Clone('h_zll')
    QCD = f.Get(folder+'/QCD').Clone('h_qcd')

    utils.copyHist(ZTT,h_ztt)
    utils.copyHist(ZLL,h_zll)
    utils.copyHist(TTT,h_ttt)
    utils.copyHist(TTL,h_ttl)
    utils.copyHist(VVT,h_vvt)
    utils.copyHist(VVL,h_vvl)
    utils.copyHist(WJ,h_wjets)
    utils.copyHist(QCD,h_qcd)
    
    h_lep = h_zll.Clone('h_lep')
    h_tau = h_ztt.Clone('h_tau')

    h_lep.Add(h_lep,h_ttl,1.,1.)
    h_lep.Add(h_lep,h_vvl,1.,1.)

    h_tau.Add(h_tau,h_ttt,1.,1.)
    h_tau.Add(h_tau,h_vvt,1.,1.)

    h_top = h_ttt.Clone('h_top')
    h_vv = h_vvt.Clone('h_vv')

    h_top.Add(h_top,h_ttl,1.,1.)
    h_vv.Add(h_vv,h_vvl,1.,1.)
    
    styles.InitData(h_data)

    xtitle = 'm_{vis} (GeV)'
    styles.InitHist(h_ztt,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_vv,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    styles.InitHist(h_qcd,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)

    styles.InitHist(h_tau,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_lep,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_wjets,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)

    x_data = h_data.GetSumOfWeights()
    
    x_ztt   = h_ztt.GetSumOfWeights()
    x_zll   = h_zll.GetSumOfWeights() 
    x_top   = h_top.GetSumOfWeights()
    x_vv    = h_vv.GetSumOfWeights() 
    x_wjets = h_wjets.GetSumOfWeights()
    x_qcd   = h_qcd.GetSumOfWeights() 
    x_tot   = x_ztt + x_zll + x_top + x_vv + x_wjets + x_qcd

    x_tau = h_tau.GetSumOfWeights()
    x_lep = h_lep.GetSumOfWeights()
        
    print('')
    print('%s : %s : binPt%s_binEta%s : %s : %s ->'%(era,chan,binPt,binEta,region,fit_suffix))
    if per_process:
        print('Ztautau    : %7.0f'%(x_ztt))
        print('Zll        : %7.0f'%(x_zll))
        print('Top        : %7.0f'%(x_top))
        print('VV         : %7.0f'%(x_vv))
        print('WJets      : %7.0f'%(x_wjets))
        print('QCD        : %7.0f'%(x_qcd))
    else:
        print('Prompt lep : %7.0f'%(x_lep))
        print('tau -> lep : %7.0f'%(x_tau))
        print('jet -> tau : %7.0f'%(x_wjets))
        print('QCD        : %7.0f'%(x_qcd))
        x_tot = x_tau + x_lep + x_wjets + x_qcd
        
    print('Total      : %7.0f'%(x_tot))
    print('Data       : %7.0f'%(x_data))
    print('')

    h_vv.Add(h_vv,h_wjets,1.,1.)
    h_vv.Add(h_vv,h_qcd,1.,1.)
    h_top.Add(h_top,h_vv,1.,1.)
    h_zll.Add(h_zll,h_top,1.,1.)
    h_ztt.Add(h_ztt,h_zll,1.,1.)

    h_wjets.Add(h_wjets,h_qcd,1.,1.)
    h_lep.Add(h_lep,h_wjets,1.,1.)
    h_tau.Add(h_tau,h_lep,1.,1.)
    
    h_tot = None
    if per_process: h_tot = h_ztt.Clone("total")
    else: h_tot = h_tau.Clone("total")

    print('Total (cross-check) : %7.0f'%(h_tot.GetSumOfWeights()))
    print('')
    
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_ztt)
    utils.zeroBinErrors(h_zll)
    utils.zeroBinErrors(h_top)
    utils.zeroBinErrors(h_vv)
    utils.zeroBinErrors(h_wjets)
    utils.zeroBinErrors(h_qcd)

    utils.zeroBinErrors(h_tau)
    utils.zeroBinErrors(h_lep)
    
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
    if per_process:
        h_ztt.Draw('hsame')
        h_zll.Draw('hsame')
        h_top.Draw('hsame')
        h_vv.Draw('hsame')
        h_qcd.Draw('hsame')
    else:
        h_tau.Draw('hsame')
        h_lep.Draw('hsame')
        h_wjets.Draw('hsame')
        h_qcd.Draw('hsame')
        
    h_data.Draw('e1same')
    h_tot.Draw('e2same')

    leg = ROOT.TLegend(0.65,0.4,0.85,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.042)
    leg.AddEntry(h_data,'data','lp')
    if per_process:
        leg.AddEntry(h_ztt,'Z#rightarrow#tau#tau','f')
        if chan=='mt': leg.AddEntry(h_zll,'Z#rightarrow#mu#mu','f')
        else: leg.AddEntry(h_zll,'Z#rightarrowee','f')
        leg.AddEntry(h_top,'t#bar{t}','f')
        leg.AddEntry(h_vv,'electroweak','f')
        leg.AddEntry(h_qcd,'QCD','f')
    else:
        if chan=='mt':
            leg.AddEntry(h_tau,'#tau#rightarrow#mu','f')
            leg.AddEntry(h_lep,'prompt #mu','f')
        else:
            leg.AddEntry(h_tau,'#tau#rightarrowe','f')
            leg.AddEntry(h_lep,'prompt e','f')
        leg.AddEntry(h_wjets,'simulated fakes','f')
        leg.AddEntry(h_qcd,'QCD','f')
        
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
    
    binPtEta = 'binPt%s_binEta%s'%(binPt,binEta)
    basedir = '%s/src/CPHiggs/IP/figures'%(os.getenv('CMSSW_BASE'))
    outputGraphics = '%s/TP_%s_%s_%s_%s_%s_%s.png'%(basedir,chan,era,binPtEta,region,fit_suffix,typ_suffix)    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-perType','--perType', dest='perType', action='store_true')
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.501)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.499)
    parser.add_argument('-postFit','--postFit', dest='postFit',action='store_true')
    
    args = parser.parse_args()
    era = args.era
    chan = args.channel
    ymin = args.ymin
    ymax = args.ymax

    plotLegend = True
    postFit = args.postFit
    
    proc = True
    if args.perType: proc = False
    
    nbinsPt = 4
    nbinsEta = 2
    if chan=='mt': nbinsPt = 5
    region_labels = ['pass','fail']
    basedir='%s/src/CPHiggs/IP/datacards'%(os.getenv('CMSSW_BASE'))
    for iPt in range(1,nbinsPt+1):
        binPt='%1i'%(iPt)
        for iEta in range(1,nbinsEta+1):
            binEta='%1i'%(iEta)
            fileName = '%s/%s_%s_binPt%s_binEta%s_fit.root'%(basedir,chan,era,binPt,binEta)            
            inputFile = ROOT.TFile(fileName,'read')
            for region in region_labels:
                datacardsFileName = '%s/%s_%s_binPt%s_binEta%s_%s.root'%(basedir,chan,era,binPt,binEta,region)
                datacardsFile = ROOT.TFile(datacardsFileName,'read')
                Plot(inputFile,
                     datacardsFile,
                     era=era,
                     channel=chan,
                     per_process=proc,
                     region=region,
                     binPt=binPt,
                     binEta=binEta,
                     postFit=postFit)
                datacardsFile.Close()

            inputFile.Close()
    
