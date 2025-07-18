#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro for Z->tautau  selection
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

XTitle = {
    'mt': {
        'mt_1'  : "m_{T} (GeV)",
        'pt_1'  : "muon p_{T} (GeV)",
        'eta_1' : "muon #eta",
        'pt_2'  : "tau p_{T} (GeV)",
        'eta_2' : "tau #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{vis} (GeV)",
        'ipsig_1': "muon IP sig",
        'ipsig_2': "tau IP sig",
        'n_jets': "number of jets",
        'n_bjets': "number of b-jets",
        'jpt_1': "leading jet p_{T} (GeV)",
        'jpt_2': "trailing jet p_{T} (GeV)",
        'mjj': "dijet mass (GeV)",
        'jdeta': "#Delta#eta(j,j)"
    },
    'et': {
        'mt_1'  : "m_{T} (GeV)",
        'pt_1'  : "electron p_{T} (GeV)",
        'eta_1' : "electron #eta",
        'pt_2'  : "tau p_{T} (GeV)",
        'eta_2' : "tau #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{vis} (GeV)",
	'ipsig_1': "electron IP sig",
        'ipsig_2': "tau IP sig",
        'n_jets': "number of jets",
        'n_bjets': "number of b-jets",
        'jpt_1': "leading jet p_{T} (GeV)",
        'jpt_2': "trailing jet p_{T} (GeV)",
        'mjj': "dijet mass (GeV)",
        'jdeta': "#Delta#eta(j,j)"
    }
}

def QCDEstimate(hists,perbin):
    histIsoSS = hists['iso_SS']
    histQCD = histIsoSS.Clone('QCD') 
    nbins = histQCD.GetNbinsX()
    histInvIsoOS = hists['inviso_OS']
    histInvIsoSS = hists['inviso_SS']
    normOS = histInvIsoOS.GetSumOfWeights()
    normSS = histInvIsoSS.GetSumOfWeights()
    ratio = 1.0
    #    ratio = normOS/normSS
    #    print('ratio = %5.3f'%(ratio))


    for ib in range(1,nbins+1):
        xOS = histInvIsoOS.GetBinContent(ib)
        eOS = histInvIsoOS.GetBinError(ib)
        xSS = histInvIsoSS.GetBinContent(ib)
        eSS = histInvIsoSS.GetBinError(ib)
        xData = histIsoSS.GetBinContent(ib)
        eData = histIsoSS.GetBinError(ib)
        xQCD = ratio * xData
        eQCD = ratio * eData 
        if perbin:
            if xOS>0 and xSS>0:
                xratio = xOS/xSS
                if xratio>0.2 and xratio<5.0:
                    xQCD = xratio*xData
                    rOS = eOS/xOS
                    rSS = eSS/xSS
                    rData = eData/xData
                    rQCD = math.sqrt(rOS*rOS+rSS*rSS+rData*rData)
                    
        if xQCD<0:
            xQCD=0
            eQCD=0
        histQCD.SetBinContent(ib,xQCD)
        histQCD.SetBinError(ib,eQCD)
        
    return histQCD

def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3_2022EE')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    suffix = kwargs.get('suffix','')
    per_process = kwargs.get('per_process',True)
    plotLegend = kwargs.get('plotLegend',True)
    ymin = kwargs.get('ymin',0.701)
    ymax = kwargs.get('ymax',1.299)

    # histograms
    h_data = hists['data_'+var+'_os_iso_all'].Clone('h_data')
    h_ztt = hists['dy_'+var+'_os_iso_tau'].Clone('h_ztt')
    h_zll = hists['dy_'+var+'_os_iso_lep'].Clone('h_zll')
    h_zll.Add(h_zll,hists['dy_'+var+'_os_iso_had'],1.,1.)
    h_top = hists['top_'+var+'_os_iso_all'].Clone('h_top')
    h_vv = hists['vv_'+var+'_os_iso_all'].Clone('h_vv')
    h_wjets = hists['wjets_'+var+'_os_iso_all'].Clone('h_wjets')

    h_lep = hists['dy_'+var+'_os_iso_lep'].Clone('h_lep')
    h_tau = hists['dy_'+var+'_os_iso_tau'].Clone('h_tau')
    h_had = hists['dy_'+var+'_os_iso_had'].Clone('h_had')

    for mc_sample in ['top','vv']:
        h_lep.Add(h_lep,hists[mc_sample+'_'+var+'_os_iso_lep'],1.,1.)
        h_tau.Add(h_tau,hists[mc_sample+'_'+var+'_os_iso_tau'],1.,1.)
        h_had.Add(h_had,hists[mc_sample+'_'+var+'_os_iso_had'],1.,1.)
    
    # qcd estimation
    hists_qcd = {}
    hists_qcd['iso_SS'] = hists['data_'+var+'_ss_iso_all'].Clone('iso_SS')
    hists_qcd['inviso_SS'] = hists['data_'+var+'_ss_antiiso_all'].Clone('inviso_SS')
    hists_qcd['inviso_OS'] = hists['data_'+var+'_os_antiiso_all'].Clone('inviso_OS')
    for mc_sample in ['dy','top','vv','wjets']:
        hists_qcd['iso_SS'].Add(hists_qcd['iso_SS'],hists[mc_sample+'_'+var+'_ss_iso_all'],1.,-1.)
        hists_qcd['inviso_SS'].Add(hists_qcd['inviso_SS'],hists[mc_sample+'_'+var+'_ss_antiiso_all'],1.,-1.)
        hists_qcd['inviso_OS'].Add(hists_qcd['inviso_OS'],hists[mc_sample+'_'+var+'_os_antiiso_all'],1.,-1.)

    h_qcd = QCDEstimate(hists_qcd,False)

    styles.InitData(h_data)

    xtitle = XTitle[chan][var]
    styles.InitHist(h_ztt,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor("#9999CC"),1001)
    styles.InitHist(h_vv,"","",ROOT.TColor.GetColor("#DE5A6A"),1001)
    styles.InitHist(h_qcd,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)

    styles.InitHist(h_tau,"","",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_lep,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_had,"","",ROOT.TColor.GetColor("#9999CC"),1001)
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
    x_had = h_had.GetSumOfWeights()
        
    print('')
    print('Yields ->')
    if per_process:
        print('Ztautau    : %7.0f'%(x_ztt))
        print('Zll        : %7.0f'%(x_zll))
        print('TTbar      : %7.0f'%(x_top))
        print('VV+ST      : %7.0f'%(x_vv))
        print('WJets      : %7.0f'%(x_wjets))
        print('QCD        : %7.0f'%(x_qcd))
    else:
        print('Prompt lep : %7.0f'%(x_lep))
        print('tau -> lep : %7.0f'%(x_tau))
        print('had -> lep : %7.0f'%(x_had))
        print('jet -> tau : %7.0f'%(x_wjets))
        print('QCD        : %7.0f'%(x_qcd))
        x_tot = x_tau + x_lep + x_had + x_wjets + x_qcd
        
    print('Total MC   : %7.0f'%(x_tot))
    print('Data       : %7.0f'%(x_data))
    print('')

    h_vv.Add(h_vv,h_wjets,1.,1.)
    h_vv.Add(h_vv,h_qcd,1.,1.)
    h_top.Add(h_top,h_vv,1.,1.)
    h_zll.Add(h_zll,h_top,1.,1.)
    h_ztt.Add(h_ztt,h_zll,1.,1.)

    h_wjets.Add(h_wjets,h_qcd,1.,1.)
    h_had.Add(h_had,h_wjets,1.,1.)
    h_lep.Add(h_lep,h_had,1.,1.)
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
    utils.zeroBinErrors(h_had)
    
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
        h_had.Draw('hsame')
        h_wjets.Draw('hsame')
        
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
            leg.AddEntry(h_had,'jet#rightarrow#mu','f')
        else:
            leg.AddEntry(h_tau,'#tau#rightarrowe','f')
            leg.AddEntry(h_lep,'prompt e','f')
            leg.AddEntry(h_had,'jet#rightarrowe','f')
        leg.AddEntry(h_wjets,'jet#rightarrow#tau','f')
            
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
    typ_suffix = 'typ'
    if per_process: typ_suffix = 'proc' 
    outputGraphics = os.getenv('CMSSW_BASE') + '/src/CPHiggs/IP/figures/' + var + '_' + chan + '_' + era + '_' + typ_suffix + suffix + '.png'    
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix','Run3_2022All','Run3_2023All'])
    parser.add_argument('-variable' ,'--variable', dest='variable', default='m_vis')
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-perType','--perType', dest='perType', action='store_true')
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',type=int,default=0)
    parser.add_argument('-mtCut','--mtCut',dest='mtCut',type=int,default=0)
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=48)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=0.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=240.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    parser.add_argument('-applyIP','--applyIP',dest='applyIP',type=int,default=0)
    parser.add_argument('-applySF', '--applySF', dest='applySF' ,type=int,default=0)
    
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
    applyIP = args.applyIP
    applySF = args.applySF
    
    plotLegend = True
    if var=='eta_1' or var=='eta_2':
        plotLegend = False
    
    proc = True
    if args.perType: proc = False
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    suffix_mt = ''
    suffix_xtrig = ''
    suffix_ip = ''
    suffix_sf = ''
    if applyMTCut==1:
        suffix_mt = '_xtrig_mtcut'
    if useCrossTrigger==1:
        suffix_xtrig = '_xtrig'        
    if applyIP==1:
        suffix_ip = '_ipcut1'
    if applySF==1:
        suffix_sf = '_promptSF_tauSF'

    suffix = '%s%s%s%s'%(suffix_mt,suffix_xtrig,suffix_ip,suffix_sf)
    basename = '%s/src/CPHiggs/IP'%(os.getenv('CMSSW_BASE'))
    inputFileName = '%s/selection/%s_%s%s.root'%(basename,chan,era,suffix)
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
    hists = utils.extractHistos(inputFile,var,bins)
    Plot(hists,era=era,var=var,channel=chan,per_process=proc,suffix=suffix,ymin=ymin,ymax=ymax,plotLegend=plotLegend)
    
