#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Creation of datacards
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.utils as utils

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def Symmetrize(hist):
    nbins = hist.GetNbinsX()
    nbins_half = int((nbins+0.1)/2.0)
    for ib in range(1,nbins_half+1):
        x1 = hist.GetBinContent(ib)
        x2 = hist.GetBinContent(nbins+1-ib)
        x = 0.5*(x1+x2)
        hist.SetBinContent(ib,x)
        hist.SetBinContent(nbins+1-ib,x)

def Flatten(hist):
    nbins = hist.GetNbinsX()
    tot = hist.GetSumOfWeights()
    x = tot/float(nbins)
    for ib in range(1,nbins+1):
        hist.SetBinContent(ib,x)

def Plot(hists,**kwargs):

    era = kwargs.get('era','Run3')
    chan = kwargs.get('channel','mt')
    cat = kwargs.get('cat','higgs')
    suffix = kwargs.get('suffix','')
    plotLegend = kwargs.get('plotLegend',True)
    blindData = kwargs.get('blindData',False)
    smooth = kwargs.get('smooth',False)
    xtitle = kwargs.get('xtitle','BDT_{sig}')
    ymin = kwargs.get('ymin',0.401)
    ymax = kwargs.get('ymax',1.599)

    # histograms
    h_data = hists['data_obs'].Clone('h_data')
    h_ztt = hists['ZTT'].Clone('h_ztt')
    h_zll = hists['ZL'].Clone('h_zll')
    h_top = hists['TTT'].Clone('h_top')
    h_vv = hists['VVT'].Clone('h_vv')
    h_qcd = hists['JetFakes'].Clone('h_qcd')
    h_sm = hists['ggH_sm_prod_sm_htt125'].Clone('h_sm')
    h_ps = hists['ggH_ps_prod_sm_htt125'].Clone('h_ps')
    h_sm.Add(h_sm,hists['qqH_sm_htt125'],1.,1.)
    h_ps.Add(h_ps,hists['qqH_ps_htt125'],1.,1.)

    h_ps.Scale(h_sm.GetSumOfWeights()/h_ps.GetSumOfWeights())
    
    styles.InitData(h_data)

    if smooth:     
        Flatten(h_ztt)
        Flatten(h_zll)
        Symmetrize(h_top)
        Symmetrize(h_vv)
        Symmetrize(h_qcd)
        Symmetrize(h_sm)
        Symmetrize(h_ps)
    
    styles.InitHist(h_ztt,xtitle,"",ROOT.TColor.GetColor("#FFCC66"),1001)
    styles.InitHist(h_zll,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_top,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_vv,"","",ROOT.TColor.GetColor(100,192,232),1001)
    styles.InitHist(h_qcd,"","",ROOT.TColor.GetColor("#FFCCFF"),1001)
    
    styles.InitModel(h_sm,"","",2)
    styles.InitModel(h_ps,"","",4)

    h_sm_bkg = h_sm.Clone('h_sm_ratio')
    h_ps_bkg = h_ps.Clone('h_ps_ratio')
    
    x_data = h_data.GetSumOfWeights()
    
    x_ztt   = h_ztt.GetSumOfWeights()
    x_zll   = h_zll.GetSumOfWeights() 
    x_top   = h_top.GetSumOfWeights()
    x_vv    = h_vv.GetSumOfWeights() 
    x_qcd   = h_qcd.GetSumOfWeights()
    
    x_tot   = x_ztt + x_zll + x_top + x_vv + x_qcd
    x_sm = h_sm.GetSumOfWeights()
    x_ps = h_ps.GetSumOfWeights()
    
    print('')
    print('Category %s'%(cat))
    print('Ztautau    : %7.0f'%(x_ztt))
    print('Zll        : %7.0f'%(x_zll))
    print('TTbar      : %7.0f'%(x_top))
    print('VV+ST      : %7.0f'%(x_vv))
    print('Fakes      : %7.0f'%(x_qcd))
    print('')
    print('Total MC   : %7.0f'%(x_tot))
    print('Data       : %7.0f'%(x_data))
    print('')
    print('CP_even    : %7.1f'%(x_sm))
    print('CP_odd     : %7.1f'%(x_ps))
    
    h_zll.Add(h_zll,h_vv,1.,1.)
    h_top.Add(h_top,h_zll,1.,1.)
    h_qcd.Add(h_qcd,h_top,1.,1.)    
    h_ztt.Add(h_ztt,h_qcd,1.,1.)

    h_tot = h_ztt.Clone("total")

    if blindData:
        nbins = h_data.GetNbinsX()
        for ib in range(1,nbins+1):
            xbkg = h_tot.GetBinContent(ib)
            xsm = h_sm.GetBinContent(ib)
            xps = h_ps.GetBinContent(ib)
            xsig = 0.5*(xsm+xps)
            xratio = xsig/xbkg
            if xratio>0.05:
                h_data.SetBinContent(ib,10000000.)
                h_data.SetBinError(ib,0.)
                
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')
    h_sm_bkg.Add(h_sm_bkg,h_tot,1.,1.) 
    h_ps_bkg.Add(h_ps_bkg,h_tot,1.,1.) 
    h_sm_ratio = utils.histoRatio(h_sm_bkg,h_tot,'ratio')
    h_ps_ratio = utils.histoRatio(h_ps_bkg,h_tot,'ratio')
    
    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    styles.InitRatioHist(h_tot_ratio)
    h_tot_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    utils.zeroBinErrors(h_ztt)
    utils.zeroBinErrors(h_zll)
    utils.zeroBinErrors(h_top)
    utils.zeroBinErrors(h_vv)
    utils.zeroBinErrors(h_qcd)

    utils.zeroBinErrors(h_sm)
    utils.zeroBinErrors(h_ps)
    utils.zeroBinErrors(h_sm_bkg)
    utils.zeroBinErrors(h_ps_bkg)
    
    YMax = h_tot.GetMaximum()
    
    h_ztt.GetYaxis().SetRangeUser(0.5,5*YMax)
    h_ztt.GetXaxis().SetLabelSize(0)
    h_ztt.GetYaxis().SetTitle("events")
    h_tot_ratio.GetYaxis().SetTitle("obs/exp")
    h_tot_ratio.GetXaxis().SetTitle(xtitle)

    # canvas and pads
    canvas = styles.MakeCanvas('canv_'+cat,'',600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    
    
    h_ztt.Draw('h')
    h_qcd.Draw('hsame')
    h_top.Draw('hsame')
    h_data.Draw('e1same')
    h_tot.Draw('e2same')
    h_ps.Draw('hsame')
    h_sm.Draw('hsame')
    
    leg1 = ROOT.TLegend(0.20,0.78,0.45,0.88)
    styles.SetLegendStyle(leg1)
    leg1.SetTextSize(0.036)
    leg1.AddEntry(h_data,'data','lp')
    leg1.AddEntry(h_ztt,'Z#rightarrow#tau#tau','f')

    leg2 = ROOT.TLegend(0.35,0.78,0.60,0.88)
    styles.SetLegendStyle(leg2)
    leg2.SetTextSize(0.036)
    leg2.AddEntry(h_qcd,'fakes','f')
    leg2.AddEntry(h_top,'other','f')

    leg3 = ROOT.TLegend(0.50,0.78,0.75,0.88)
    styles.SetLegendStyle(leg3)
    leg3.SetTextSize(0.036)
    leg3.AddEntry(h_sm,'H#rightarrow#tau#tau (CP-even)','l')
    leg3.AddEntry(h_ps,'H#rightarrow#tau#tau (CP-odd)','l')

    if plotLegend:
        leg1.Draw()
        leg2.Draw()
        leg3.Draw()

    styles.CMS_label(upper,era=era)

    upper.SetLogy(True)
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

    h_tot_ratio.Draw('e2')
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
    outputGraphics = '%s/figures/datacards/%s/%s/cards_%s_%s'%(utils.outputFolder,chan,suffix,era,cat)
    if blindData:
        outputGraphics += '_blind'
    outputGraphics += '.png'
    canvas.Print(outputGraphics)

    yields = {}
    yields['ztt'] = x_ztt
    yields['zll'] = x_zll
    yields['top'] = x_top
    yields['vv'] = x_vv
    yields['fakes'] = x_qcd
    yields['data'] = x_data
    yields['higgs_sm'] = x_sm
    yields['higgs_ps'] = x_ps
    return yields


if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3',choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-bVeto','--bVeto', dest='bVeto',action='store_true')
    parser.add_argument('-useCrossTrigger','--useCrossTrigger', dest='useCrossTrigger',action='store_true')
    parser.add_argument('-blindData','--blindData', dest='blindData', action='store_true')
    parser.add_argument('-cf','--cf', dest='cf', action='store_true')
    parser.add_argument('-ymin','--ymin',dest='ymin',type=float,default=0.401)
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.599)
    args = parser.parse_args()

    procs = ['data_obs','ZTT','ZL','TTT','VVT','JetFakes','ggH_ps_prod_sm_htt125','ggH_sm_prod_sm_htt125','qqH_sm_htt125','qqH_ps_htt125']

    plotLegend = True 
    era = args.era
    bVeto = args.bVeto
    useCrossTrigger = args.useCrossTrigger
    chan = args.channel
    ymin = args.ymin
    ymax = args.ymax
    blindData = args.blindData
    cats_bdt = ['higgs','tau','fake']
    modes = ['pi','rho','a1','a11pr']
    chan_lep = {
        'mt' : 'mu',
        'et' : 'e',
    }
    nbdt_bins = 5
    cf = args.cf
    suffix_out = 'x'
    if bVeto:
        suffix_out += '_bveto'
    if useCrossTrigger:
        suffix_out += '_xtrig'
    
    suffix = ''
    if bVeto:
        suffix += '_bveto'     
    if useCrossTrigger:
        suffix += '_xtrig'

    fileName = utils.outputFolder+'/datacardsPhiCP/'+chan+'_'+era+suffix+'.root'
    if cf:
        suffix_out = 'CF'
        nbdt_bins = 3
        cats_bdt.remove('higgs')
        fileName = utils.outputFolder+'/datacardsPhiCP/'+chan+'_2022_2023.root'

    if os.path.isfile(fileName):
        print('')
        print('Loading ROOT file %s'%(fileName))
    else:
        print('')
        print('ROOT file %s is not found'%(fileName))
        print('Quitting')
        exit()

    inputfile = ROOT.TFile(fileName,'READ')

    ztt = {}
    zll = {}
    top = {}
    vv  = {}
    fakes = {}
    data = {}
    higgs_sm = {}
    higgs_ps = {}
    for cat in cats_bdt:
        cat_name = '%s_mva_%s'%(chan,cat)
        hists = {}
        for proc in procs:
            hists[proc] = inputfile.Get(cat_name+'/'+proc)
#            print(proc,hists[proc])

        xtitle=''
        smooth = False
        if cat=='higgs':
            xtitle='BDT_{higgs}'
        elif cat=='tau':
            xtitle='BDT_{ditau}'
        else:
            xtitle='BDT_{fakes}'
        yields = Plot(hists,era=era,channel=chan,cat=cat,suffix=suffix_out,
                      plotLegend=plotLegend,blindData=blindData,xtitle=xtitle,
                      smooth=smooth,ymin=ymin,ymax=ymax)
        data[cat] = yields['data'] 
        ztt[cat] = yields['ztt'] 
        zll[cat] = yields['zll'] 
        top[cat] = yields['top'] 
        vv[cat] = yields['vv'] 
        fakes[cat] = yields['fakes'] 
        higgs_sm[cat] = yields['higgs_sm'] 
        higgs_ps[cat] = yields['higgs_ps'] 

    for dm in modes:
        n_data = 0.
        n_ztt  = 0.
        n_zll  = 0.
        n_vv = 0.
        n_top = 0.
        n_fakes = 0.
        n_higgs_sm = 0.
        n_higgs_ps = 0.
        for ib in range(0,nbdt_bins):
            bdt_bin = '%1i'%(ib)
            cat_name = '%s_mva_higgs_cat%s_%s%s'%(chan,bdt_bin,chan_lep[chan],dm)
            hists = {}
            for proc in procs:
                hists[proc] = inputfile.Get(cat_name+'/'+proc)
            cat = 'cat%s_%s%s'%(bdt_bin,chan_lep[chan],dm)
            xtitle='#phi_{CP} [deg]'
            smooth = True
            yields = Plot(hists,era=era,channel=chan,cat=cat,suffix=suffix_out,
                          plotLegend=plotLegend,blindData=blindData,xtitle=xtitle,
                          smooth=smooth,ymin=ymin,ymax=ymax)
            n_data += yields['data']
            n_ztt += yields['ztt']
            n_zll += yields['zll']
            n_top += yields['top']
            n_vv += yields['vv']
            n_fakes += yields['fakes']
            n_higgs_sm += yields['higgs_sm']
            n_higgs_ps += yields['higgs_ps']
            
        data[dm] = n_data
        ztt[dm] = n_ztt
        zll[dm] = n_zll
        top[dm] = n_top
        vv[dm] = n_vv
        fakes[dm] = n_fakes
        higgs_sm[dm] = n_higgs_sm
        higgs_ps[dm] = n_higgs_ps
 
    print('')
    frmt = {
        'tau'   : 'tau  ',
        'fake'  : 'fake ',
        'pi'    : 'pi   ',
        'rho'   : 'rho  ',
        'a1'    : 'a1   ',
        'a11pr' : 'a11pr',
    }
    for cat in ['tau','fake','pi','rho','a1','a11pr']:
        other = zll[cat]+top[cat]+vv[cat]
        total = zll[cat]+top[cat]+vv[cat]+ztt[cat]+fakes[cat]
        print('%s  %6.0f %6.0f %6.0f %6.0f %3.0f'%(frmt[cat],data[cat],total,ztt[cat],fakes[cat],higgs_sm[cat]))
