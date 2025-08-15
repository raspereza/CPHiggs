#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Tag-and-Probe measurement of SF with
# Z->ee(mumu) samples
# ########################################
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def FitFuncSF(x,par):
    b = par[0]+par[1]*x[0]+par[2]*x[0]*x[0]
    return b

def FitFunc(x,par):
    a = (x[0]-par[1])/par[2]
    b = (x[0]-par[1])/par[3]
    if x[0]>par[1]:
        a = (x[0]-par[1])/par[4]
        b = (x[0]-par[1])/par[5]
    signal = par[6]*ROOT.TMath.Exp(-0.5*a*a)+(1-par[6])*ROOT.TMath.Exp(-0.5*b*b)
    signal = par[0]*signal
    bkg = par[7]*ROOT.TMath.Exp(-x[0]*par[8])
    return signal + bkg

def FitFuncAlt(x,par):
    a = (x[0]-par[1])/par[2]
    b = (x[0]-par[1])/par[3]
    if x[0]>par[1]:
        a = (x[0]-par[1])/par[4]
        b = (x[0]-par[1])/par[5]
    signal = par[6]*ROOT.TMath.Exp(-0.5*a*a)+(1-par[6])*ROOT.TMath.Exp(-0.5*b*b)
    signal = par[0]*signal
    bkg = par[7]+x[0]*par[8]
    return signal + bkg

def SignalFunc(x,par):
    a = (x[0]-par[1])/par[2]
    b = (x[0]-par[1])/par[3]
    if x[0]>par[1]:
        a = (x[0]-par[1])/par[4]
        b = (x[0]-par[1])/par[5]
    signal = par[6]*ROOT.TMath.Exp(-0.5*a*a)+(1-par[6])*ROOT.TMath.Exp(-0.5*b*b)
    signal = par[0]*signal
    return signal

# computing efficiency
def ComputeEff(nCountPass,xIntegralPass,eIntegralPass,
               nCountFail,xIntegralFail,eIntegralFail):
    sysPass = abs(xIntegralPass-nCountPass)
    sysFail = abs(xIntegralFail-nCountFail)
    errPass = math.sqrt(sysPass*sysPass+eIntegralPass*eIntegralPass)
    errFail = math.sqrt(sysFail*sysFail+eIntegralFail*eIntegralFail)
    if errPass>0.02*xIntegralPass: errPass=0.02*xIntegralPass
    if errFail>0.02*xIntegralFail: errFail=0.02*xIntegralFail
    xTot = xIntegralFail+xIntegralPass    
    eff = xIntegralPass/xTot
    relFail = xIntegralPass*errFail/(xTot*xTot)
    relPass = xIntegralFail*errPass/(xTot*xTot)
    effErr = math.sqrt(relFail*relFail+relPass*relPass)
    return eff,effErr
    

def PlotTagProbe(eff_data,eff_mc,**kwargs):
    era = kwargs.get('era','Run3_2022EE')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    binEta = kwargs.get('binEta','1')
    ymin = kwargs.get('ymin',0.5)
    ymax = kwargs.get('ymax',1.5)
    secondLep = kwargs.get('secondLep',False)
    generator = kwargs.get('generator','amcatnlo')
    
    header = 'prompt #mu, |#eta|<1.0'
    if chan=='mm':
        if binEta=='1':
            header = 'prompt #mu, |#eta|<1.0'
        elif binEta=='2':
            header = 'prompt #mu, 1.0<|#eta|<1.6'
        else:
            header = 'prompt #mu, 1.6<|#eta|<2.4'
    else:
        if binEta=='1':
            header = 'prompt e, |#eta|<1.0'
        elif binEta=='2':
            header = 'prompt e, 1.0<|#eta|<1.6'
        else:
            header = 'prompt e, 1.6<|#eta|<2.1'

    
    styles.InitData(eff_data)
    styles.InitData(eff_mc)
    eff_data.SetMarkerStyle(20)
    eff_data.SetMarkerColor(ROOT.kRed)
    eff_data.SetLineColor(ROOT.kRed)
    eff_mc.SetMarkerStyle(21)
    eff_mc.SetMarkerColor(ROOT.kBlue)
    eff_mc.SetLineColor(ROOT.kBlue)

    sf = eff_data.Clone('sf')
    nbins = sf.GetNbinsX()
    fmin = sf.GetXaxis().GetBinLowEdge(1)
    fmax = sf.GetXaxis().GetBinLowEdge(nbins+1)

    print('fmin=%2.0f  fmax=%3.0f'%(fmin,fmax))
    
    for ib in range(1,nbins+1):
        xdata = eff_data.GetBinContent(ib)
        edata = eff_data.GetBinError(ib)
        rdata = edata/xdata
        xmc = eff_mc.GetBinContent(ib)
        emc = eff_mc.GetBinError(ib)
        rmc = emc/xmc
        xsf = xdata/xmc
        rsf = math.sqrt(rdata*rdata+rmc*rmc)
        esf = xsf*rsf
        sf.SetBinContent(ib,xsf)
        sf.SetBinError(ib,esf)

    eff_data.GetYaxis().SetRangeUser(0.,1.1)
    eff_data.GetXaxis().SetLabelSize(0)
    eff_data.GetYaxis().SetTitle('Efficiency')
    eff_data.GetXaxis().SetRangeUser(fmin,fmax)
    
    sf.GetYaxis().SetTitle('SF')
    sf.GetXaxis().SetTitle('p_{T} (GeV)')
    sf.GetYaxis().SetRangeUser(ymin,ymax)
    sf.GetXaxis().SetRangeUser(fmin,fmax)
    styles.InitRatioHist(sf)

    print 
    fitSF = ROOT.TF1('fitFuncSF',FitFuncSF,fmin,fmax,3)
    fitSF.SetLineColor(4)
    fitSF.SetLineWidth(2)
    fitSF.SetParName(0,'a0')
    fitSF.SetParName(1,'a1')
    fitSF.SetParName(2,'a2')
    fitSF.SetParameter(0,1.)
    fitSF.SetParameter(1,0.)
    fitSF.SetParameter(2,0.)
    
    xtit = ''
    ytit = 'SF (IP sig. cut)'
    if chan=='mm':
        xtit = 'muon p_{T} (GeV)'
    else:
        ytit = 'electron p_{T} (GeV)'
    postfix = 'binEta%s'%(binEta)
    if secondLep:
        postfix += '_2'
    dummy = styles.MakeCanvas('dummy_'+postfix,'',400,400)
    sf.Fit(fitSF,"R")
    hfit = ROOT.TH1D("hfit_"+postfix,"",200,fmin,fmax)
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hfit,0.68)
    styles.InitModel(hfit,xtit,ytit,4)
    hfit.SetFillColor(ROOT.kCyan)
    hfit.SetFillStyle(1001)
    hfit.SetLineWidth(0)
    hfit.SetLineColor(4)
    hfit.SetMarkerSize(0)
    hfit.SetMarkerStyle(0)
    if chan=='mm':
        hfit.GetXaxis().SetTitle("muon p_{T} [GeV]")
    else:
        hfit.GetXaxis().SetTitle("electron p_{T} [GeV]")
    hfit.GetYaxis().SetTitle("scale factor (IP sig.)")

    
    # canvas and pads
    canvas = styles.MakeCanvas("canv_"+postfix,"",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    

    eff_data.Draw('e')
    eff_mc.Draw('esame')

    leg = ROOT.TLegend(0.25,0.5,0.5,0.75)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.05)
    leg.SetHeader(header)
    leg.AddEntry(eff_data,'data','lp')
    leg.AddEntry(eff_mc,'simulation','lp')
    if plotLegend: leg.Draw()

    styles.CMS_label(upper,era=era)

    upper.SetLogx(True)
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

    sf.Draw('e1')
    hfit.Draw('e2same')
    sf.Draw('e1same')

    xmin = sf.GetXaxis().GetBinLowEdge(1)    
    xmax = sf.GetXaxis().GetBinLowEdge(nbins+1)
    line = ROOT.TLine(xmin,1.,xmax,1.)
    line.SetLineStyle(1)
    line.SetLineWidth(2)
    line.SetLineColor(4)
    line.Draw()
    lower.SetLogx(True)
    lower.Modified()
    lower.RedrawAxis()

    canvas.cd()
    canvas.Modified()
    canvas.cd()
    canvas.SetSelected(canvas)
    canvas.Update()

    print('')
    suffix = generator
    if secondLep:
        suffix += '_2'
    outputGraphics = '%s/figures/SF_%s_binEta%s_%s_%s.png'%(utils.outputFolder,chan,binEta,era,suffix)    
    canvas.Print(outputGraphics)
    return hfit
    
def ComputeYield(hist):

    nLow = hist.FindBin(61.)
    nSignalLow  = hist.FindBin(71.)
    nSignalHigh = hist.FindBin(109.)
    nHigh = hist.FindBin(119.)

    xBkgLeft = hist.Integral(nLow,nSignalLow-1)
    xBkgRight = hist.Integral(nSignalHigh+1,nHigh)
    xBkg = 2.0*(xBkgLeft+xBkgRight)
    xPeak = hist.Integral(nSignalLow,nSignalHigh)
    xSignal = xPeak - xBkg
    return xSignal
    

def RunTagProbe(hists,**kwargs):

    era = kwargs.get('era','Run3_2022')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    xmin = kwargs.get('xmin',60.)
    xmax = kwargs.get('xmax',120.)
    binPt = kwargs.get('binPt','1')
    binEta = kwargs.get('binEta','1')
    sample = kwargs.get('sample','data')
    region = kwargs.get('region','pass')
    secondLep = kwargs.get('secondLep',False)
    generator = kwargs.get('generator','amcatnlo')
    if sample not in ['data','dy']:
        print('Unknow sample %s : quitting'%(sample))
        exit()

    postfixLep = ''
    if secondLep:
        postfixLep = '_2'
    # histogram name
    name = '%s_m_vis_%s_%s_%s_os_iso_all'%(sample,region,binPt,binEta)
    hist = hists[name].Clone(name+'_fit')    
    
    styles.InitData(hist)
    hist.SetMarkerSize(1.0)

    xtitle = 'm_{#mu#mu} (GeV)'
    header = 'Z#rightarrow#mu#mu'
    if chan=='ee':
        xtitle = 'm_{ee} (GeV)'
        header = 'Z#rightarrowee'

    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle('Events')

    fitFunc = ROOT.TF1('fitFunc',FitFunc,xmin,xmax,9)
    fitFunc.SetLineColor(2)
    # initialize parameters
    # background exponent
    x1 = xmin + 1.
    x2 = 125.
    y1 = hist.GetBinContent(hist.FindBin(x1))
    y2 = hist.GetBinContent(hist.FindBin(x2))
    if y1<2: y1=2
    if y2<1: y2=1
    b = (ROOT.TMath.Log(y1)-ROOT.TMath.Log(y2))/(x2-x1)
    a = y1*ROOT.TMath.Exp(b*x1)
    # signal : two gaussians
    x0 = 90.
    sigmaLeft = 2.0
    sigmaRight = 2.0
    sigma2Left = 5.0
    sigma2Right = 4.0
    norm = hist.GetBinContent(hist.FindBin(x0)) # - a*ROOT.TMath.Exp(-b*x0)
    fraction = 0.5
    fitFunc.SetParameter(0,norm)
    fitFunc.SetParameter(1,x0)
    fitFunc.SetParameter(2,sigmaLeft)
    fitFunc.SetParameter(3,sigma2Left)
    fitFunc.SetParameter(4,sigmaRight)
    fitFunc.SetParameter(5,sigma2Right)
    fitFunc.SetParameter(6,fraction)
    fitFunc.SetParameter(7,a)
    fitFunc.SetParameter(8,b)

    # alternative func
    fitFuncAlt = ROOT.TF1('fitFuncAlt',FitFuncAlt,xmin,xmax,9)
    fitFuncAlt.SetLineColor(2)
    bAlt = (y2-y1)/(x2-x1)
    aAlt = y1 - bAlt*x1
    fitFuncAlt.SetParameter(0,norm)
    fitFuncAlt.SetParameter(1,x0)
    fitFuncAlt.SetParameter(2,sigmaLeft)
    fitFuncAlt.SetParameter(3,sigma2Left)
    fitFuncAlt.SetParameter(4,sigmaRight)
    fitFuncAlt.SetParameter(5,sigma2Right)
    fitFuncAlt.SetParameter(6,fraction)
    fitFuncAlt.SetParameter(7,aAlt)
    fitFuncAlt.SetParameter(8,bAlt)
    
    # canvas and pads
    postfix = '%s_%s_%s_%s'%(region,binPt,binEta,sample)
    canvas = styles.MakeCanvas("canv_"+postfix,"",700,600)
    sigFunc = ROOT.TF1('sigFunc',SignalFunc,xmin,xmax,7)
    hist.Fit('fitFuncAlt','LR')
    for npar in range(0,7):
        sigFunc.SetParameter(npar,fitFuncAlt.GetParameter(npar))
    nCount = sigFunc.Integral(70.,110.)
    hist.Fit('fitFunc','LR')
    for npar in range(0,7):
        sigFunc.SetParameter(npar,fitFunc.GetParameter(npar))
    nIntegral = sigFunc.Integral(70.,110.)
    eIntegral = nIntegral*fitFunc.GetParError(0)/fitFunc.GetParameter(0)

    hist.Draw('e1')    
    
    leg = ROOT.TLegend(0.65,0.55,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.042)
    leg.SetHeader(header)
    if sample=='data': leg.AddEntry(hist,'data','lp')
    else: leg.AddEntry(hist,'simulation','lp')
    if plotLegend: leg.Draw()

    styles.CMS_label(canvas,era=era)

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()
    print('')

    suffix = generator
    if secondLep:
        suffix += '_2'
    outputGraphics = utils.outputFolder + '/figures/tag_probe/' + sample + '_' + chan + '_ptBin' + binPt + '_etaBin' + binEta + '_' + era + '_' + region + '_' + suffix + '.png'
    canvas.Print(outputGraphics)

    del canvas
    return nCount,nIntegral,eIntegral
    
if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2023'])
    parser.add_argument('-channel','--channel', dest='channel',default='mm',choices=['mm','ee'])
    parser.add_argument('-secondLep','--secondLep', dest='secondLep',action='store_true')
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=60)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=60.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=120.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.501)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.499)
    parser.add_argument('-binEta','--binEta',dest='binEta',default='1',choices=['1','2','3'])
    parser.add_argument('-generator', '--generator',dest='generator',default='amcatnlo',choices=['amcatnlo','MG','powheg'])
    
    args = parser.parse_args()

    labels = {
        'mm': 'PromptMu',
        'ee': 'PromptE'
    }
    
    era = args.era
    chan = args.channel
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax
    generator = args.generator
    secondLep = args.secondLep
    binEta = args.binEta
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    inputFileName = '%s/selection/%s_%s_x.root'%(utils.outputFolder,chan,era)
    inputFile = ROOT.TFile(inputFileName,'read')
    print('')
    print(inputFile)
    print('')
    hists = utils.extractTagProbeHistos(inputFile,bins,generator,era,chan,secondLep)
    histPtBins = inputFile.Get('ptBins')
    histEtaBins = inputFile.Get('etaBins')
    nbinsPt = histPtBins.GetNbinsX()
    nbinsEta = histEtaBins.GetNbinsX()

    ptBins = []
    for iPt in range(1,nbinsPt+2):
        ptBins.append(histPtBins.GetBinLowEdge(iPt))
        
    etaBins = []
    for iEta in range(1,nbinsEta+2):
        etaBins.append(histEtaBins.GetBinLowEdge(iEta))

    print('nbinsPt',nbinsPt,'ptBins',ptBins)
    print('nbinsEta',nbinsEta,'etaBins',etaBins)
        
    effData_2D = ROOT.TH2D('effData','',nbinsPt,array('d',list(ptBins)),nbinsEta,array('d',list(etaBins)))
    effMC_2D = ROOT.TH2D('effMC','',nbinsPt,array('d',list(ptBins)),nbinsEta,array('d',list(etaBins)))

    hfit = {}
    effData_1D = ROOT.TH1D('effData_eta'+binEta,'',nbinsPt,array('d',list(ptBins)))
    effMC_1D = ROOT.TH1D('effMC_eta'+binEta,'',nbinsPt,array('d',list(ptBins)))
    for iPt in range(1,nbinsPt+1):
        binPt = '%1i'%(iPt)
        # tag-and-probe data
        nPassData,xPassData,ePassData = RunTagProbe(hists,
                                                    era=era,
                                                    channel=chan,
                                                    xmin=xmin,
                                                    xmax=xmax,
                                                    binPt=binPt,
                                                    binEta=binEta,
                                                    sample='data',
                                                    region='pass',
                                                    generator=generator,
                                                    secondLep=secondLep)
        nFailData,xFailData,eFailData = RunTagProbe(hists,
                                                    era=era,
                                                    channel=chan,
                                                    xmin=xmin,
                                                    xmax=xmax,
                                                    binPt=binPt,
                                                    binEta=binEta,
                                                    sample='data',
                                                    region='fail',
                                                    generator=generator,
                                                    secondLep=secondLep
                                                    )
        
        effData,effErrData = ComputeEff(nPassData,xPassData,ePassData,
                                        nFailData,xFailData,eFailData)
        print('pass : %1.0f (%1.0f) +/- %1.0f - fail : %1.0f (%1.0f) +/- %1.0f'%(nPassData,xPassData,ePassData,nFailData,xFailData,eFailData))
        effData_1D.SetBinContent(iPt,effData)
        effData_1D.SetBinError(iPt,effErrData)
        effData_2D.SetBinContent(iPt,iEta,effData)
        effData_2D.SetBinError(iPt,iEta,effErrData)

        nPassMC,xPassMC,ePassMC = RunTagProbe(hists,
                                              era=era,
                                              channel=chan,
                                              xmin=xmin,
                                              xmax=xmax,
                                              binPt=binPt,
                                              binEta=binEta,
                                              sample='dy',
                                              region='pass',
                                              generator=generator,
                                              secondLep=secondLep)
        nFailMC,xFailMC,eFailMC = RunTagProbe(hists,
                                              era=era,
                                              channel=chan,
                                              xmin=xmin,
                                              xmax=xmax,
                                              binPt=binPt,
                                              binEta=binEta,
                                              sample='dy',
                                              region='fail',
                                              generator=generator,
                                              secondLep=secondLep)
        effMC,effErrMC = ComputeEff(nPassMC,xPassMC,ePassMC,
                                    nFailMC,xFailMC,eFailMC)
        effMC_1D.SetBinContent(iPt,effMC)
        effMC_1D.SetBinError(iPt,effErrMC)
        effMC_2D.SetBinContent(iPt,iEta,effMC)
        effMC_2D.SetBinError(iPt,iEta,effErrMC)
    hfit[binEta] = PlotTagProbe(effData_1D,effMC_1D,era=era,channel=chan,binEta=binEta,ymin=ymin,ymax=ymax,generator=generator,secondLep=secondLep)

    suffix = generator
    if secondLep:
        suffix += '_2'
    outputFileName = '%s/ScaleFactors/SF_%s_%s_%s_%s.root'%(utils.outputFolder,labels[chan],binEta,era,suffix)
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    postfixLep = ''
    if secondLep:
        postfixLep = '_2'
    effData_2D.Write('effData'+postfixLep)
    effMC_2D.Write('effMC'+postfixLep)
    for h in hfit:
        hfit[h].Write('hfit_binEta'+h+postfixLep)
    outputFile.Close()
