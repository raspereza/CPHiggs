#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Tag-and-Probe measurement of SF with Z->ee(mumu) samples
#
import ROOT
import math
from array import array
import os

import CPHiggs.IP.styles as styles
import CPHiggs.IP.utils as utils

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
    
    header = 'prompt #mu, |#eta|<1.2'
    if chan=='mm':
        if binEta=='1':
            header = 'prompt #mu, |#eta|<1.2'
        else:
            header = 'prompt #mu, |#eta|>1.2'
    else:
        if binEta=='1':
            header = 'prompt e, |#eta|<1.48'
        else:
            header = 'prompt e, |#eta|>1.48'

    
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
    sf.GetYaxis().SetTitle('SF')
    sf.GetXaxis().SetTitle('p_{T} (GeV)')
    sf.GetYaxis().SetRangeUser(ymin,ymax)
    styles.InitRatioHist(sf)
    
    # canvas and pads
    canvas = styles.MakeCanvas("canv","",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.31,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    

    eff_data.Draw('e')
    eff_mc.Draw('esame')

    leg = ROOT.TLegend(0.25,0.55,0.45,0.8)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.05)
    leg.SetHeader(header)
    leg.AddEntry(eff_data,'data','lp')
    leg.AddEntry(eff_mc,'simulation','lp')
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

    sf.Draw('e1')

    xmin = sf.GetXaxis().GetBinLowEdge(1)    
    xmax = sf.GetXaxis().GetBinLowEdge(nbins+1)
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
    outputGraphics = '%s/src/CPHiggs/IP/figures/SF_%s_binEta%s_%s.png'%(os.getenv('CMSSW_BASE'),chan,binEta,era)    
    canvas.Print(outputGraphics)

    
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

    era = kwargs.get('era','Run3_2022EE')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    xmin = kwargs.get('xmin',60.)
    xmax = kwargs.get('xmax',130.)
    binPt = kwargs.get('binPt','1')
    binEta = kwargs.get('binEta','1')
    sample = kwargs.get('sample','data')
    region = kwargs.get('region','pass')
    if sample not in ['data','dy']:
        print('Unknow sample %s : quitting'%(sample))
        exit()
    
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
    canvas = styles.MakeCanvas("canv","",700,600)
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
    outputGraphics = os.getenv('CMSSW_BASE') + '/src/CPHiggs/IP/figures/' + sample + '_' + chan + '_ptBin' + binPt + '_etaBin' + binEta + '_' + era + '.png'    
    canvas.Print(outputGraphics)

    del canvas
    return nCount,nIntegral,eIntegral
    
if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mm',choices=['mm','ee'])
    parser.add_argument('-nbins','--nbins', dest='nbins', type=int, default=60)
    parser.add_argument('-xmin','--xmin', dest='xmin', type=float, default=60.0)
    parser.add_argument('-xmax','--xmax', dest='xmax', type=float, default=120.)
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    nbins = args.nbins
    xmin = args.xmin
    xmax = args.xmax
    ymin = args.ymin
    ymax = args.ymax

    basedir = '%s/src/CPHiggs/IP'%(os.getenv('CMSSW_BASE'))
    
    bins = []
    width = (xmax-xmin)/float(nbins)
    for i in range(0,nbins+1):
        xb = xmin + width*float(i)
        bins.append(xb)

    inputFileName = '%s/selection/%s_%s.root'%(basedir,chan,era)
    inputFile = ROOT.TFile(inputFileName,'read')
    print('')
    print(inputFile)
    print('')
    hists = utils.extractTagProbeHistos(inputFile,bins)
    histPtBins = inputFile.Get('ptBins')
    histEtaBins = inputFile.Get('etaBins')
    nbinsPt = histPtBins.GetNbinsX()
    nbinsEta = histEtaBins.GetNbinsX()

    ptBins = []
    for iPt in range(1,nbinsPt+1):
        ptBins.append(histPtBins.GetBinLowEdge(iPt))
    ptBins.append(100.)

    etaBins = []
    for iEta in range(1,nbinsEta+2):
        etaBins.append(histEtaBins.GetBinLowEdge(iEta))

    print('ptBins',ptBins)
    print('etaBins',etaBins)
        
    effData_2D = ROOT.TH2D('effData','',nbinsPt,array('d',list(ptBins)),nbinsEta,array('d',list(etaBins)))
    effMC_2D = ROOT.TH2D('effMC','',nbinsPt,array('d',list(ptBins)),nbinsEta,array('d',list(etaBins)))

    for iEta in range(1,nbinsEta+1):
        binEta = '%1i'%(iEta)
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
                                                        region='pass')
            nFailData,xFailData,eFailData = RunTagProbe(hists,
                                                        era=era,
                                                        channel=chan,
                                                        xmin=xmin,
                                                        xmax=xmax,
                                                        binPt=binPt,
                                                        binEta=binEta,
                                                        sample='data',
                                                        region='fail')
            
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
                                                  region='pass')
            nFailMC,xFailMC,eFailMC = RunTagProbe(hists,
                                                  era=era,
                                                  channel=chan,
                                                  xmin=xmin,
                                                  xmax=xmax,
                                                  binPt=binPt,
                                                  binEta=binEta,
                                                  sample='dy',
                                                  region='fail')
            effMC,effErrMC = ComputeEff(nPassMC,xPassMC,ePassMC,
                                        nFailMC,xFailMC,eFailMC)
            effMC_1D.SetBinContent(iPt,effMC)
            effMC_1D.SetBinError(iPt,effErrMC)
            effMC_2D.SetBinContent(iPt,iEta,effMC)
            effMC_2D.SetBinError(iPt,iEta,effErrMC)
        PlotTagProbe(effData_1D,effMC_1D,era=era,channel=chan,binEta=binEta,ymin=ymin,ymax=ymax)

    outputFileName = '%s/ScaleFactors/SF_%s_%s.root'%(basedir,chan,era)
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    effData_2D.Write('effData')
    effMC_2D.Write('effMC')
    outputFile.Close()
