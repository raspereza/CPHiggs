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

def ComputeEff(f):

    histPass = {}
    histFail = {}
    fits = ['shapes_prefit','shapes_fit_s']
    histNames = ['TTT','VVT']
    for fit in fits:
        histPass[fit] = f.Get('%s/ch1/ZTT_pass'%(fit))
        histFail[fit] = f.Get('%s/ch2/ZTT_fail'%(fit))
        print('histPass',histPass[fit])
        print('histFail',histFail[fit])
        
    for fit in fits:
        for name in histNames:
            histPass[fit].Add(histPass[fit],f.Get('%s/ch1/%s_pass'%(fit,name)),1.,1.)
            histFail[fit].Add(histFail[fit],f.Get('%s/ch2/%s_fail'%(fit,name)),1.,1.)

    passMC = histPass['shapes_prefit'].GetSumOfWeights()
    failMC = histFail['shapes_prefit'].GetSumOfWeights()
    passData = histPass['shapes_fit_s'].GetSumOfWeights()
    failData = histFail['shapes_fit_s'].GetSumOfWeights()

    r_pass = 1.0
    e_pass = 0.2
    r_fail = 1.0
    e_fail = 0.2
    fitResult = f.Get('fit_s')
    pars = fitResult.floatParsFinal()

    for par in pars:
        parname =  par.GetName()
        if parname=='r_pass':
            r_pass = par.getVal()
            e_pass = par.getError()
        if parname=='r_fail':
            r_fail = par.getVal()
            e_fail = par.getError()
    
    
    totMC = passMC+failMC
    effMC = passMC/totMC
    errMC = 0.

    totData = passData+failData
    effData = passData/totData
    errPass = passData*e_pass/r_pass
    errFail = failData*e_fail/r_fail

    relFail = passData*errFail/(totData*totData)
    relPass = failData*errPass/(totData*totData)
    errData = math.sqrt(relFail*relFail+relPass*relPass)
    
    return effData,errData,effMC,errMC
    

def PlotTagProbe(eff_data,eff_mc,**kwargs):
    era = kwargs.get('era','Run3_2022EE')
    chan = kwargs.get('channel','mm')
    plotLegend = kwargs.get('plotLegend',True)
    binEta = kwargs.get('binEta','1')
    ymin = kwargs.get('ymin',0.5)
    ymax = kwargs.get('ymax',1.5)
    
    header = '#tau#rightarrow#mu, |#eta|<1.2'
    if chan=='mt':
        if binEta=='1':
            header = '#tau#rightarrow#mu, |#eta|<1.0'
        elif binEta=='2':
            header = '#tau#rightarrow#mu, 1.0<|#eta|<1.6'
        else:
            header = '#tau#rightarrow#mu, 1.6<|#eta|<2.4'
    else:
        if binEta=='1':
            header = '#tau#rightarrowe, |#eta|<1.0'
        elif binEta=='2':
            header = '#tau#rightarrowe, 1.0<|#eta|<1.6'
        else:
            header = '#tau#rightarrowe, 1.6<|#eta|<2.1'

    
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

    leg = ROOT.TLegend(0.55,0.25,0.8,0.5)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.055)
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
    
if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix','Run3_2022All','Run3_2023All'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-ymin','--ymin', dest='ymin', type=float, default=0.701)
    parser.add_argument('-ymax','--ymax', dest='ymax', type=float, default=1.299)
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    ymin = args.ymin
    ymax = args.ymax
    
    basedir = '%s/src/CPHiggs/IP'%(os.getenv('CMSSW_BASE'))
    
    inputFileName = '%s/selection/%s_%s_xtrig_promptSF.root'%(basedir,chan,era)
    inputFile = ROOT.TFile(inputFileName,'read')
    print('')
    print(inputFile)
    print('')
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
    print('')
    
    effData = {}
    errData = {}
    effData['1'] = []
    errData['1'] = []
    effData['2'] = []
    errData['2'] = []
    effData['3'] = []
    errData['3'] = []
    effMC = {}
    errMC = {}
    effMC['1'] = []
    errMC['1'] = []
    effMC['2'] = []
    errMC['2'] = []
    effMC['3'] = []
    errMC['3'] = []
    
    cards_folder = '%s/src/CPHiggs/IP/datacards'%(os.getenv('CMSSW_BASE'))
    for iEta in range(1,nbinsEta+1):
        binEta = '%1i'%(iEta)
        effData_1D = ROOT.TH1D('effData_eta'+binEta,'',nbinsPt,array('d',list(ptBins)))
        effMC_1D = ROOT.TH1D('effMC_eta'+binEta,'',nbinsPt,array('d',list(ptBins)))
        for iPt in range(1,nbinsPt+1):
            binPt = '%1i'%(iPt)
            # tag-and-probe data
            inputFileName = '%s/%s_%s_binPt%s_binEta%s_fit.root'%(cards_folder,chan,era,binPt,binEta)
            inputFile = ROOT.TFile(inputFileName,'read')
            print('%s'%(inputFile))
            xData,eData,xMC,eMC = ComputeEff(inputFile)
            effData[binEta].append(xData)
            errData[binEta].append(eData)
            effMC[binEta].append(xMC)
            errMC[binEta].append(eMC)
            print('binPt%s_binEta%s -> Data : %4.2f+/-%4.2f : MC = %4.2f'%(binPt,binEta,xData,eData,xMC))
            inputFile.Close()

    effData_2D = ROOT.TH2D('effData','',nbinsPt,array('d',list(ptBins)),nbinsEta,array('d',list(etaBins)))
    effMC_2D = ROOT.TH2D('effMC','',nbinsPt,array('d',list(ptBins)),nbinsEta,array('d',list(etaBins)))
    for iEta in range(1,nbinsEta+1):
        binEta = '%1i'%(iEta)
        effData_1D = ROOT.TH1D('effData_eta'+binEta,'',nbinsPt,array('d',list(ptBins)))
        effMC_1D = ROOT.TH1D('effMC_eta'+binEta,'',nbinsPt,array('d',list(ptBins)))
        for iPt in range(1,nbinsPt+1):

            effData_1D.SetBinContent(iPt,effData[binEta][iPt-1])
            effData_1D.SetBinError(iPt,errData[binEta][iPt-1])

            effData_2D.SetBinContent(iPt,iEta,effData[binEta][iPt-1])
            effData_2D.SetBinError(iPt,iEta,errData[binEta][iPt-1])

            effMC_1D.SetBinContent(iPt,effMC[binEta][iPt-1])
            effMC_1D.SetBinError(iPt,errMC[binEta][iPt-1])

            effMC_2D.SetBinContent(iPt,iEta,effMC[binEta][iPt-1])
            effMC_2D.SetBinError(iPt,iEta,errMC[binEta][iPt-1])
            
        PlotTagProbe(effData_1D,effMC_1D,era=era,channel=chan,binEta=binEta,ymin=ymin,ymax=ymax)


    lepType = 'TauE'
    if chan=='mt':
        lepType = 'TauMu'
        
    outputFileName = '%s/ScaleFactors/SF_%s_%s.root'%(basedir,lepType,era)
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    effData_2D.Write('effData')
    effMC_2D.Write('effMC')
    outputFile.Close()
