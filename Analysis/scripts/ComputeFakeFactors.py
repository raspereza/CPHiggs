#! /usr/bin/env python3
# Author: Alexei Raspereza (October 2025)
# Measurement of fake factors
###########################################
import ROOT
import math
from array import array
import os
import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

stop_point = 70.

def extractHistos(f,regions):
    hists = {}
    for sample in utils.samples:
        for region in regions:
            for dm in utils.dm_labels:
                for njets in utils.njets_labels:
                    for eta in utils.eta_labels:
                        for tauid in utils.tauid_labels:
                            for typ in utils.type_labels:
                                name = '%s_%s_%s_%s_%s_%s_%s'%(sample,region,dm,njets,eta,tauid,typ)
                                nameInput = '%s_pt_2_%s_%s_%s_%s_%s_%s'%(sample,region,dm,njets,eta,tauid,typ)
                                hists[name] = f.Get(nameInput)
    return hists
    
# Fitting function: pol2
def FitFuncPol2(x,par):
    b = par[0]+par[1]*x[0]+par[2]*x[0]*x[0]
    return b

# Fitting function: (pol3)
def FitFunc(x,par):
    arg = x[0]
    if arg>stop_point:
        arg=stop_point
    b = par[0]+par[1]*arg+par[2]*arg*arg+par[3]*arg*arg*arg
    return b

# uncertainty band
def Band(x,par):
    W = 0
    xx = 1
    arg = x[0]
    if arg>stop_point:
        arg=stop_point
    for i in range(0,7):
        W += par[i]*xx
        xx *= arg
    return ROOT.TMath.Sqrt(W)

# uncertainty band (pol2)
def BandPol2(x,par):
    W = 0
    xx = 1
    for i in range(0,5):
        W += par[i]*xx
        xx *= x[0]
    return ROOT.TMath.Sqrt(W)

# upper uncertainty 
def FitFuncPlus(x,par):
    parBand = []
    # shift index parameter by 4
    for i in range(0,7):
        parBand.append(par[i+4])
    return FitFunc(x,par)+Band(x,parBand)

# lower uncertainty 
def FitFuncMinus(x,par):
    parBand = []
    # shift index parameter by 4
    for i in range(0,7):
        parBand.append(par[i+4])
    ff = FitFunc(x,par)-Band(x,parBand)
    if ff<0: ff=0
    return ff

# upper uncertainty 
def FitFuncPlusPol2(x,par):
    parBand = []
    # shift index parameter by 3
    for i in range(0,5):
        parBand.append(par[i+3])
    return FitFunc(x,par)+BandPol2(x,parBand)

# lower uncertainty 
def FitFuncMinusPol2(x,par):
    parBand = []
    # shift index parameter by 3
    for i in range(0,5):
        parBand.append(par[i+3])
    return FitFunc(x,par)-BandPol3(x,parBand)

def PlotFF(hists,bins,**kwargs):

    era = kwargs.get('era','Run3')
    chan = kwargs.get('channel','mt')
    region = kwargs.get('region','qcd')
    dm = kwargs.get('dm','pi')
    njets = kwargs.get('njets','njets0')
    eta = kwargs.get('eta','barrel')
    sampleToProcess = kwargs.get('sample','data')
    ipcut = kwargs.get('ipcut',True)
    
    suffix = '%s_%s_%s_%s'%(region,dm,njets,eta)

    data_num = None
    data_den = None
    n_data_num = 0.
    n_data_den = 0.
    if sampleToProcess=='data':
    
        data_num = hists['data_'+suffix+'_nominal_all'].Clone('data_num'+suffix)
        data_den = hists['data_'+suffix+'_inverted_all'].Clone('data_den'+suffix)

        # subtracting leptons
        for sample in ['zll','ztt','top','vv','wjets']:
            for typ in ['lep','tau']:
                data_num.Add(data_num,hists[sample+'_'+suffix+'_nominal_'+typ],1.,-1.)
                data_den.Add(data_den,hists[sample+'_'+suffix+'_inverted_'+typ],1.,-1.)

        mc_had_num = hists['zll_'+suffix+'_nominal_had'].Clone('mc_had_num')
        mc_had_den = hists['zll_'+suffix+'_inverted_had'].Clone('mc_had_den')
        for sample in ['ztt','top','vv','wjets']:
            name = '%s_%s_nominal_had'%(sample,suffix)
            mc_had_num.Add(mc_had_num,hists[name],1.,1.)
            name = '%s_%s_inverted_had'%(sample,suffix)
            mc_had_den.Add(mc_had_num,hists[name],1.,1.)

        n_mc_num = mc_had_num.GetSumOfWeights()
        n_mc_den = mc_had_den.GetSumOfWeights()
        n_data_num = data_num.GetSumOfWeights()
        n_data_den = data_den.GetSumOfWeights()
        f_mc_num = n_mc_num/max(n_data_num,1.0)
        f_mc_den = n_mc_den/max(n_data_den,1.0)
    else:
        data_num = hists[sampleToProcess+'_'+suffix+'_nominal_had'].Clone(sampleToProcess+'_num'+suffix)
        data_den = hists[sampleToProcess+'_'+suffix+'_inverted_had'].Clone(sampleToProcess+'_den'+suffix)
        n_data_num = data_num.GetSumOfWeights()
        n_data_den = data_den.GetSumOfWeights()
        
    
    num_x = utils.rebinHisto(data_num,bins,'rebinned')
    den_x = utils.rebinHisto(data_den,bins,'rebinned')
    newbins = []
    nb = len(bins)-1
    for ib in range(0,len(bins)-1):
        newbins.append(bins[ib])
    newbins.append(130.)

    num = ROOT.TH1D('num'+suffix,'',nb,array('d',list(newbins)))
    den = ROOT.TH1D('den'+suffix,'',nb,array('d',list(newbins)))
    for ib in range(1,nb+1):
        num.SetBinContent(ib,num_x.GetBinContent(ib))
        num.SetBinError(ib,num_x.GetBinError(ib))
        den.SetBinContent(ib,den_x.GetBinContent(ib))
        den.SetBinError(ib,den_x.GetBinError(ib))
    
    ff_hist = utils.divideHistos(num,den,'ff_'+suffix+'_'+sampleToProcess)
    
    styles.InitData(ff_hist)

    nbins = ff_hist.GetNbinsX()
    fmin = ff_hist.GetBinLowEdge(1)
    fmax = ff_hist.GetBinLowEdge(nbins+1)

    average = 0.
    maximum = 0.
    for ib in range(1,nbins+1):
        x = ff_hist.GetBinContent(ib)
        e = ff_hist.GetBinError(ib)
        s = x+e
        if s>maximum: maximum=s
        average += x
        
    average /= float(nbins)
    ff_hist.GetYaxis().SetRangeUser(0.,2*maximum)
    
    fitSF = ROOT.TF1('fitFunc_'+suffix+'_'+sampleToProcess,FitFunc,fmin,fmax,4)
    fitSF.SetLineColor(ROOT.kBlue)
    fitSF.SetParameter(0,average)
    fitSF.SetParameter(1,0.0)
    fitSF.SetParameter(2,0.0)
    fitSF.SetParameter(3,0.0)
    fitSF.SetParName(0,'a0')
    fitSF.SetParName(1,'a1')
    fitSF.SetParName(2,'a2')
    fitSF.SetParName(3,'a3')

    print('')
    print('+++++++++++++++++++++++++++++++++++++++')
    if sampleToProcess=='data':
        print('%s %s %s %s'%(region,dm,njets,eta))
        print('direct  : data = %6.0f  mc = %6.0f  f = %5.3f'%(n_data_num,n_mc_num,f_mc_num))
        print('iverted : data = %6.0f  mc = %6.0f  f = %5.3f'%(n_data_den,n_mc_den,f_mc_den))
    else:
        print('%s %s %s %s'%(region,dm,njets,eta))
        print('direct  : %s = %6.0f'%(sampleToProcess,n_data_num))
        print('iverted : %s = %6.0f'%(sampleToProcess,n_data_den))
    
    dummy = styles.MakeCanvas('dummy','',400,400)
    fitRes = ff_hist.Fit(fitSF,"S")

    cov = fitRes.GetCovarianceMatrix()
    err_par0 = fitSF.GetParError(0)
    err_par1 = fitSF.GetParError(1)
    err_par2 = fitSF.GetParError(2)
    err_par3 = fitSF.GetParError(3)
    corr00 = cov(0,0)/(err_par0*err_par0)
    corr01 = cov(0,1)/(err_par0*err_par1)
    corr02 = cov(0,2)/(err_par0*err_par2)
    corr03 = cov(0,3)/(err_par0*err_par3)
    corr11 = cov(1,1)/(err_par1*err_par1)
    corr12 = cov(1,2)/(err_par1*err_par2)
    corr13 = cov(1,3)/(err_par1*err_par3)
    corr22 = cov(2,2)/(err_par2*err_par2)
    corr23 = cov(2,3)/(err_par2*err_par3)
    corr33 = cov(3,3)/(err_par3*err_par3)

    # upper uncertainty
    # 4 + 7 = 11 parameters
    fitSF_up = ROOT.TF1('FitFuncPlus_'+suffix+'_'+sampleToProcess,FitFuncPlus,fmin,fmax,11)
    fitSF_up.SetLineColor(ROOT.kRed)
    fitSF_up.SetParameter(0,fitSF.GetParameter(0))
    fitSF_up.SetParameter(1,fitSF.GetParameter(1))
    fitSF_up.SetParameter(2,fitSF.GetParameter(2))
    fitSF_up.SetParameter(3,fitSF.GetParameter(3))
    fitSF_up.SetParameter(4,cov(0,0))
    fitSF_up.SetParameter(5,2*cov(0,1))
    fitSF_up.SetParameter(6,2*cov(0,2)+cov(1,1))
    fitSF_up.SetParameter(7,2*(cov(0,3)+cov(1,2)))
    fitSF_up.SetParameter(8,2*cov(1,3)+cov(2,2))
    fitSF_up.SetParameter(9,2*cov(2,3))
    fitSF_up.SetParameter(10,cov(3,3))
    
    # lower uncertainty
    fitSF_down = ROOT.TF1('FitFuncMinus_'+suffix+'_'+sampleToProcess,FitFuncMinus,fmin,fmax,11)
    fitSF_down.SetLineColor(ROOT.kRed)
    for i in range(0,11):
        fitSF_down.SetParameter(i,fitSF_up.GetParameter(i))
    
    #    print('')    
    #    print('Correlation matrix ->')
    #    print('%8.5f  %8.5f  %8.5f  %8.5f'%(corr00,corr01,corr02,corr03))
    #    print('%8.5f  %8.5f  %8.5f  %8.5f'%(corr01,corr11,corr12,corr13))
    #    print('%8.5f  %8.5f  %8.5f  %8.5f'%(corr02,corr12,corr22,corr23))
    #    print('%8.5f  %8.5f  %8.5f  %8.5f'%(corr03,corr13,corr23,corr33))
    
    hfit = ROOT.TH1D("hfit_"+suffix+"_"+sampleToProcess,"",200,fmin,fmax)
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hfit,0.68)
    styles.InitModel(hfit,"p_{T} (GeV)","Fake Factors",4)
    hfit.SetFillColor(ROOT.kCyan)
    hfit.SetFillStyle(1001)
    hfit.SetLineWidth(0)
    hfit.SetLineColor(4)
    hfit.SetMarkerSize(0)
    hfit.SetMarkerStyle(0)
    hfit.GetYaxis().SetTitle("fake factor")
    ff_hist.GetYaxis().SetTitle('fake factor')
    ff_hist.GetXaxis().SetTitle('p_{T} (GeV)')
    ff_hist.GetXaxis().SetMoreLogLabels(True)
    ff_hist.GetXaxis().SetNoExponent()
    ff_hist.GetXaxis().SetTitleOffset(1.2)

    # canvas and pads
    canvas = styles.MakeCanvas("canv_"+suffix,"",600,700)

    ff_hist.Draw('e1')
    hfit.Draw('e2same')
    fitSF.Draw('lsame')
    fitSF_up.Draw('lsame')
    fitSF_down.Draw('lsame')
    ff_hist.Draw('e1same')
    
    leg = ROOT.TLegend(0.2,0.8,0.5,0.9)
    styles.SetLegendStyle(leg)
    leg.SetTextSize(0.04)
    leg.SetHeader(chan+' '+suffix)

    leg.Draw()
    styles.CMS_label(canvas,era=era)
    canvas.SetLogx(True)
    canvas.Update()

    subfolder = 'FakeFactors'
    if ipcut:
        subfolder = 'FakeFactors_ipcut'
    if sampleToProcess=='data':
        figure_folder = '/eos/home-r/rasp/php-plots/plots/%s/%s/data/%s'%(subfolder,chan,region)
    else:
        figure_folder = '/eos/home-r/rasp/php-plots/plots/%s/%s/mc/%s'%(subfolder,chan,sampleToProcess)
    outputGraphics = '%s/FF_%s_%s_%s.png'%(figure_folder,dm,njets,eta)    
    canvas.Print(outputGraphics)
    return hfit,fitSF,fitSF_up,fitSF_down
    
if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3', choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et'])
    parser.add_argument('-ipcut','--ipcut', dest='ipcut', action='store_true')
    parser.add_argument('-mc','--mc',dest='mc',action='store_true')
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    
    bins_pt = [20, 25, 30, 35, 40, 50, 70, 200.]
    
    basedir = '%s'%(utils.outputFolder)
    
    suffixFile = 'x'
    ipcut = False
    if args.ipcut:
        suffixFile = 'x_ipcut1'
        ipcut = True
    suffixOutFile = suffixFile
    if args.mc:
        suffixOutFile += '_mc'
        
    inputFileName = '%s/selection/jetFakes/%s_%s_%s.root'%(basedir,chan,era,suffixFile)
    inputFile = ROOT.TFile(inputFileName,'read')
    print('')
    print(inputFile)
    print('')
    regions = ['qcd','wj','ss_antiiso','top','os_antiiso']
    hists = extractHistos(inputFile,regions)

    hfit = {}
    fitSF = {}
    fitSF_up = {}
    fitSF_down = {}

    if args.mc:
        sample='wjets'
        region='wj'
        for dm in utils.dm_labels:
            for njets in utils.njets_labels:
                for eta in utils.eta_labels:
                    suffix = 'mc_%s_%s_%s_%s'%(region,dm,njets,eta)
                    hfit[suffix],fitSF[suffix],fitSF_up[suffix],fitSF_down[suffix] = PlotFF(hists,bins_pt,
                                                                                            era=era,
                                                                                            channel=chan,
                                                                                            region=region,
                                                                                            dm=dm,
                                                                                            njets=njets,
                                                                                            sample=sample,
                                                                                            eta=eta,
                                                                                            ipcut=ipcut)
                
        sample = 'top'
        region = 'top'
        for dm in utils.dm_labels:
            for njets in utils.njets_labels:
                for eta in utils.eta_labels:
                    suffix = 'mc_%s_%s_%s_%s'%(region,dm,njets,eta)
                    hfit[suffix],fitSF[suffix],fitSF_up[suffix],fitSF_down[suffix] = PlotFF(hists,bins_pt,
                                                                                            era=era,
                                                                                            channel=chan,
                                                                                            region=region,
                                                                                            dm=dm,
                                                                                            njets=njets,
                                                                                            sample=sample,
                                                                                            eta=eta,
                                                                                            ipcut=ipcut)
    else:
        sample = 'data'
        for region in ['os_antiiso']:
#        for region in ['qcd','wj','ss_antiiso','os_antiiso']:

            for dm in utils.dm_labels:
                for njets in utils.njets_labels:
                    for eta in utils.eta_labels:
                        suffix = '%s_%s_%s_%s'%(region,dm,njets,eta)
                        hfit[suffix],fitSF[suffix],fitSF_up[suffix],fitSF_down[suffix] = PlotFF(hists,bins_pt,
                                                                                                era=era,
                                                                                                channel=chan,
                                                                                                region=region,
                                                                                                dm=dm,
                                                                                                njets=njets,
                                                                                                sample=sample,
                                                                                                eta=eta,
                                                                                                ipcut=ipcut)

                        
    outputFolder='%s/src/IPcorrectionsRun3/FakeFactors/data'%(os.getenv('CMSSW_BASE'))
    outputFileName = '%s/FF_%s_%s_%s.root'%(outputFolder,era,chan,suffixOutFile)
    print('opening file %s'%(outputFileName))
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    for h in hfit:
        hfit[h].Write('hfit_'+h)
        fitSF[h].Write('fitFunc_'+h+'_nom')
        fitSF_up[h].Write('fitFunc_'+h+'_up')
        fitSF_down[h].Write('fitFunc_'+h+'_down')

    outputFile.Close()
