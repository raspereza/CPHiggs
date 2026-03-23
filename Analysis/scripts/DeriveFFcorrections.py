#! /usr/bin/env python3
# Author: Alexei Raspereza (November 2025)
# Determination of systematic uncertainties
# in jet->tau fake background
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def getBinning(chan,closure):
    xmin = 0.3
    xmax = 1.0
    nbins = 7
    if chan=='et' and closure=='QCD':
        nbins = 3
    return nbins,xmin,xmax


fit_dict = {
    'mt' : {'QCD': {'bdt_ditau':0,  'bdt_fakes':1, 'bdt_signal':2},
            'W'  : {'bdt_ditau':1,  'bdt_fakes':3, 'bdt_signal':1},
            'MC' : {'bdt_ditau':-1, 'bdt_fakes':3, 'bdt_signal':2},
            },
    'et' : {'QCD': {'bdt_ditau':1,  'bdt_fakes':1, 'bdt_signal':1},
            'W'  : {'bdt_ditau':1,  'bdt_fakes':1, 'bdt_signal':0},
            'MC' : {'bdt_ditau':-1, 'bdt_fakes':3, 'bdt_signal':-1},
            },
}

def FitFuncPol0(x,par):
    return par[0]

def FitFuncPol1(x,par):
    b = par[0]+par[1]*x[0]
    return b

def FitFuncPol2(x,par):
    b = par[0]+par[1]*x[0]+par[2]*x[0]*x[0]
    return b

def FitFuncPol3(x,par):
    b = par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]
    return b
    
def ClosureMC(hists):
    
    hist_out = {}

    hist_out['Num'] = hists['data_wj_ff_wj_all'].Clone('Num_MC')
    hist_out['Den'] = hists['wjets_wj_ff_mc_wj_had'].Clone('Den_MC')
    totalData_AR = hists['data_wj_ff_ar_all'].GetSumOfWeights()
    totalW_AR = hists['wjets_wj_ff_ar_had'].GetSumOfWeights()
    for s in utils.bkg_samples:
        for t in ['lep','tau']:
            name = f'{s}_wj_ff_wj_{t}'
            hist_out['Num'].Add(hist_out['Num'],hists[name],1.,-1.)
            name = f'{s}_wj_ff_ar_{t}'
            totalData_AR -= hists[name].GetSumOfWeights()
    fraction = totalW_AR/totalData_AR
    hist_out['Num'].Scale(fraction)

    return hist_out
    
def ClosureQCD(hists):

    hist_out = {}
    hist_out['Num'] = hists['data_qcd_closure_ff_qcd_all'].Clone('Num_QCD')
    hist_out['Den'] = hists['data_qcd_closure_ff_all'].Clone('Den_QCD')
    for s in utils.bkg_samples:
        for t in ['lep','tau']:
            name = f'{s}_qcd_closure_ff_qcd_{t}'
            hist_out['Num'].Add(hist_out['Num'],hists[name],1.,-1.)
            name =  f'{s}_qcd_closure_ff_{t}'
            hist_out['Den'].Add(hist_out['Den'],hists[name],1.,-1.)
    return hist_out

def ClosureW(hists):
    hist_out = {}
    hist_out['Den'] = hists['wjets_lowmt_os_iso_mc_wj_had'].Clone('Num_W')
    hist_out['Num'] = hists['wjets_lowmt_os_iso_had'].Clone('Den_W')

    return hist_out

def ExtractHistosFF(f,var):
    hists = {}

    for sample in utils.samples:
        for region in utils.region_labels:
            for typ in utils.type_labels:
                nameInput='%s_%s_%s_%s'%(sample,var,region,typ)
                name='%s_%s_%s'%(sample,region,typ)
                hists[name]=f.Get(nameInput)
                for ff in utils.ff_labels:
                    nameInput='%s_%s_%s_%s_%s'%(sample,var,region,ff,typ)
                    name='%s_%s_%s_%s'%(sample,region,ff,typ)
                    hists[name] = f.Get(nameInput)

    return hists
    
def Plot(hists,bins,**kwargs):

    era = kwargs.get('era','Run3')
    var = kwargs.get('var','m_vis')
    chan = kwargs.get('channel','mt')
    ymin = kwargs.get('ymin',0.501)
    ymax = kwargs.get('ymax',1.499)
    closure = kwargs.get('closure','QCD')
    
    h_data_x = hists['Num'].Clone('h_data')
    h_model_x = hists['Den'].Clone('h_model')
    h_data = utils.rebinHisto(h_data_x,bins,'rebinned')
    h_model = utils.rebinHisto(h_model_x,bins,'rebinned')
    
    nbins = h_model.GetNbinsX()
    err2_data = 0.
    err2_model = 0.
    x_data = 0.
    x_model = 0.    
    fmin = h_model.GetBinLowEdge(1)
    fmax = h_model.GetBinLowEdge(nbins+1)
    for ib in range(nbins):
        e_data = h_data.GetBinError(ib+1)
        e_model = h_model.GetBinError(ib+1)
        err2_data += e_data*e_data
        err2_model += e_model*e_model
        x_data += h_data.GetBinContent(ib+1)
        x_model += h_model.GetBinContent(ib+1)
    e_data = ROOT.TMath.Sqrt(err2_data)
    e_model = ROOT.TMath.Sqrt(err2_model)
    
    ratio = x_data/x_model
    r_data = e_data/x_data
    r_model = e_model/x_model
    r_ratio = ROOT.TMath.Sqrt(r_model*r_model+r_data*r_data)
    eratio = ratio*r_ratio

    print('')
    print('Yields %s ->'%(var))
    print('Num   : %7.0f'%(x_data))
    print('Den   : %7.0f'%(x_model))
    print('')

    
    styles.InitData(h_data)
    xtitle = utils.XTitle[chan][var]
    styles.InitData(h_model)
    h_model.SetMarkerColor(2)
    h_model.SetLineColor(2)
    h_model.SetMarkerStyle(21)
    
    h_tot = h_model.Clone("total")
    styles.InitTotalHist(h_tot)
    
    h_ratio = utils.histoErrRatio(h_data,h_tot,'ratio')
    h_tot_ratio = utils.createUnitHisto(h_tot,'tot_ratio')

    styles.InitRatioHist(h_ratio)
    h_ratio.GetYaxis().SetRangeUser(ymin,ymax)
    
    YMax = h_data.GetMaximum()
    if h_tot.GetMaximum()>YMax: YMax = h_tot.GetMaximum()
    h_data.GetYaxis().SetRangeUser(0.,1.5*YMax)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetYaxis().SetTitle("events")
    
    h_ratio.GetYaxis().SetTitle("ratio")
    h_ratio.GetXaxis().SetTitle(xtitle)
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetYaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetLabelSize(0.07)
    h_ratio.GetYaxis().SetLabelSize(0.07)
#    h_ratio.GetYaxis().SetTitleOffset(1.)
    
    # canvas and pads
    name_canv = 'canv_'+closure+'_'+var
    canvas = styles.MakeCanvas(name_canv,"",600,700)
    # upper pad
    upper = ROOT.TPad("upper", "pad",0,0.41,1,1)
    upper.Draw()
    upper.cd()
    styles.InitUpperPad(upper)    
    
    h_data.Draw('e1')
    h_model.Draw('e1same')
    h_data.Draw('e1same')
    
    leg = ROOT.TLegend(0.25,0.65,0.45,0.85)
    styles.SetLegendStyle(leg)
    leg.SetHeader(closure+' closure')
    leg.SetTextSize(0.045)
    leg.AddEntry(h_data,'num','lp')
    leg.AddEntry(h_model,'den','lp')
    leg.Draw()

    styles.CMS_label(upper,era=era)

    upper.Draw("SAME")
    upper.RedrawAxis()
    upper.Modified()
    upper.Update()
    canvas.cd()

    # lower pad
    lower = ROOT.TPad("lower", "pad",0,0,1,0.40)
    lower.Draw()
    lower.cd()
    styles.InitLowerPad(lower)
    
    fitFunc = None
    opt = fit_dict[chan][closure][var]
    funcName = f'fitFunc_{var}_{closure}'
    print(f'Fitting with option {opt}')
    if opt<=0:
        fitFunc = ROOT.TF1(funcName,FitFuncPol0,fmin,fmax,1)
    elif opt==1:
        fitFunc = ROOT.TF1(funcName,FitFuncPol1,fmin,fmax,2)
    elif opt==2:
        fitFunc = ROOT.TF1(funcName,FitFuncPol2,fmin,fmax,3)
    elif opt==3:
        fitFunc = ROOT.TF1(funcName,FitFuncPol3,fmin,fmax,4)
        
    fitFunc.SetLineColor(ROOT.kBlue)
    fitFunc.SetParameter(0,ratio)
    fitFunc.SetParameter(1,0.0)
    fitFunc.SetParameter(2,0.0)
    fitFunc.SetParameter(3,0.0)
    fitFunc.SetParameter(4,0.0)
    
    hfit = ROOT.TH1D("hfit_"+var+"_"+closure,"",200,fmin,fmax)
    if opt>=0:
        h_ratio.Fit(fitFunc)
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hfit,0.68)
    else:
        for ib in range(200):
            hfit.SetBinContent(ib+1,ratio)
            hfit.SetBinError(ib+1,eratio)

    styles.InitModel(hfit,xtitle,"ratio",4)
    hfit.SetFillColor(ROOT.kCyan)
    hfit.SetFillStyle(1001)
    hfit.SetLineWidth(0)
    hfit.SetLineColor(4)
    hfit.SetMarkerSize(0)
    hfit.SetMarkerStyle(0)

    h_ratio.Draw('e1')
    hfit.Draw('e2same') 
    fitFunc.Draw('lsame')
    h_ratio.Draw('e1same')

    lower.Modified()
    lower.RedrawAxis()

    canvas.cd()
    canvas.Modified()
    canvas.cd()
    canvas.SetSelected(canvas)
    canvas.Update()
    print('')
    outputFolder = '/eos/home-r/rasp/php-plots/plots/FFcorrections_%s'%(chan)
    outputGraphics = '%s/%s_%s.png'%(outputFolder,var,closure)
    canvas.Print(outputGraphics)
    return fitFunc

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era','--era',dest='era',default='Run3',choices=['Run3_2022','Run3_2023','Run3'])
    parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et'])
    parser.add_argument('-ymin','--ymin',dest='ymin',type=float,default=0.201)
    parser.add_argument('-ymax','--ymax',dest='ymax',type=float,default=1.999)
    parser.add_argument('-suffix','--suffix',dest='suffix',default='x_ipcut1_ff_ipcut')
    
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    ymin = args.ymin
    ymax = args.ymax
    suffix = args.suffix
    plotLegend = True

    cmssw_base = os.getenv('CMSSW_BASE')
    basedir = utils.outputFolder+'/selection'


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

    variables = {
        'bdt_fakes': 'fakes',
        'bdt_ditau': 'ditau',
        'bdt_signal': 'signal'
    }
    closures = {
        'QCD': 'qcd',
        'MC': 'mc_top',
        'W' : 'wj'
    }

    funcs = {}
    for var in variables:
        hists = ExtractHistosFF(inputFile,var)
        for closure in closures:
            hist_out = {}
            if closure=='QCD':
                hist_out = ClosureQCD(hists)
            if closure=='W':
                hist_out = ClosureW(hists)
            if closure=='MC':
                hist_out = ClosureMC(hists)
            nbins,xmin,xmax = getBinning(chan,closure)
            bins = utils.createBins(nbins,xmin,xmax)
            ff_name = closures[closure]
            cat_name = variables[var]
            name = f'{ff_name}_{cat_name}'
            funcs[name] = Plot(hist_out,bins,era=era,var=var,channel=chan,
                               ymin=ymin,ymax=ymax,closure=closure)

    outputFileName = f'{cmssw_base}/src/IPcorrectionsRun3/FakeFactors/data/FF_closure_{chan}.root'
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputFile.cd('')
    for func in funcs:
        funcs[func].Write(func)
    outputFile.Close()
