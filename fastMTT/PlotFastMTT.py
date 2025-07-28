#! /usr/bin/env python3
# Author: Alexei Raspereza (December 2024)
# Plotting macro to test FastMTT
import ROOT
import math
from array import array
import os

import CPHiggs.Analysis.styles as styles
import CPHiggs.Analysis.utils as utils

def Plot(hist1,hist2,**kwargs):

    sample = kwargs.get('sample','higgs')
    era = kwargs.get('era','Run3_2022')
    chan = kwargs.get('chan','mt')
    isDY = kwargs.get('isDY',True) 
    plot = kwargs.get('plot',0) # 0 - mass, 1 = dpt1, 2 = dpt2
    isIC = kwargs.get('isIC',False) # is IC tuple ?
    isMass = plot==0
    
    # histograms
    h_even = hist1
    h_odd = hist2

    xtitle = 'p_{T,1}^{rec}/p_{T,1}^{gen}'
    leg_even = 'w/o m_{H}'
    leg_odd = 'with m_{H}'
    if isMass:
        xtitle = 'mass (GeV)'
        leg_even = 'm_{vis}'
        leg_odd  = 'm_{#tau#tau}'
    else:
        if plot==2:
            xtitle = 'p_{T,2}^{rec}/p_{T,2}^{gen}'

    header = 'H#rightarrow #tau_{#mu}#tau_{h}'
    if chan=='tt':
        header = 'H#rightarrow#tau_{h}#tau_{h}'
    if isDY:
        header = 'Z#rightarrow#tau_{#mu}#tau_{h}'
        if chan=='tt':
            header = 'Z#rightarrow#tau_{h}#tau_{h}'
        
    ytitle = 'normalised'

    styles.InitModel(h_even,xtitle,ytitle,2)
    styles.InitModel(h_odd,xtitle,ytitle,4)

    utils.zeroBinErrors(h_even)
    utils.zeroBinErrors(h_odd)
    
    YMax = h_even.GetMaximum()
    if h_odd.GetMaximum()>YMax: YMax = h_odd.GetMaximum()
    
    h_even.GetYaxis().SetRangeUser(0.,1.1*YMax)

    # canvas and pads
    canv_name = 'canv_mass'
    if plot==1:
        canv_name = 'canv_dpt1'
    if plot==2:
        canv_name = 'canv_dpt2'
    canvas = styles.MakeCanvas(canv_name,"",800,700)
    
    h_even.Draw('h')
    h_odd.Draw('hsame')

    leg = ROOT.TLegend(0.7,0.5,0.9,0.7)
    styles.SetLegendStyle(leg)
    leg.SetHeader(header)
    leg.SetTextSize(0.05)
    leg.AddEntry(h_even,leg_even,'l')
    leg.AddEntry(h_odd,leg_odd,'l')
    leg.Draw()

    styles.CMS_label(canvas,era=era,extraText='Simulation')

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    suffix = 'mass'
    if plot==1:
        suffix = 'dpt1'
    if plot==2:
        suffix = 'dpt2'
    outputGraphics = '%s_%s_%s_%s'%(sample,era,chan,suffix)
    if isIC:
        outputGraphics += '_tuple'
    outputGraphics += '.png'
    canvas.Print(outputGraphics)

if __name__ == "__main__":

    styles.InitROOT()
    styles.SetStyle()
    
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-era' ,'--era', dest='era', default='Run3_2022', choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
    parser.add_argument('-channel','--channel', dest='channel', default='mt',choices=['mt','et','tt'])
    parser.add_argument('-sample','--sample',dest='sample',default='higgs',choices=['higgs','dy'])
    parser.add_argument('-ICTuple','--ICTuple',dest='icTuple',action='store_true')
    args = parser.parse_args()

    era = args.era
    chan = args.channel
    sample = args.sample
    isDY = sample=='dy'
    isIC = args.icTuple
    
    name_mvis = 'mvis'
    name_mtt = 'mtt'
    name_dpt1 = 'dpt1'
    name_dpt1_BW = 'dpt1_BW'
    name_dpt2 = 'dpt2'
    name_dpt2_BW = 'dpt2_BW'
    if args.icTuple:
        name_mtt += '_nom'
        name_dpt1 += '_nom'
        name_dpt1_BW += '_nom'
        name_dpt2 += '_nom'
        name_dpt2_BW += '_nom'
        
    filename = '%s_%s_%s.root'%(sample,era,chan)
    inputfile = ROOT.TFile(filename,'read')
    h_mvis = inputfile.Get('mvis') 
    h_mtt = inputfile.Get(name_mtt)
    h_dpt1 = inputfile.Get(name_dpt1)
    h_dpt1_BW = inputfile.Get(name_dpt1_BW)
    h_dpt2 = inputfile.Get(name_dpt2)
    h_dpt2_BW = inputfile.Get(name_dpt2_BW)
    
    h_mvis.Scale(1./h_mvis.GetSumOfWeights())
    h_mtt.Scale(1./h_mtt.GetSumOfWeights())
    h_dpt1.Scale(1./h_dpt1.GetSumOfWeights())
    h_dpt1_BW.Scale(1./h_dpt1_BW.GetSumOfWeights())
    h_dpt2.Scale(1./h_dpt2.GetSumOfWeights())
    h_dpt2_BW.Scale(1./h_dpt2_BW.GetSumOfWeights())

    mean_dpt1 = h_dpt1.GetMean()
    rms_dpt1 = h_dpt1.GetRMS()
    mean_dpt1_BW = h_dpt1_BW.GetMean()
    rms_dpt1_BW = h_dpt1_BW.GetRMS()

    mean_dpt2 = h_dpt2.GetMean()
    rms_dpt2 = h_dpt2.GetRMS()
    mean_dpt2_BW = h_dpt2_BW.GetMean()
    rms_dpt2_BW = h_dpt2_BW.GetRMS()

    print('')
    print('pt1  :   w/o mH : with mH  ')
    print('Mean :   %5.3f  :  %5.3f'%(mean_dpt1,mean_dpt1_BW))
    print('RMS  :   %5.3f  :  %5.3f'%(rms_dpt1,rms_dpt1_BW))
    print('')
    print('pt2  :   w/o mH : with mH  ')
    print('Mean :   %5.3f  :  %5.3f'%(mean_dpt2,mean_dpt2_BW))
    print('RMS  :   %5.3f  :  %5.3f'%(rms_dpt2,rms_dpt2_BW))
    print('')

    Plot(h_mvis,h_mtt,sample=sample,era=era,isDY=isDY,chan=chan,plot=0,isIC=isIC)
    Plot(h_dpt1,h_dpt1_BW,sample=sample,era=era,isDY=isDY,chan=chan,plot=1,isIC=isIC)
    Plot(h_dpt2,h_dpt2_BW,sample=sample,era=era,isDY=isDY,chan=chan,plot=2,isIC=isIC)
    
    
