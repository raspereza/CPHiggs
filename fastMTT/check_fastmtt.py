###########################################
#
# Author: Alexei Raspereza (June 2025)
# Code to validate FastMTT programme
###########################################

import numpy as np
#from fastmtt_cpp import fastmtt_cpp
import fastmtt_cpp
import ROOT 
import CPHiggs.Analysis.utils as utils
from argparse import ArgumentParser

def fill_hist(hist, array):
    [hist.Fill(x) for x in array]


################
# Main routine #
################
parser = ArgumentParser()
parser.add_argument('-era','--era',dest='era',default='Run3_2022',choices=['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'])
parser.add_argument('-channel','--channel',dest='channel',default='mt',choices=['mt','et','tt'])
parser.add_argument('-sample','--sample',dest='sample',default='GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay',choices=["GluGluHTo2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay","DYto2L_M_50_madgraphMLM","GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay"])
parser.add_argument('-nevts','--nevts',dest='nevts',type=int,default=10000)

print('Beginning')
    
args = parser.parse_args()
dirname='/eos/cms/store/group/phys_tau/ksavva/For_Aliaksei/files/testingzpt'
filename=dirname+'/'+args.era+'/'+args.channel+"/"+args.sample+"/nominal/merged.root"

print('opening %s'%(filename))    
df = ROOT.RDataFrame("ntuple",filename)

print('defining cuts')
cuts = 'os>0.5&&idDeepTau2018v2p5VSe_2>=6&&idDeepTau2018v2p5VSmu_2>=4&&idDeepTau2018v2p5VSjet_2>=7&&pt_2>20.&&fabs(eta_2)<2.3'

if args.channel=='mt':
    cuts += '&&iso_1<0.15&&pt_1>25&&fabs(eta_1)<2.4'
if args.channel=='et':
    cuts += '&&iso_1<0.10&&pt_1>32&&fabs(eta_1)<2.1'
if args.channel=='tt':
    cuts += '&&idDeepTau2018v2p5VSe_1>=6&&idDeepTau2018v2p5VSmu_1>=4&&idDeepTau2018v2p5VSjet_1>=7&&pt_1>40.&&pt_2>40.&&fabs(eta_1)<2.3'
    
print('reading tuple as numpy columns')
cols = df.Filter(cuts).AsNumpy(["pt_1","pt_2","eta_1","eta_2","phi_1","phi_2","mass_1","mass_2","puppimet","puppimetphi","puppimetcov00","puppimetcov01","puppimetcov11","m_vis","FastMTT_mass","FastMTT_pt_1","FastMTT_pt_2","FastMTT_pt_1_constraint","FastMTT_pt_2_constraint","genPart_pt_1","genPart_pt_2"])

print('Length of %: \n',len(cols["pt_1"]))

# decay_type = 0 : tau -> electron
#            = 1 : tau -> muon
#            = 2 : tau -> hadrons

decay_type_1 = 1*np.ones(len(cols["pt_1"]),dtype=np.int32)
decay_type_2 = 2*np.ones(len(cols["pt_1"]),dtype=np.int32)
if args.channel=='et':
    decay_type_1 = 0*np.ones(len(cols["pt_1"]),dtype=np.int32)
if args.channel=='tt':
    decay_type_1 = 2*np.ones(len(cols["pt_1"]),dtype=np.int32)
                             
N = min(args.nevts,len(cols['pt_1']))

verbosity = -1
delta = 1.0/1.15
reg_order = 6.0
mX = 125.10
widthX = 2.5
            
results = fastmtt_cpp.fastmtt_cpp(int(N),
                                  cols["pt_1"],
                                  cols["eta_1"],
                                  cols["phi_1"],
                                  cols["mass_1"],
                                  decay_type_1,
                                  cols["pt_2"],
                                  cols["eta_2"],
                                  cols["phi_2"],
                                  cols["mass_2"],
                                  decay_type_2,
                                  cols["puppimet"],
                                  cols["puppimetphi"],
                                  cols["puppimetcov00"],
                                  cols["puppimetcov01"],
                                  cols["puppimetcov11"],
                                  verbosity,
                                  delta,
                                  reg_order,
                                  mX,
                                  widthX)

mass = results['mass']
x1 = results['x1']
x2 = results['x2']
x1_BW = results['x1_BW']
x2_BW = results['x2_BW']
x1_cons = reults['x1_cons']
x2_cons = reults['x2_cons']

pt1 = np.where(x1>0.01,cols['pt_1']/x1,cols['pt_1'])
pt2 = np.where(x2>0.01,cols['pt_2']/x2,cols['pt_2'])
pt1_BW = np.where(x1_BW>0.01,cols['pt_1']/x1_BW,cols['pt_1'])
pt2_BW = np.where(x2_BW>0.01,cols['pt_2']/x2_BW,cols['pt_2'])
pt1_cons = np.where(x1_cons>0.01,cols['pt_1']/x1_cons,cols['pt_1'])
pt2_cons = np.where(x2_cons>0.01,cols['pt_2']/x2_cons,cols['pt_2'])

dpt1_nom = np.where(cols['genPart_pt_1']>0.001,cols['FastMTT_pt_1']/cols['genPart_pt_1'],0.0001) 
dpt1_BW_nom = np.where(cols['genPart_pt_1']>0.001,cols['FastMTT_pt_1_constraint']/cols['genPart_pt_1'],0.0001) 

dpt1 = np.where(cols['genPart_pt_1']>0.001,pt1/cols['genPart_pt_1'],0.0001) 
dpt1_BW = np.where(cols['genPart_pt_1']>0.001,pt1_BW/cols['genPart_pt_1'],0.0001) 
dpt1_cons = np.where(cols['genPart_pt_1']>0.001,pt1_cons/cols['genPart_pt_1'],0.0001) 
    
dpt2_nom = np.where(cols['genPart_pt_2']>0.001,cols['FastMTT_pt_2']/cols['genPart_pt_2'],0.0001) 
dpt2_BW_nom = np.where(cols['genPart_pt_2']>0.001,cols['FastMTT_pt_2_constraint']/cols['genPart_pt_2'],0.0001) 

dpt2 = np.where(cols['genPart_pt_2']>0.001,pt2/cols['genPart_pt_2'],0.0001) 
dpt2_BW = np.where(cols['genPart_pt_2']>0.001,pt2_BW/cols['genPart_pt_2'],0.0001) 
dpt2_cons = np.where(cols['genPart_pt_2']>0.001,pt2_cons/cols['genPart_pt_2'],0.0001) 
    
f = ROOT.TFile("%s.root"%(args.sample),"recreate")
f.cd('')
hist_mvis = ROOT.TH1D("mvis","mvis",40,0.,200.)
hist_mtt = ROOT.TH1D("mtt","mtt",60,0.,300.)
hist_mtt_nom = ROOT.TH1D("mtt_nom","mtt",60,0.,300.)

hist_x1 = ROOT.TH1D('x1','x1',51,-0.01,1.01)
hist_x2 = ROOT.TH1D('x2','x2',51,-0.01,1.01)
hist_x1_BW = ROOT.TH1D('x1_BW','x1',51,-0.01,1.01)
hist_x2_BW = ROOT.TH1D('x2_BW','x2',51,-0.01,1.01)
hist_x1_cons = ROOT.TH1D('x1_cons','x1',51,-0.01,1.01)
hist_x2_cons = ROOT.TH1D('x2_cons','x2',51,-0.01,1.01)

hist_dpt1_nom = ROOT.TH1D("dpt1_nom","dpt",40,0.,2.)
hist_dpt1_BW_nom = ROOT.TH1D("dpt1_BW_nom","dpt",40,0.,2.)
hist_dpt1 = ROOT.TH1D("dpt1","dpt",40,0.,2.)
hist_dpt1_BW = ROOT.TH1D("dpt1_BW","dpt",40,0.,2.)
hist_dpt1_cons = ROOT.TH1D("dpt1_cons","dpt",40,0.,2.)

hist_dpt2_nom = ROOT.TH1D("dpt2_nom","dpt",40,0.,2.)
hist_dpt2_BW_nom = ROOT.TH1D("dpt2_BW_nom","dpt",40,0.,2.)
hist_dpt2 = ROOT.TH1D("dpt2","dpt",40,0.,2.)
hist_dpt2_BW = ROOT.TH1D("dpt2_BW","dpt",40,0.,2.)
hist_dpt2_cons = ROOT.TH1D("dpt2_cons","dpt",40,0.,2.)

fill_hist(hist_mvis,cols["m_vis"])
fill_hist(hist_mtt,mass)
fill_hist(hist_mtt_nom,cols["FastMTT_mass"])

fill_hist(hist_x1,x1)
fill_hist(hist_x2,x2)
fill_hist(hist_x1_BW,x1_BW)
fill_hist(hist_x2_BW,x2_BW)
fill_hist(hist_x1_cons,x1_cons)
fill_hist(hist_x2_cons,x2_cons)
    
fill_hist(hist_dpt1,dpt1)
fill_hist(hist_dpt1_BW,dpt1_BW)
fill_hist(hist_dpt1_cons,dpt1_const)
fill_hist(hist_dpt1_nom,dpt1_nom)
fill_hist(hist_dpt1_BW_nom,dpt1_BW_nom)
    
fill_hist(hist_dpt2,dpt2)
fill_hist(hist_dpt2_BW,dpt2_BW)
fill_hist(hist_dpt2_cons,dpt2_const)
fill_hist(hist_dpt2_nom,dpt2_nom)
fill_hist(hist_dpt2_BW_nom,dpt2_BW_nom)

f.cd('')
hist_mvis.Write("mvis")
hist_mtt.Write("mtt")
hist_mtt_nom.Write("mtt_nom")

hist_x1.Write("x1")
hist_x2.Write("x2")
hist_x1_BW.Write("x1_BW")
hist_x2_BW.Write("x2_BW")
hist_x1_cons.Write("x1_cons")
hist_x2_cons.Write("x2_cons")

hist_dpt1.Write("dpt1")
hist_dpt1_BW.Write("dpt1_BW")
hist_dpt1_cons.Write("dpt1_cons")
hist_dpt1_nom.Write("dpt1_nom")
hist_dpt1_BW_nom.Write("dpt1_BW_nom")

hist_dpt2.Write("dpt2")
hist_dpt2_BW.Write("dpt2_BW")
hist_dpt2_cons.Write("dpt2_cons")
hist_dpt2_nom.Write("dpt2_nom")
hist_dpt2_BW_nom.Write("dpt2_BW_nom")

f.Close()
