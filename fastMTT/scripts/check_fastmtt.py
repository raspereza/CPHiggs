#! /usr/bin/env python3
# Author: Alexei Raspereza (June 2025)
# Code to validate FastMTT programme

import awkward as ak
import numpy as np
from fastmtt_cpp import fastmtt_cpp
import ROOT 

def fill_hist(hist, array):
    [hist.Fill(x) for x in array]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-isDY','--sample',dest='sample',default='GluGluHToTauTauUncorrDecays_M125',choices=["GluGluHToTauTauUncorrDecays_M125","DYJetsToLL_M-50",""])
args = parser.parse_args()

dirname_v1=''
dirname_v2=

#
#'/eos/cms/store/group/phys_tau/lrussell/forAliaksei/Run3_2022EE/mt/GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay/nominal/merged.root'

df = ROOT.RDataFrame("ntuple",filename)

cols = df.AsNumpy(["pt_1","pt_2","eta_1","eta_2","phi_1","phi_2","m_1","m_2","puppimet","puppimetphi","puppimetcov00","puppimetcov01","puppimetcov11","m_vis","gen_Higgs_pt","gen_Higgs_eta","gen_Higgs_phi"])

print('Length of %: \n',len(cols["pt_1"]))

decay_type_1 = 1*np.ones(len(cols["pt_1"]),dtype=np.int32)
decay_type_2 = 2*np.ones(len(cols["pt_1"]),dtype=np.int32)
N = 100000

verbosity = -1
delta = 1.0/1.15
reg_order = 6.0
mX = 125.10
GammaX = 0.004

results = fastmtt(int(N),
                  cols["pt_1"],
                  cols["eta_1"],
                  cols["phi_1"],
                  cols["m_1"],
                  decay_type_1,
                  cols["pt_2"],
                  cols["eta_2"],
                  cols["phi_2"],
                  cols["m_2"],
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
                  GammaX)

eta_1 = cols["eta_1"]
eta_2 = cols["eta_2"]

phi_1 = cols["phi_1"]
phi_2 = cols["phi_2"]

x1 = np.where(results["x1"]>0.01,results["x1"],0.5)
x2 = np.where(results["x2"]>0.01,results["x2"],0.5)

pt_1 = np.where(results["x1"]>0.01,cols["pt_1"]/x1,cols["pt_1"])
pt_2 = np.where(results["x2"]>0.01,cols["pt_2"]/x2,cols["pt_2"])

x1_cons = np.where(results["x1_cons"]>0.01,results["x1_cons"],0.5)
x2_cons = np.where(results["x2_cons"]>0.01,results["x2_cons"],0.5)

pt_1_cons = np.where(results["x1_cons"]>0.01,cols["pt_1"]/x1_cons,cols["pt_2"])
pt_2_cons = np.where(results["x2_cons"]>0.01,cols["pt_2"]/x2_cons,cols["pt_2"])

mass = np.where(results["x1"]>0.01,results["mass"],-9999999.)
mass_cons = np.where(results["x1"]>0.01,results["mass_cons"],-9999999.)

px_H = pt_1 * np.cos(phi_1) + pt_2 * np.cos(phi_2)
py_H = pt_1 * np.sin(phi_1) + pt_2 * np.sin(phi_2)
pz_H = pt_1 * np.sinh(eta_1) + pt_2 * np.sinh(eta_2)
pt_H = np.sqrt(px_H*px_H+py_H*py_H)

px_H_cons = pt_1_cons * np.cos(phi_1) + pt_2_cons * np.cos(phi_2)
py_H_cons = pt_1_cons * np.sin(phi_1) + pt_2_cons * np.sin(phi_2)
pz_H_cons = pt_1_cons * np.sinh(eta_1) + pt_2_cons * np.sinh(eta_2)
pt_H_cons = np.sqrt(px_H_cons*px_H_cons+py_H_cons*py_H_cons)

dpt = pt_H - cols["gen_Higgs_pt"]
dpt_BW = pt_H_cons - cols["gen_Higgs_pt"]

f = ROOT.TFile("%s.root"%(args.sample),"recreate")
hist_mvis = ROOT.TH1D("mvis","mvis",40,0.,200.)
hist_mtt = ROOT.TH1D("mtt","mtt",60,0.,300.)
hist_dpt = ROOT.TH1D("dpt","dpt",100,-100.,100.)
hist_dpt_BW = ROOT.TH1D("dpt_cons","dpt_cons",100,-100.,100.)
fill_hist(hist_mvis,cols["m_vis"])
fill_hist(hist_mtt,mass)
fill_hist(hist_dpt,dpt)
fill_hist(hist_dpt_BW,dpt_BW)
f.cd("")
hist_mvis.Write("mvis")
hist_mtt.Write("mtt")
hist_dpt.Write("dpt")
hist_dpt_BW.Write("dpt_BW")
f.Close()
