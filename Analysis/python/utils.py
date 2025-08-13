import ROOT 
import math
from array import array
import numpy as np
import os

#############################
##### General settings ######
#############################

##################################
# constants
##################################
muon_mass = 0.105658
electron_mass = 0.000511
tau_mass = 1.777
pi_mass = 0.140
pi0_mass = 0.13498
rho_mass = 0.7755
ctau = 0.0087

##################################
# folder where tuples are stored #
##################################
tupleFolderMG='/eos/cms/store/group/phys_tau/ksavva/For_Aliaksei/files/testingzpt'
tupleFolderPOWHEG='/eos/cms/store/group/phys_tau/ksavva/For_Aliaksei/files/production_22_06_2025'
tupleFolderV1='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/old/oldTuples'
tupleFolderV2='/eos/cms/store/group/phys_tau/lrussell/forAliaksei'

#########################################
# folder with outputs (ROOT histograms) #
#########################################
outputFolder = '/afs/cern.ch/work/r/rasp/CPHiggs/Analysis'

periods = {'Run3_2022preEE' : ['Run3_2022'] ,
           'Run3_2022postEE' : ['Run3_2022EE'],
           'Run3_2023preBPix' : ['Run3_2023'],
           'Run3_2023postBPix' : ['Run3_2023BPix'],
           'Run3_2022' : ['Run3_2022','Run3_2022EE'],
           'Run3_2023' : ['Run3_2023','Run3_2023BPix'],
           'Run3' : ['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix'],
           }
    
eraLumi = {
    "Run3_2022preEE"    : 7980.4,
    "Run3_2022postEE"   : 26671.7,
    "Run3_2023preBPix"  : 17794.0,
    "Run3_2023postBPix" : 9451.0,
    "Run3_2022"  : 34652.1,
    "Run3_2023"  : 27245.0,
    "Run3"       : 61897.1,
}

###############
# MC samples ##
###############

zmm_powheg_samples = {
    'Run3_2022'    : ['DYto2Mu_MLL_50to120_powheg','DYto2Mu_MLL_120to200_powheg'],
    'Run3_2022EE'  : ['DYto2Mu_MLL_50to120_powheg','DYto2Mu_MLL_120to200_powheg'],
    'Run3_2023'    : ['DYto2Mu_MLL_50to120_powheg','DYto2Mu_MLL_120to200_powheg'],
    'Run3_2023BPix': ['DYto2Mu_MLL_50to120_powheg','DYto2Mu_MLL_120to200_powheg'],
}

zee_powheg_samples = {
    'Run3_2022'    : ['DYto2E_MLL_50to120_powheg','DYto2E_MLL_120to200_powheg'],
    'Run3_2022EE'  : ['DYto2E_MLL_50to120_powheg','DYto2E_MLL_120to200_powheg'],
    'Run3_2023'    : ['DYto2E_MLL_50to120_powheg','DYto2E_MLL_120to200_powheg'],
    'Run3_2023BPix': ['DYto2E_MLL_50to120_powheg','DYto2E_MLL_120to200_powheg'],
}

ztt_powheg_samples = {
    'Run3_2022'    : ['DYto2Tau_MLL_50to120_powheg','DYto2Tau_MLL_120to200_powheg'],
    'Run3_2022EE'  : ['DYto2Tau_MLL_50to120_powheg','DYto2Tau_MLL_120to200_powheg'],
    'Run3_2023'    : ['DYto2Tau_MLL_50to120_powheg','DYto2Tau_MLL_120to200_powheg'],
    'Run3_2023BPix': ['DYto2Tau_MLL_50to120_powheg','DYto2Tau_MLL_120to200_powheg'],
}

dy_incl_MG_samples = {
    'Run3_2022'    : ['DYto2L_M_50_madgraphMLM'],
    'Run3_2022EE'  : ['DYto2L_M_50_madgraphMLM'],
    'Run3_2023'    : ['DYto2L_M_50_madgraphMLM'],
    'Run3_2023BPix': ['DYto2L_M_50_madgraphMLM'],
}

dy_ext_MG_samples = {
    'Run3_2022'    : ['DYto2L_M_50_madgraphMLM_ext1'],
    'Run3_2022EE'  : ['DYto2L_M_50_madgraphMLM_ext1'],
    'Run3_2023'    : [],
    'Run3_2023BPix': [],
}

dy_1j_MG_samples = {
    'Run3_2022'    : ['DYto2L_M_50_1J_madgraphMLM'],
    'Run3_2022EE'  : ['DYto2L_M_50_1J_madgraphMLM'],
    'Run3_2023'    : ['DYto2L_M_50_1J_madgraphMLM'],
    'Run3_2023BPix': ['DYto2L_M_50_1J_madgraphMLM'],
}

dy_2j_MG_samples = {
    'Run3_2022'    : ['DYto2L_M_50_2J_madgraphMLM'],
    'Run3_2022EE'  : ['DYto2L_M_50_2J_madgraphMLM'],
    'Run3_2023'    : ['DYto2L_M_50_2J_madgraphMLM'],
    'Run3_2023BPix': ['DYto2L_M_50_2J_madgraphMLM'],
}

dy_3j_MG_samples = {
    'Run3_2022'    : ['DYto2L_M_50_3J_madgraphMLM'],
    'Run3_2022EE'  : ['DYto2L_M_50_3J_madgraphMLM'],
    'Run3_2023'    : ['DYto2L_M_50_3J_madgraphMLM'],
    'Run3_2023BPix': ['DYto2L_M_50_3J_madgraphMLM'],
}

dy_4j_MG_samples = {
    'Run3_2022'    : ['DYto2L_M_50_4J_madgraphMLM'],
    'Run3_2022EE'  : ['DYto2L_M_50_4J_madgraphMLM'],
    'Run3_2023'    : ['DYto2L_M_50_4J_madgraphMLM'],
    'Run3_2023BPix': ['DYto2L_M_50_4J_madgraphMLM'],
}


zll_0j_samples = {
    'Run3_2022'    : ['DYto2L_M_50_0J_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2L_M_50_0J_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2L_M_50_0J_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2L_M_50_0J_amcatnloFXFX'],
}

zll_1j_samples = {
    'Run3_2022'    : ['DYto2L_M_50_1J_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2L_M_50_1J_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2L_M_50_1J_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2L_M_50_1J_amcatnloFXFX'],
}

zll_2j_samples = {
    'Run3_2022'    : ['DYto2L_M_50_2J_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2L_M_50_2J_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2L_M_50_2J_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2L_M_50_2J_amcatnloFXFX'],
}

zll_incl_samples = {
    'Run3_2022'    : ['DYto2L_M_50_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2L_M_50_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2L_M_50_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2L_M_50_amcatnloFXFX'],
}

zll_ext_samples = {
    'Run3_2022'    : ['DYto2L_M_50_amcatnloFXFX_ext1'],
    'Run3_2022EE'  : ['DYto2L_M_50_amcatnloFXFX_ext1'],
    'Run3_2023'    : [],
    'Run3_2023BPix': [],
    
}

ztt_0j_samples = {
    'Run3_2022'    : ['DYto2Tau_MLL_50_0J_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2Tau_MLL_50_0J_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2Tau_MLL_50_0J_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2Tau_MLL_50_0J_amcatnloFXFX'],
}

ztt_1j_samples = {
    'Run3_2022'    : ['DYto2Tau_MLL_50_1J_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2Tau_MLL_50_1J_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2Tau_MLL_50_1J_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2Tau_MLL_50_1J_amcatnloFXFX'],
}

ztt_2j_samples = {
    'Run3_2022'    : ['DYto2Tau_MLL_50_2J_amcatnloFXFX'],
    'Run3_2022EE'  : ['DYto2Tau_MLL_50_2J_amcatnloFXFX'],
    'Run3_2023'    : ['DYto2Tau_MLL_50_2J_amcatnloFXFX'],
    'Run3_2023BPix': ['DYto2Tau_MLL_50_2J_amcatnloFXFX'],
}

top_samples = {
    'Run3_2022'    : ['TTto2L2Nu','TTto2L2Nu_ext1',
                      'TTtoLNu2Q','TTtoLNu2Q_ext1',
                      'TTto4Q','TTto4Q_ext1'],
    'Run3_2022EE'  : ['TTto2L2Nu','TTto2L2Nu_ext1',
                      'TTtoLNu2Q','TTtoLNu2Q_ext1',
                      'TTto4Q','TTto4Q_ext1'],
    'Run3_2023'    : ['TTto2L2Nu',
                      'TTtoLNu2Q',
                      'TTto4Q'],
    'Run3_2023BPix': ['TTto2L2Nu',
                      'TTtoLNu2Q',
                      'TTto4Q'],
}
               
vv_samples = {
    'Run3_2022'    : ['WW','WZ','ZZ',
                      'ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1',
                      'ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1',
                      'ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1',
                      'ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1'],
    'Run3_2022EE'  : ['WW','WZ','ZZ',
                      'ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1',
                      'ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1',
                      'ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1',
                      'ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1'],
    'Run3_2023'    : ['WW','WZ','ZZ',
                      'ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu',
                      'ST_tW_antitop_2L2Nu',
                      'ST_tW_top_LNu2Q',
                      'ST_tW_antitop_LNu2Q'],
    'Run3_2023BPix': ['WW','WZ','ZZ',
                      'ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu',
                      'ST_tW_antitop_2L2Nu',
                      'ST_tW_top_LNu2Q',
                      'ST_tW_antitop_LNu2Q'],
}

wjets_samples = {
    'Run3_2022'    : ['WtoLNu_madgraphMLM','WtoLNu_madgraphMLM_ext1',
                      'WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM',
                      'WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM'],
    'Run3_2022EE'  : ['WtoLNu_madgraphMLM','WtoLNu_madgraphMLM_ext1',
                      'WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM',
                      'WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM'],
    'Run3_2023'    : ['WtoLNu_madgraphMLM',
                      'WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM',
                      'WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM'],
    'Run3_2023BPix': ['WtoLNu_madgraphMLM',
                      'WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM',
                      'WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM'],
}

ggH_even_samples = ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay']
ggH_odd_samples = ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay']
ggH_maxmix_samples = ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay']
qqH_samples = ['VBFHToTauTau_UncorrelatedDecay_Filtered']
HWplus_samples = ['WplusHToTauTau_UncorrelatedDecay_Filtered']
HWminus_samples = ['WminusHToTauTau_UncorrelatedDecay_Filtered']
ZH_samples = ['ZHToTauTau_UncorrelatedDecay_Filtered']

################
# Data samples #
################

muonSamples = {'Run3_2022': ['SingleMuon_Run2022C','Muon_Run2022C','Muon_Run2022D'],
               'Run3_2022EE': ['Muon_Run2022E','Muon_Run2022F','Muon_Run2022G'],
               'Run3_2023': ['Muon0_Run2023C_v1','Muon0_Run2023C_v2','Muon0_Run2023C_v3','Muon0_Run2023C_v4','Muon1_Run2023C_v1','Muon1_Run2023C_v2','Muon1_Run2023C_v3','Muon1_Run2023C_v4'],
               'Run3_2023BPix': ['Muon0_Run2023D_v1','Muon0_Run2023D_v2','Muon1_Run2023D_v1','Muon1_Run2023D_v2'] 
               }


elecSamples = {'Run3_2022': ['EGamma_Run2022C','EGamma_Run2022D'],
               'Run3_2022EE': ['EGamma_Run2022E','EGamma_Run2022F','EGamma_Run2022G'],
               'Run3_2023': ['EGamma0_Run2023C_v1','EGamma0_Run2023C_v2','EGamma0_Run2023C_v3','EGamma0_Run2023C_v4','EGamma0_Run2023C_v1','EGamma0_Run2023C_v2','EGamma1_Run2023C_v3','EGamma1_Run2023C_v4'],
               'Run3_2023BPix': ['EGamma0_Run2023D_v1','EGamma0_Run2023D_v2','EGamma1_Run2023D_v1','EGamma1_Run2023D_v2']}

tauVsEleWPs = {
    'VVVLoose': "1",
    'VVLoose' : "2",
    'VLoose'  : "3",
    'Loose'   : "4",
    'Medium'  : "5",
    'Tight'   : "6",
    'VTight'  : "7",
    'VVTight' : "8"
}

tauVsEleIntWPs = {
    'VVVLoose': 1,
    'VVLoose' : 2,
    'VLoose'  : 3,
    'Loose'   : 4,
    'Medium'  : 5,
    'Tight'   : 6,
    'VTight'  : 7,
    'VVTight' : 8
}

tauVsMuWPs = {
    'VLoose'  : "1",
    'Loose'   : "2",
    'Medium'  : "3",
    'Tight'   : "4"
}

tauVsMuIntWPs = {
    'VLoose'  : 1,
    'Loose'   : 2,
    'Medium'  : 3,
    'Tight'   : 4
}

tauWPs = {
    'Loose': "4",
    'Medium': "5",
    'Tight': "6",
    'VTight': "7",
    'VVTight': "8"
}

tauIntWPs = {
    'Loose': 4,
    'Medium': 5,
    'Tight': 6,
    'VTight': 7,
    'VVTight': 8    
}

##################################
### Settings for SF measurements #
##################################

decayModeCuts = {
    '1prong'    : 'dm_2==0',
    '1prongPi0' : '(dm_2==1||dm_2==2)',
    '3prong'    : 'dm_2==10',
    '3prongPi0' : 'dm_2==11'
}

decayProngCuts = {
    '1prong' : '(dm_2==0||dm_2==1||dm_2==2)',
    '2prong' : '(dm_2==5||dm_2==6)',
    '3prong' : '(dm_2==10||dm_2==11)'
}

lib_histos = {
    'm_vis' : [240,0.,240.],
    'mt_1'  : [250,0.,250.],
    'met'   : [250,0.,250.],
    'pt_1'  : [200,0.,200.],
    'pt_2'  : [200,0.,200.],
    'eta_1' : [50,-2.5,2.5],
    'eta_2' : [50,-2.5,2.5],
    'ipsig_1': [200,0.,10.],
    'ipsig_2': [200,0.,10.],
    'n_jets' : [10,-0.5,9.5],
    'n_bjets' : [10,-0.5,9.5],
    'jpt_1' : [60,0.,300.],
    'jpt_2' : [60,0.,300.],
    'mjj' : [200,0.,2000.],
    'jdeta' : [80,0.,8.],
}

lib_datacards_histos = {
    'bdt_tautau' : [20,0.0,1.0],
    'bdt_jetFakes': [20,0.0,1.0],
    'bdt_signal': [20,0.0,1.0],
}
lib_datacards_1Dhistos = {
    'bdt_tautau' : [20,0.0,1.0],
    'bdt_jetFakes': [20,0.0,1.0],
    'bdt_signal': [20,0.0,1.0],
}
lib_datacards_2Dhistos = {
    'lep_pi'    : [20,0.0,1.0,40,0,360],
    'lep_rho'   : [20,0.0,1.0,40,0,360],
    'lep_a1_1pr': [20,0.0,1.0,40,0,360],
    'lep_a1_3pr': [20,0.0,1.0,40,0,360],
}

lib_phiCP_histos = {    
    'alpha_lep_pi': [90,0.,90.],
    'alpha_lep_rho': [90,0.,90.],
    'alpha_lep_a1': [90,0.,90.],

    'aco_lep_pi_plus': [8,0.,360.],
    'aco_lep_pi_minus': [8,0.,360.],
    'aco_lep_pi': [8,0.,360.],
    'aco_lep_piIP': [8,0.,360.],
    'aco_lepIP_pi': [8,0.,360.],
    'aco_lepIP_piIP': [8,0.,360.],    
    
    'aco_lep_rho_plus': [8,0.,360.],
    'aco_lep_rho_minus': [8,0.,360.],
    'aco_lep_rho': [8,0.,360.],
    'aco_lepIP_rho': [8,0.,360.],

    'aco_lep_rhoECut': [8,0.,360.],
    'aco_lep_rhoGen': [8,0.,360.],
    'aco_lep_rhoReco': [8,0.,360.],
    'aco_lep_rhoRecoECut': [8,0.,360.],
    'aco_lep_rhoCollinear': [8,0.,360],
    'aco_lep_rhoRecoGen': [8,0.,360.],
    'aco_lep_rhoRecoGenECut': [8,0.,360.],

    'aco_lep_rhoRecoIP1p0': [8,0.,360.],
    'aco_lep_rhoRecoGenIP1p0': [8,0.,360.],
    'aco_lep_rhoRecoIP1p2': [8,0.,360.],
    'aco_lep_rhoRecoGenIP1p2': [8,0.,360.],
    'aco_lep_rhoRecoIP1p5': [8,0.,360.],
    'aco_lep_rhoRecoGenIP1p5': [8,0.,360.],

    'aco_lep_rhoRecoIP1p0ECut': [8,0.,360.],
    'aco_lep_rhoRecoGenIP1p0ECut': [8,0.,360.],
    'aco_lep_rhoRecoIP1p2ECut': [8,0.,360.],
    'aco_lep_rhoRecoGenIP1p2ECut': [8,0.,360.],
    'aco_lep_rhoRecoIP1p5ECut': [8,0.,360.],
    'aco_lep_rhoRecoGenIP1p5ECut': [8,0.,360.],

    'aco_lep_a1_plus': [8,0.,360.],
    'aco_lep_a1_minus': [8,0.,360.],
    'aco_lep_a1': [8,0.,360.],
    'aco_lepIP_a1': [8,0.,360.],
    'aco_lep_a1DP': [8,0.,360.],
    'aco_lep_a1PVGen': [8,0.,360.],
    'aco_lepIP_a1PVGen': [8,0.,360.],
    'aco_lep_a1PVDESY': [8,0.,360.],
    'aco_lepIP_a1PVDESY': [8,0.,360.],
    'aco_lep_a1PVIC': [8,0.,360.],
    

#    'aco_pi_pi': [8,0.,360.],
#    'aco_pi_rho': [8,0.,360.],
#    'aco_pi_rhoPV': [8,0.,360.],
#    'aco_pi_a1': [8,0.,360.],
#    'aco_pi_a1DP': [8,0.,360.],
#    'aco_pi_a1PV': [8,0.,360.],
    
#    'aco_rho_rho': [8,0.,360.],
#    'aco_rhoPV_rhoPV': [8,0.,360.],
#    'aco_rho_a1': [8,0.,360.],
#    'aco_rho_a1DP': [8,0.,360.],
#    'aco_rho_a1PV': [8,0.,360.],
#    'aco_a1_a1': [8,0.,360.]
    
}

# (pt,eta) binning for scale factors of IP significance cut 
ptbins = {
    'ee' : [30.,35.,40.,50.,60.,80.,110.,    150.],
    'mm' : [25.,30.,35.,40.,50.,60.,80.,110.,150.],
    'et' : [24.,30.,35.,40.,50.,80.,         150.],
    'mt' : [21.,25.,30.,35.,40.,50.,80.,     150.],
}

etabins = {
    'ee' : [0.,1.0,1.6,2.1],
    'mm' : [0.,1.0,1.6,2.4],
    'et' : [0.,1.0,1.6,2.1],
    'mt' : [0.,1.0,1.6,2.4]
}    
    

#######################################
# Creating shape systematic templates #
#######################################
def ComputeSystematics(h_central, h_sys, name):
    h_up = h_central.Clone(name+"Up")
    h_down = h_central.Clone(name+"Down")
    nbins = h_central.GetNbinsX()
    for i in range(1,nbins+1):
        x_up = h_sys.GetBinContent(i)
        x_central = h_central.GetBinContent(i)
        x_down = max(0,2.0*x_central-x_up)
        h_up.SetBinContent(i,x_up)
        h_down.SetBinContent(i,x_down)
    return h_up, h_down


##############################
# Utilities manipulating     #
# with arrays and histograms #
##############################

def createBins(nbins,xmin,xmax):
    binwidth = (xmax-xmin)/float(nbins)
    bins = []
    for i in range(0,nbins+1):
        xb = xmin + float(i)*binwidth
        bins.append(xb)
    return bins

def removeNegativeBins(hist):
    nbins = hist.GetNbinsX()
    for i in range(1,nbins+1):
        x = hist.GetBinContent(i)
        if x<0:
            hist.SetBinError(i,0.)
            hist.SetBinContent(i,0.)
            

def zeroBinErrors(hist):
    nbins = hist.GetNbinsX()
    for i in range(1,nbins+1):
        hist.SetBinError(i,0.)

def createUnitHisto(hist,histName):
    nbins = hist.GetNbinsX()
    unitHist = hist.Clone(histName)
    for i in range(1,nbins+1):
        x = hist.GetBinContent(i)
        e = hist.GetBinError(i)
        if x>0:
            rat = e/x
            unitHist.SetBinContent(i,1.)
            unitHist.SetBinError(i,rat)
    return unitHist

def dividePassProbe(passHist,failHist,histName):
    nbins = passHist.GetNbinsX()
    hist = passHist.Clone(histName)
    for i in range(1,nbins+1):
        xpass = passHist.GetBinContent(i)
        epass = passHist.GetBinError(i)
        xfail = failHist.GetBinContent(i)
        efail = failHist.GetBinError(i)
        xprobe = xpass+xfail
        ratio = 1
        eratio = 0
        if xprobe>1e-4:
            ratio = xpass/xprobe
            dpass = xfail*epass/(xprobe*xprobe)
            dfail = xpass*efail/(xprobe*xprobe)
            eratio = math.sqrt(dpass*dpass+dfail*dfail)
        hist.SetBinContent(i,ratio)
        hist.SetBinError(i,eratio)
    return hist

def divideHistos(numHist,denHist,histName):
    nbins = numHist.GetNbinsX()
    hist = numHist.Clone(histName)
    for i in range(1,nbins+1):
        xNum = numHist.GetBinContent(i)
        eNum = numHist.GetBinError(i)
        xDen = denHist.GetBinContent(i)
        eDen = denHist.GetBinError(i)
        ratio = 1
        eratio = 0
        if xNum>1e-7:
            ratio = xNum/xDen
            rNum = eNum/xNum
            rDen = eDen/xDen
            rratio = math.sqrt(rNum*rNum+rDen*rDen)
            eratio = rratio * ratio
        else:
            ratio = 0.5*eNum/xDen
            eratio = ratio
            
        hist.SetBinContent(i,ratio)
        hist.SetBinError(i,eratio)
    return hist

def histoRatio(numHist,denHist,histName):
    nbins = numHist.GetNbinsX()
    hist = numHist.Clone(histName)
    for i in range(1,nbins+1):
        xNum = numHist.GetBinContent(i)
        eNum = numHist.GetBinError(i)
        xDen = denHist.GetBinContent(i)
        ratio = 1
        eratio = 0
        if xNum>1e-7 and xDen>1e-7:
            ratio = xNum/xDen
            eratio = eNum/xDen
        hist.SetBinContent(i,ratio)
        hist.SetBinError(i,eratio)
    return hist

def rebinHisto(hist,bins,suffix):
    nbins = hist.GetNbinsX()
    newbins = len(bins)-1
    name = hist.GetName()+"_"+suffix
    newhist = ROOT.TH1D(name,"",newbins,array('d',list(bins)))
    for ib in range(1,nbins+1):
        centre = hist.GetBinCenter(ib)
        bin_id = newhist.FindBin(centre)
        xbin = hist.GetBinContent(ib)
        ebin = hist.GetBinError(ib)
        xnew = newhist.GetBinContent(bin_id)
        enew = newhist.GetBinError(bin_id)
        x_update = xbin + xnew;
        e_update = math.sqrt(ebin*ebin + enew*enew);
        newhist.SetBinContent(bin_id,x_update)
        newhist.SetBinError(bin_id,e_update)
    return newhist

def copyHist(inputHist,outputHist):
    nbins = inputHist.GetNbinsX()
    for ib in range(1,nbins+1):
        outputHist.SetBinContent(ib,inputHist.GetBinContent(ib))
        outputHist.SetBinError(ib,inputHist.GetBinError(ib))

# extracting Tag-and-Probe histos from ROOT file created by RunSelection.py macro
def extractTagProbeHistos(f,bins,generator,era,channel,isSecond):
    var = 'm_vis'
    if isSecond:
        var += '_2'
    is2022 = era=='Run3_2022preEE' or era=='Run3_2022postEE' or era=='Run3_2022'
    samples = ['data','dy']
    if channel=='mt' or channel=='et':
        samples.append('top')
        samples.append('vv')
        samples.append('wjets')
    dy_samples = ['zll_0j','zll_1j','zll_2j','ztt_0j','ztt_1j','ztt_2j']
    if is2022:
        dy_samples.append('zll_ext')
    dy_base = 'zll_incl'
    if generator=='MG':
        dy_samples = ['dy_1j','dy_2j','dy_3j','dy_4j']
        if is2022:
            dy_samples.append('dy_ext')
        dy_base = 'dy_incl'
    if generator=='powheg':
        dy_samples = []
        dy_base = 'zll_powheg'
    sign_labels = ['os','ss']
    iso_labels = ['iso','antiiso']
    typ_labels = ['lep','tau','had','all']
    region_labels = ['pass','fail','incl']
    hists = {}
    histPtBins = f.Get('ptBins')
    histEtaBins = f.Get('etaBins')
    nbinsPt = histPtBins.GetNbinsX()
    nbinsEta = histEtaBins.GetNbinsX()
    for sample in samples:
        for sign in sign_labels:
            for iso in iso_labels:
                for typ in typ_labels:
                    for binPt in range(1,nbinsPt+1):
                        for binEta in range(1,nbinsEta+1):
                            label = '%1i_%1i'%(binPt,binEta)
                            for region in region_labels:
                                if sample=='dy':
                                    name = 'dy_m_vis_%s_%s_%s_%s_%s'%(region,label,sign,iso,typ)
                                    nameInput = '%s_%s_%s_%s_%s_%s_%s'%(dy_base,var,region,label,sign,iso,typ)
                                    histo = f.Get(nameInput).Clone(name)
                                    for dy_sample in dy_samples:
                                        dy_name = '%s_%s_%s_%s_%s_%s_%s'%(dy_sample,var,region,label,sign,iso,typ)
                                        dy_histo = f.Get(dy_name)
                                        histo.Add(histo,dy_histo,1.,1.)
                                    hists[name] = rebinHisto(histo,bins,'rebinned')
                                else:
                                    name = '%s_m_vis_%s_%s_%s_%s_%s'%(sample,region,label,sign,iso,typ)
                                    nameInput = '%s_%s_%s_%s_%s_%s_%s'%(sample,var,region,label,sign,iso,typ)
                                    histo = f.Get(nameInput).Clone(name)
                                    hists[name] = rebinHisto(histo,bins,'rebinned')
#                                for unc in unc_labels:
#                                    name = '%s_m_vis_%s_%s_%s_%s_%s_%s'%(sample,unc,region,label,sign,iso,typ)
#                                    histo = f.Get(name)
#                                    hists[name] = rebinHisto(histo,bins,'rebinned')
    return hists

# extracting histos from ROOT file created by RunSelection.py 
def extractHistos(f,var,bins,generator,era,channel):
    is2022 = era=='Run3_2022preEE' or era=='Run3_2022postEE' or era=='Run3_2022'
    samples = ['data','dy','top','vv','wjets']
    dy_samples = []
    dy_base = ''
    if generator=='amcatnlo':
        dy_samples.extend(['zll_0j','zll_1j','zll_2j','ztt_0j','ztt_1j','ztt_2j'])
        if is2022:
            dy_samples.append('zll_ext')
        dy_base = 'zll_incl'
    if generator=='MG':
        dy_samples.extend(['dy_1j','dy_2j','dy_3j','dy_4j'])
        if is2022:
            dy_samples.append('dy_ext')
        dy_base = 'dy_incl'
    if generator=='powheg':
        if channel=='mt' or channel=='et':
            dy_samples.append('ztt_powheg')
        dy_base = 'zll_powheg'
    sign_labels = ['os','ss']
    iso_labels = ['iso','antiiso']
    typ_labels = ['lep','tau','had','all']
    hists = {}
    for sample in samples:
        for sign in sign_labels:
            for iso in iso_labels:
                for typ in typ_labels:
                    if sample=='dy':
                        name = 'dy_%s_%s_%s_%s'%(var,sign,iso,typ)
                        nameInput = '%s_%s_%s_%s_%s'%(dy_base,var,sign,iso,typ)
                        histo = f.Get(nameInput).Clone(name)
                        for dy_sample in dy_samples:
                            dy_name = '%s_%s_%s_%s_%s'%(dy_sample,var,sign,iso,typ)
                            dy_histo = f.Get(dy_name)
                            histo.Add(histo,dy_histo,1.,1.)
                        hists[name] = rebinHisto(histo,bins,'rebinned')
                    else:
                        name = '%s_%s_%s_%s_%s'%(sample,var,sign,iso,typ)
                        histo = f.Get(name)
                        hists[name] = rebinHisto(histo,bins,'rebinned')
    return hists

