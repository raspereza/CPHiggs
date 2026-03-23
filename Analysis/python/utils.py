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
massH = 125.10

##################################
# folder where tuples are stored #
##################################
#tupleFolderPhiCP='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/oldTuples/old0407'
#tupleFolderPhiCP='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/oldTuples/old0808'
tupleFolderPhiCP='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/CPSignalStudies'
tupleFolderV2='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/ForFakeFactors'

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
    "Run3_2023preBPix"  : 18063.0,
    "Run3_2023postBPix" : 9693.0,
    "Run3_2022"  : 34652.1,
    "Run3_2023"  : 27245.0,
    "Run3"       : 61897.1,
}

eraLumiOld = {
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

top_2l2v_samples = {
    'Run3_2022'    : ['TTto2L2Nu'],
    'Run3_2022EE'  : ['TTto2L2Nu'],
    'Run3_2023'    : ['TTto2L2Nu'],
    'Run3_2023BPix': ['TTto2L2Nu'],
}

top_2l2v_ext_samples = {
    'Run3_2022'    : ['TTto2L2Nu_ext1'],
    'Run3_2022EE'  : ['TTto2L2Nu_ext1'],
    'Run3_2023'    : [],
    'Run3_2023BPix': [],
}

top_lv2q_samples = {
    'Run3_2022'    : ['TTtoLNu2Q'],
    'Run3_2022EE'  : ['TTtoLNu2Q'],
    'Run3_2023'    : ['TTtoLNu2Q'],
    'Run3_2023BPix': ['TTtoLNu2Q'],
}

top_lv2q_ext_samples = {
    'Run3_2022'    : ['TTtoLNu2Q_ext1'],
    'Run3_2022EE'  : ['TTtoLNu2Q_ext1'],
    'Run3_2023'    : [],
    'Run3_2023BPix': [],
}

vv_samples = {
    'Run3_2022'    : ['WW','WZ','ZZ'],
    'Run3_2022EE'  : ['WW','WZ','ZZ'],
    'Run3_2023'    : ['WW','WZ','ZZ'],
    'Run3_2023BPix': ['WW','WZ','ZZ'],
}

st_samples = {
    'Run3_2022'    : ['ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1',
                      'ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1',
                      'ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1',
                      'ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1'],
    'Run3_2022EE'  : ['ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1',
                      'ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1',
                      'ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1',
                      'ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1'],
    'Run3_2023'    : ['ST_t_channel_top_4f_InclusiveDecays',
                      'ST_t_channel_antitop_4f_InclusiveDecays',
                      'ST_tW_top_2L2Nu',
                      'ST_tW_antitop_2L2Nu',
                      'ST_tW_top_LNu2Q',
                      'ST_tW_antitop_LNu2Q'],
    'Run3_2023BPix': ['ST_t_channel_top_4f_InclusiveDecays',
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

ggH_even_samples = {
    'Run3_2022' : ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
    'Run3_2022EE' : ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
    'Run3_2023' : ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
    'Run3_2023BPix' : ['GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'],
}
ggH_odd_samples = {
    'Run3_2022' : ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
    'Run3_2022EE' : ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
    'Run3_2023' : ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
    'Run3_2023BPix' : ['GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay'],
}    

ggH_maxmix_samples = {
    'Run3_2022' : ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
    'Run3_2022EE' : ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
    'Run3_2023' : ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
    'Run3_2023BPix' : ['GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay'],
}

qqH_samples = {
    'Run3_2022' : ['VBFHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2022EE' : ['VBFHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023' : ['VBFHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023BPix' : ['VBFHToTauTau_UncorrelatedDecay_Filtered'],
}

HWplus_samples = {
    'Run3_2022' : ['WplusHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2022EE' : ['WplusHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023' : ['WplusHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023BPix' : ['WplusHToTauTau_UncorrelatedDecay_Filtered'],
}

HWminus_samples = {
    'Run3_2022' : ['WminusHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2022EE' : ['WminusHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023' : ['WminusHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023BPix' : ['WminusHToTauTau_UncorrelatedDecay_Filtered'],
}

ZH_samples = {
    'Run3_2022' : ['ZHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2022EE' : ['ZHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023' : ['ZHToTauTau_UncorrelatedDecay_Filtered'],
    'Run3_2023BPix' : ['ZHToTauTau_UncorrelatedDecay_Filtered'],
}


samplesDict = {
    'zll_incl'     : zll_incl_samples,
    'zll_ext'      : zll_ext_samples,
    'zll_0j'       : zll_0j_samples,
    'zll_1j'       : zll_1j_samples,
    'zll_2j'       : zll_2j_samples,
    'ztt_0j'       : ztt_0j_samples,
    'ztt_1j'       : ztt_1j_samples,
    'ztt_2j'       : ztt_2j_samples,
    'top_2l2v'     : top_2l2v_samples,
    'top_2l2v_ext' : top_2l2v_ext_samples,
    'top_lv2q'     : top_lv2q_samples,
    'top_lv2q_ext' : top_lv2q_ext_samples,
    'vv'           : vv_samples,
    'st'           : st_samples,
    'wjets'        : wjets_samples,
    'ggH_sm'       : ggH_even_samples,
    'ggH_ps'       : ggH_odd_samples,
    'ggH_mm'       : ggH_maxmix_samples,
    'qqH'          : qqH_samples,
    'HWplus'       : HWplus_samples, 
    'HWminus'      : HWminus_samples,
    'ZH'           : ZH_samples,
}

samplesNameDict = {
    'data'          : 'data',
    'zll_incl'      : 'zll',
    'zll_ext'       : 'zll',
    'zll_0j'        : 'zll',
    'zll_1j'        : 'zll',
    'zll_2j'        : 'zll',
    'ztt_0j'        : 'ztt',
    'ztt_1j'        : 'ztt',
    'ztt_2j'        : 'ztt',
    'top_2l2v'      : 'top',
    'top_2l2v_ext'  : 'top',
    'top_lv2q'      : 'top',
    'top_lv2q_ext'  : 'top',
    'vv'            : 'vv',
    'st'            : 'vv',
    'wjets'         : 'wjets',
    'ggH_sm'        : 'ggH_sm',
    'ggH_ps'        : 'ggH_ps',
    'ggH_mm'        : 'ggH_mm',
    'qqH'           : 'qqH',
    'HWplus'        : 'HWplus',
    'HWminus'       : 'HWminus',
    'ZH'            : 'ZH',
}

signal_samples = ['ggH_sm','ggH_ps', 'ggH_mm', 'qqH', 'HWplus', 'HWminus', 'ZH']
bkg_samples = ['zll','ztt','top','vv','wjets']
samples = ['data','zll','ztt','top','vv','wjets']

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
               'Run3_2023': ['EGamma0_Run2023C_v1','EGamma0_Run2023C_v2','EGamma0_Run2023C_v3','EGamma0_Run2023C_v4','EGamma1_Run2023C_v1','EGamma1_Run2023C_v2','EGamma1_Run2023C_v3','EGamma1_Run2023C_v4'],
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
    'm_vis' :   [240,0.,240.],
    'm_FastMTT' : [400,0.,400.],
    'pt_tt' :   [400,0.,400.],
    'pt_FastMTT' : [400,0.,400.],
    'mt_1'  :   [250,0.,250.],
    'met'   :   [250,0.,250.],
    'pt_1'  :  [200,0.,200.],
    'pt_2'  :  [200,0.,200.],
    'eta_1' :  [50,-2.5,2.5],
    'eta_2' :  [50,-2.5,2.5],
    'dm_2'  :  [12,-0.5,11.5],
    'ipsig_1': [200,0.,10.],
    'ipsig_2': [200,0.,10.],
    'CMetQCD' : [100,-5.0,5.0],
    'CMetW'   : [100,-5.0,5.0],
    'n_jets' : [10,-0.5,9.5],
    'n_bjets': [10,-0.5,9.5],
    'jpt_1' :   [60,0.,300.],
    'jpt_2' : [60,0.,300.],
    'jeta_1': [60,-6.0,6.0],
    'jeta_2': [60,-6.0,6.0],
    'mjj'   : [200,0.,2000.],
    'dR'    : [100,0.,10.],
    'jdeta' : [80,0.,8.],
    'aco_lep_pi': [360, 0., 360.],
    'aco_lep_rho': [360, 0., 360. ],
    'aco_lep_a1_1pr': [360, 0., 360.],
    'aco_lep_a1_3pr': [360, 0., 360.],
    'bdt_ditau'  : [100,0.,1.],
    'bdt_fakes'  : [100,0.,1.],
    'bdt_signal' : [100,0.,1.],
}

sign_labels = ['os','ss']
iso_labels = ['iso','antiiso']
tautype_labels = ['lepton','fake','all']
type_labels = ['tau','lep','had','all']
tauid_labels = ['inverted','nominal']
ff_labels = ['ar','qcd','wj','mc_top','mc_wj','qcd_looseIso']
njets_labels = ['njets0','njets1','njets2']
dm_labels = ['pi','rho','a1_1pr','a1_3pr']
eta_labels = ['barrel','endcap','all']
mt_labels = ['lowmt','highmt']
cp_hypotheses = ['sm','ps','mm']
region_labels = ['lowmt_os_iso',
                 'lowmt_ss_iso',
                 'qcd_closure_ff',
                 'wj_ff',
                 'qcd_ff',
                 'qcd_looseIso_ff']

lib_datacards_1Dhistos = ['bdt_ditau','bdt_fakes','bdt_signal']

lib_phiCP_histos = {
    
    'aco_lep_pi': [8,0.,360.],
    'aco_lep_rho': [8,0.,360.],
    
    'aco_lep_a1': [8,0.,360.],
    'aco_lep_a1_chi5': [8,0.,360.],
    'aco_lep_a1_chi10': [8,0.,360.],
    'aco_lep_a1_chi20': [8,0.,360.],
    
    'aco_lep_a1_KinFit': [8,0.,360.],
    'aco_lep_a1_chi5_KinFit': [8,0.,360.],    
    'aco_lep_a1_chi10_KinFit': [8,0.,360.],    
    'aco_lep_a1_chi20_KinFit': [8,0.,360.],    

    'chi2_lep_a1' : [1000,0.,100.],
    
    'aco_lep_a1_GenPV': [8,0.,360.],
    'aco_lep_a1_GenTheta': [8,0.,360.],

    'dPt': [500,0.,5.],
    'dPt_KinFit': [500,0.,5.],

    'dTheta': [500,-0.1,0.1],
    'dTheta_KinFit': [500,-0.1,0.1],
    'dTheta_KinFitSV': [500,-0.1,0.1],

    'dAlpha': [200,0.,0.1],
    'dAlpha_KinFit': [200,0.,0.1],
    'dAlpha_KinFitSV': [200,0.,0.1],
    
    'theta_GJ': [200,0.,0.2],
    'theta_GJ_max': [200,0.,0.2],
    
    'aco_pi_pi': [8,0.,360.],
    'aco_pi_rho': [8,0.,360.],
    'aco_rho_rho': [8,0.,360.],

    'aco_pi_a1': [8,0.,360.],
    'aco_rho_a1': [8,0.,360.],
    'aco_a1_a1': [8,0.,360.],

    'aco_pi_a1_chi5': [8,0.,360.],
    'aco_rho_a1_chi5': [8,0.,360.],
    'aco_a1_a1_chi5': [8,0.,360.],
    
    'aco_pi_a1_chi10': [8,0.,360.],
    'aco_rho_a1_chi10': [8,0.,360.],
    'aco_a1_a1_chi10': [8,0.,360.],

    'aco_pi_a1_chi20': [8,0.,360.],
    'aco_rho_a1_chi20': [8,0.,360.],
    'aco_a1_a1_chi20': [8,0.,360.],
    
    'aco_pi_a1_KinFit': [8,0.,360.],
    'aco_rho_a1_KinFit': [8,0.,360.],
    'aco_a1_a1_KinFit': [8,0.,360.],

    'aco_pi_a1_chi5_KinFit': [8,0.,360.],
    'aco_rho_a1_chi5_KinFit': [8,0.,360.],
    'aco_a1_a1_chi5_KinFit': [8,0.,360.],

    'aco_pi_a1_chi10_KinFit': [8,0.,360.],
    'aco_rho_a1_chi10_KinFit': [8,0.,360.],
    'aco_a1_a1_chi10_KinFit': [8,0.,360.],

    'aco_pi_a1_chi20_KinFit': [8,0.,360.],
    'aco_rho_a1_chi20_KinFit': [8,0.,360.],
    'aco_a1_a1_chi20_KinFit': [8,0.,360.],
    
    'aco_pi_a1_GenPV': [8,0.,360.],
    'aco_rho_a1_GenPV': [8,0.,360.],
    'aco_a1_a1_GenPV': [8,0.,360.],

    'aco_pi_a1_GenTheta': [8,0.,360.],
    'aco_rho_a1_GenTheta': [8,0.,360.],
    'aco_a1_a1_GenTheta': [8,0.,360.],

    'mass': [300, 0., 300.],
    'mass_KinFit': [300,0.,300.],
    
    'mass_a1_a1': [300, 0., 300.],
    'mass_a1_a1_KinFit': [300,0.,300.],

    'chi2_pi_a1'  : [1000,0.,100.],
    'chi2_rho_a1' : [1000,0.,100.],
    'chi2_a1_a1'  : [1000,0.,100.],
    
    
}

nbins_bdt_signal = {
    'tt' : 5,
    'mt' : 5,
    'et' : 3,
}

# (pt,eta) binning for scale factors of IP significance cut 
ptbins = {
    'ee' : [30.,35.,40.,50.,60.,80.,110.,    150.],
    'mm' : [25.,30.,35.,40.,50.,60.,80.,110.,150.],
    'et' : [31.,35.,40.,50.,80.,             150.],
    'mt' : [25.,30.,35.,40.,50.,80.,         150.],
    'tt' : [25.,30.,35.,40.,50.,80.,         150.],
}

jtau_ptbins = {
    'mt' : [20.,25.,30.,40.,50.,60.,80.,100.,150.],
    'et' : [20.,25.,30.,40.,50.,60.,80.,100.,150.],
    'tt' : [20.,25.,30.,40.,50.,60.,80.,100.,150.],
}

etabins = {
    'ee' : [0.,1.0,1.6,2.1],
    'mm' : [0.,1.0,1.6,2.4],
    'et' : [0.,1.0,1.6,2.1],
    'mt' : [0.,1.0,1.6,2.4],
    'tt' : [0.,1.0,1.6,2.4],
}    

XTitle = {
    'mt': {
        'mt_1'  : "m_{T} (GeV)",
        'pt_1'  : "muon p_{T} (GeV)",
        'eta_1' : "muon #eta",
        'pt_2'  : "tau p_{T} (GeV)",
        'eta_2' : "tau #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{vis} (GeV)",
        'ipsig_1': "muon IP sig",
        'ipsig_2': "tau IP sig",
        'n_jets': "number of jets",
        'n_bjets': "number of b-jets",
        'jpt_1': "leading jet p_{T} (GeV)",
        'jpt_2': "trailing jet p_{T} (GeV)",
        'mjj': "dijet mass (GeV)",
        'jdeta': "#Delta#eta(j,j)",
        'CMetQCD' : "C_{MET}^{QCD}",
        'CMetW' : "C_{MET}^{W}",
        'dR' : '#DeltaR(#mu,#tau)',
        'bdt_signal' : 'BDT_{sig}',
        'bdt_ditau' : 'BDT_{#tau#tau}',
        'bdt_fakes' : 'BDT_{fakes}',
        'aco_lep_pi': '#phi_{CP} (#pi)',
        'aco_lep_rho': '#phi_{CP} (#rho)',
        'aco_lep_a1_1pr': '#phi_{CP} (a_{1}(1-prong))',
        'aco_lep_a1_3pr': '#phi_{CP} (a_{1}(3-prong))',
    },
    'et': {
        'mt_1'  : "m_{T} (GeV)",
        'pt_1'  : "electron p_{T} (GeV)",
        'eta_1' : "electron #eta",
        'pt_2'  : "tau p_{T} (GeV)",
        'eta_2' : "tau #eta",
        'met': "E_{T}^{mis} (GeV)",
        'm_vis': "m_{vis} (GeV)",
	'ipsig_1': "electron IP sig",
        'ipsig_2': "tau IP sig",
        'n_jets': "number of jets",
        'n_bjets': "number of b-jets",
        'jpt_1': "leading jet p_{T} (GeV)",
        'jpt_2': "trailing jet p_{T} (GeV)",
        'mjj': "dijet mass (GeV)",
        'jdeta': "#Delta#eta(j,j)",
        'CMetQCD' : "C_{MET}^{QCD}",
        'CMetW' : "C_{MET}^{W}",
        'dR' : '#DeltaR(e,#tau)',
        'bdt_signal' : 'BDT_{sig}',
        'bdt_ditau' : 'BDT_{#tau#tau}',
        'bdt_fakes' : 'BDT_{fakes}',
        'aco_lep_pi': '#phi_{CP} (#pi)',
        'aco_lep_rho': '#phi_{CP} (#rho)',
        'aco_lep_a1_1pr': '#phi_{CP} (a_{1}(1-prong))',
        'aco_lep_a1_3pr': '#phi_{CP} (a_{1}(3-prong))',
    }
}

#################################
# Symmetrization of histogram   #
#################################
def symmetrize(hist):
    nbins = hist.GetNbinsX()
    bins2 = nbins//2
    for ib in range(1,bins2+1):
        ib2 = nbins + 1 - ib
        x = 0.5*(hist.GetBinContent(ib)+hist.GetBinContent(ib2))
        err = math.sqrt(hist.GetBinError(ib)**2+hist.GetBinError(ib2**2))/2
        hist.SetBinContent(ib,x)
        hist.SetBinContent(ib2,x)
        hist.SetBinError(ib,err)
        hist.SetBinError(ib2,err)



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
        if xDen>1e-5 and xNum>1e-5:
            ratio = xNum/xDen
            rNum = eNum/xNum
            rDen = eDen/xDen
            rratio = math.sqrt(rNum*rNum+rDen*rDen)
            eratio = rratio * ratio
        else:
            ratio = eNum
            eratio = eNum
            
        hist.SetBinContent(i,ratio)
        hist.SetBinError(i,eratio)
    return hist

def histoErrRatio(numHist,denHist,histName):
    nbins = numHist.GetNbinsX()
    hist = numHist.Clone(histName)
    for i in range(1,nbins+1):
        xNum = numHist.GetBinContent(i)
        eNum = numHist.GetBinError(i)
        xDen = denHist.GetBinContent(i)
        eDen = denHist.GetBinError(i)
        ratio = 1
        eratio = 0.1
        if xNum>1e-7 and xDen>1e-7:
            ratio = xNum/xDen
            rNum = eNum/xNum
            rDen = eDen/xDen
            eratio = ratio*math.sqrt(rNum*rNum+rDen*rDen)
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

def slice2DHisto(hist,xmin,xmax,bins,suffix):
    nbinsX = hist.GetNbinsX()
    nbinsY = hist.GetNbinsY()
    newbins = len(bins)-1
    name = hist.GetName()+'_'+suffix
    newhist = ROOT.TH1D(name,"",newbins,array('d',list(bins)))
    for xb in range(1,nbinsX+1):
        for yb in range(1,nbinsY+1):
            xcenter = hist.GetXaxis().GetBinCenter(xb)
            if xcenter>xmin and xcenter<xmax:
                ycenter = hist.GetYaxis().GetBinCenter(yb)
                bin_id = newhist.GetXaxis().FindBin(ycenter)
                xnew = newhist.GetBinContent(bin_id)
                enew = newhist.GetBinError(bin_id)
                xbin = hist.GetBinContent(xb,yb)
                ebin = hist.GetBinError(xb,yb)
                x_update = xnew + xbin
                e_update = math.sqrt(enew*enew+ebin*ebin)
                newhist.SetBinContent(bin_id,x_update)
                newhist.SetBinError(bin_id,e_update)
    return newhist
                

def rebin2DHisto(hist,binsX,binsY,suffix):
    nbinsX = hist.GetNbinsX()
    nbinsY = hist.GetNbinsY()
    newbinsX = len(binsX)-1
    newbinsY = len(binsY)-1
    name = hist.GetName()+"_"+suffix
    newhist = ROOT.TH2D(name,"",newbinsX,array('d',list(binsX)),newbinsY,array('d',list(binsY)))
    for xb in range(1,nbinsX+1):
        for yb in range(1,nbinsY+1):
            xcenter = hist.GetXaxis().GetBinCenter(xb)
            ycenter = hist.GetYaxis().GetBinCenter(yb)
            binx = newhist.GetXaxis().FindBin(xcenter)
            biny = newhist.GetYaxis().FindBin(ycenter)
            xbin = hist.GetBinContent(xb,yb)
            ebin = hist.GetBinError(xb,yb)
            xnew = newhist.GetBinContent(binx,biny)
            enew = newhist.GetBinError(binx,biny)
            x_update = xbin + xnew;
            e_update = math.sqrt(ebin*ebin + enew*enew);
            newhist.SetBinContent(binx,biny,x_update)
            newhist.SetBinError(binx,biny,e_update)
    return newhist

def ReplicaHist(hist,name):
    nbins = hist.GetNbinsX()
    xbins = []
    for ib in range(1,nbins+2):
        xbins.append(hist.GetXaxis().GetBinLowEdge(ib))
    histOutput = ROOT.TH1D(name,'',nbins,array('d',list(xbins)))
    return histOutput

def copyHist(inputHist,outputHist):
    nbins = inputHist.GetNbinsX()
    for ib in range(1,nbins+1):
        outputHist.SetBinContent(ib,inputHist.GetBinContent(ib))
        outputHist.SetBinError(ib,inputHist.GetBinError(ib))

# extracting Tag-and-Probe histos from ROOT file created by RunSelection.py macro
def extractTagProbeHistos(f,bins,era,channel,isSecond):
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
    region_labels = ['pass','fail','incl']
    hists = {}
    histPtBins = f.Get('ptBins')
    histEtaBins = f.Get('etaBins')
    nbinsPt = histPtBins.GetNbinsX()
    nbinsEta = histEtaBins.GetNbinsX()
    for sample in samples:
        for sign in sign_labels:
            for iso in iso_labels:
                for typ in type_labels:
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

    return hists


