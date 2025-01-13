import ROOT 
import math
from array import array
import numpy as np
import os

#############################
##### General settings ######
#############################

#########################
# folder with tuples    #
#########################
tupleFolder='/eos/cms/store/group/phys_tau/ksavva/For_Aliaksei/files/testingzpt'

############################
# base folder with outputs #
############################
baseFolder = '/afs/cern.ch/work/r/rasp/CPHiggs/IP'

eras = ['Run3_2022','Run3_2022EE','Run3_2023','Run3_2023BPix']
    
eraLumi = {
    "Run3_2022"     : 7980.4,
    "Run3_2022EE"   : 26671.7,
    "Run3_2023"     : 17794.0,
    "Run3_2023BPix" : 9451.0 
}

###############
# MC samples ##
###############

dy_samples = ['DYto2L_M_50_madgraphMLM','DYto2L_M_50_madgraphMLM_ext1','DYto2L_M_50_1J_madgraphMLM','DYto2L_M_50_2J_madgraphMLM','DYto2L_M_50_3J_madgraphMLM','DYto2L_M_50_4J_madgraphMLM']
top_samples = ['TTto2L2Nu','TTto2L2Nu_ext1','TTtoLNu2Q','TTtoLNu2Q_ext1','TTto4Q','TTto4Q_ext1']
vv_samples = ['WW','WZ','ZZ','ST_t_channel_top_4f_InclusiveDecays','ST_t_channel_antitop_4f_InclusiveDecays','ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1','ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1','ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1','ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1']
wjets_samples = ['WtoLNu_madgraphMLM','WtoLNu_madgraphMLM_ext1','WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM','WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM']

def samplesEra(era,dy_samples,top_samples,vv_samples,wjets_samples):
    
    dy = dy_samples
    top = top_samples
    vv  = vv_samples
    wjets = wjets_samples

    if era in ["Run3_2023", "Run3_2023BPix"]:
        dy.remove('DYto2L_M_50_madgraphMLM_ext1')
        top.remove('TTto2L2Nu_ext1')
        top.remove('TTtoLNu2Q_ext1')
        top.remove('TTto4Q_ext1')
        vv.remove('ST_tW_top_2L2Nu_ext1')
        vv.remove('ST_tW_antitop_2L2Nu_ext1')
        vv.remove('ST_tW_top_LNu2Q_ext1')
        if args.era == "Run3_2023":
            vv.remove('ST_tW_antitop_LNu2Q')
            vv.remove('ST_tW_antitop_LNu2Q_ext1')
            wjets.remove('WtoLNu_madgraphMLM_ext1')    
    
    return dy,top,vv,wjets
    
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
    'ipsig_2': [200,0.,10.]
}

ptbins = {
    'ee' : [25.,30.,40.,100000],
    'mm' : [21.,25.,30.,40.,100000],
    'et' : [25.,30.,40.,100000],
    'mt' : [21.,25.,30.,40.,100000]
}

etabins = {
    'ee' : [0.,1.48,2.5],
    'mm' : [0.,1.2, 2.4],
    'et' : [0.,1.48,2.5],
    'mt' : [0.,1.2, 2.4]
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


def createBins(nbins,xmin,xmax):
    binwidth = (xmax-xmin)/float(nbins)
    bins = []
    for i in range(0,nbins+1):
        xb = xmin + float(i)*binwidth
        bins.append(xb)
    return bins

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

def extractHistos(f,var,bins):
    samples = ['data','dy','top','vv','wjets']
    sign_labels = ['os','ss']
    iso_labels = ['iso','antiiso']
    typ_labels = ['lep','tau','had','all']
    hists = {}
    for sample in samples:
        for sign in sign_labels:
            for iso in iso_labels:
                for typ in typ_labels:
                    name = '%s_%s_%s_%s_%s'%(sample,var,sign,iso,typ)
                    histo = f.Get(name)
                    hists[name] = rebinHisto(histo,bins,'rebinned')
    return hists
