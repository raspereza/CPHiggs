import ROOT 
import math
from array import array
import numpy as np
import os
import CPHiggs.Analysis.utils as utils
from CPHiggs.Analysis.ScaleFactor import ScaleFactor
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1
import CPHiggs.Analysis.pv_utils as pv_utils
#from CPHiggs.fastMTT.fastmtt import fastmtt

# Cuts for Z->tau+tau and Z->ll selection
class AnalysisCuts:
    def __init__(self,**kwargs):
        self.mtCut = kwargs.get('mtCut',999999.)
        self.mvisUpperCut = kwargs.get('mvisUpperCut',999999.)
        self.mvisLowerCut = kwargs.get('mvisLowerCut',0.)
        
        self.etaLep1Cut = kwargs.get('etaLep1Cut',2.4)
        self.etaLep2Cut = kwargs.get('etaLep2Cut',2.3)
        
        self.ptLep1Cut = kwargs.get('ptLep1Cut',21.)
        self.ptLep2Cut = kwargs.get('ptLep2Cut',20.)

        self.ptSingleLepTrigger  = kwargs.get('ptSingleLepTrigger',26.) # 31 for single-e
        self.etaSingleLepTrigger = kwargs.get('etaSingleLepTrigger',2.4) # 2.4 for single-e

        self.ptLepCrossTrigger  = kwargs.get('ptLepCrossTrigger',21.)  # 25 for e-tau
        self.etaLepCrossTrigger = kwargs.get('etaLepCrossTrigger',2.1) # 2.1 for e-tau
        self.ptTauCrossTrigger  = kwargs.get('ptTauCrossTrigger',32.)  # 35 for e-tau
        self.etaTauCrossTrigger = kwargs.get('etaTauCrossTrigger',2.1) # 2.1 for e-tau

        # antiMu:
        # 1 - VLoose
        # 2 - Loose
        # 3 - Medium
        # 4 - Tight

        # antiE:
        # 1 - VVVLoose
        # 2 - VVLoose
        # 3 - VLoose
        # 4 - Loose
        # 5 - Medium
        # 6 - Tight
        # 7 - VTight
        # 8 - VVTight

        # antiJet:
        # 1 - VVVLoose
        # 2 - VVLoose
        # 3 - VLoose
        # 4 - Loose
        # 5 - Medium
        # 6 - Tight
        # 7 - VTight
        # 8 - VVTight

        # mt : 
        self.antiMu  = kwargs.get('antiMu',4)
        self.antiE   = kwargs.get('antiE',2)
        self.antiJet = kwargs.get('antiJet',7)
        self.useCrossTrigger = kwargs.get('useCrossTrigger',False)

        self.tauToRhoDE = kwargs.get('TauToRhoDE',0.2)
        
        self.isoLepCut = kwargs.get('isoLepCut',0.15)

        self.ipsigLepCut = kwargs.get('ipsigLepCut',1.0)
        self.ipsigTauCut = kwargs.get('ipsigTauCut',1.25)
        
        self.isoLepInverseLowerCut = kwargs.get('isoLepInverseLowerCut',0.20)
        self.isoLepInverseUpperCut = kwargs.get('isoLepInverseUpperCut',0.50)

        self.applyIPSigLep1Cut = kwargs.get('applyIPSigLep1Cut',False)
        self.applyIPSigLep2Cut = kwargs.get('applyIPSigLep2Cut',False)

        self.lepMomScale = kwargs.get('lepMomScale',0.002) # lepton mom scale unc. 0.002/0.025 for muons/electrons
        self.tauMomScale = kwargs.get('tauMomScale',0.02)  # tau mom scale unc. 
        self.lepTauMomScale = kwargs.get('lepTauMomScale',0.04) # lep->tau fake mom scale unc. 0.02/0.04 for muons/electrons

        print('')
        print("Setting cuts ->")

        print("mtCut",self.mtCut)
        print("mvisLowerCut",self.mvisLowerCut)
        print("mvisUpperCut",self.mvisUpperCut)
        
        print("etaLep1Cut",self.etaLep1Cut)
        print("etaLep2Cut",self.etaLep2Cut)

        print("ptLep1Cut",self.ptLep1Cut)
        print("ptLep2Cut",self.ptLep2Cut)

        print("isoLepCut",self.isoLepCut)
        print("isoLepInverseLowerCut",self.isoLepInverseLowerCut)
        print("isoLepInverseUpperCut",self.isoLepInverseUpperCut)
        
        print("ptSingleLepTrigger",self.ptSingleLepTrigger)
        print("etaSingleLepTrigger",self.etaSingleLepTrigger)

        print("ptLepCrossTrigger",self.ptLepCrossTrigger)
        print("etaLepCrossTrigger",self.etaLepCrossTrigger)
        print("ptTauCrossTrigger",self.ptTauCrossTrigger)
        print("etaTauCrossTrigger",self.etaTauCrossTrigger)
        
        print("antiMu",self.antiMu)
        print("antiE",self.antiE)
        print("antiJet",self.antiJet)
        print("useCrossTrigger",self.useCrossTrigger)

        print("ipsigLepCut",self.ipsigLepCut)
        print("ipsigTauCut",self.ipsigTauCut)
        print("applyIPSigLep1Cut",self.applyIPSigLep1Cut)
        print("applyIPSigLep2Cut",self.applyIPSigLep2Cut)
        print("TauToRhoDE",self.tauToRhoDE)
        
        print("lepMomScale",self.lepMomScale)
        print("tauMomScale",self.tauMomScale)
        print("lepTauMomScale",self.lepTauMomScale)

        print('')
        
# Run over set of samples and create histogram
def RunSamples(samples,var,cut,xbins,name,**kwargs):
    weight = kwargs.get('weight','weight')
    verbosity = kwargs.get('verbosity',False)
    if verbosity:
        print('')
        print("Running",name,var,cut)
    nbins = len(xbins)-1
    #    print(xbins)
    #    exit
    hist = ROOT.TH1D(name,"",nbins,array('d',list(xbins)))
    for sampleName in samples:
        sample = samples[sampleName]
        histsample = sample.CreateHisto(var,weight,cut,xbins,name)
        hist.Add(hist,histsample,1.,1.)
    return hist

# Run over set of samples and create histograms for Z->tautau channel
# for each sample loop over Tree entries is performed
def RunSamplesTuple(samples,name,**kwargs):

    hists = {} # discionary of histograms

    first = True
    for sampleName in samples:
        print("Running on sample %s"%(sampleName))
        sample = samples[sampleName]
        nameSample = sample.getSampleName()
        histsample = sample.CreateHistosTuple()
        if first:
            for hist in histsample:
                namehist = name + '_' + hist
                hists[namehist] = histsample[hist].Clone(namehist)
            first = False
        else:
            for hist in histsample:
                namehist = name	+ '_' +	hist
                hists[namehist].Add(hists[namehist],histsample[hist])

#    for hist in hists:
#        print(hist)
    return hists
        
class analysisSample:

    def __init__(self,basefolder,era,channel,samplename,norm,isdata,**kwargs):
        filename = basefolder + "/" + era + "/" + channel + "/" + samplename + "/nominal/merged.root"
        if not os.path.isfile(filename):
            print("")
            print('File %s is not found '%(filename))
            print('for specified era : %s'%(era))
            print('check if variable in util.py file is correctly set')
            print('or check naming of samples')
            print("")
            exit()
        self.sampleName = samplename+'_'+era
        self.sampleFile = ROOT.TFile(filename,"READ")
        self.channel = channel
        self.era = era
        self.norm = norm
        self.ismc = True
        self.isdata = isdata
        self.printout = kwargs.get('printout',False)
        self.analysisType = kwargs.get('analysisType','baseline')
        if self.analysisType not in ['baseline','ipSig','phiCP','jetFakes','datacards']:
            print('')
            print('Unknown analysis type : %s'%(self.analysisType))
            print('')
            exit()
        if isdata:
            self.norm = 1.0
            self.ismc = False

        self.sign_labels = ['os','ss']
        self.iso_labels = ['iso','antiiso','rest']
        self.type_labels = ['tau','lep','had','all']
        self.scale_unc = ['lepUp','lepDown','tauUp','tauDown','lfakeUp','lfakeDown']
        self.mt_regions = ['low_mt','high_mt']
        self.ff_tauids = ['loose','nominal']
        self.njet_labels = ['njets0','njets1','njets2']
        self.dm_labels = ['pi','rho','a1_1pr','a1_3pr']
        self.trig_labels = ['trig','nontrig']
        
        print('%s : %s : %s : norm = %7.3f'%(era,self.channel,self.sampleName,self.norm))
        
    def SetConfig(self,analysisCuts,histPtBins,histEtaBins,**kwargs):
        self.analysisCuts = analysisCuts

        self.histEtaBins = histEtaBins
        self.histPtBins = histPtBins

        self.nbinsPt = histPtBins.GetNbinsX()
        self.nbinsEta = histEtaBins.GetNbinsX()

        self.ptMin = histPtBins.GetBinLowEdge(1)
        self.ptMax = histPtBins.GetBinLowEdge(self.nbinsPt+1)

        self.etaMin = histEtaBins.GetBinLowEdge(1)
        self.etaMax = histEtaBins.GetBinLowEdge(self.nbinsEta+1)

        self.applyIPSigPromptLepSF = kwargs.get('applyIPSigPromptLepSF',False)
        self.applyIPSigTauLepSF = kwargs.get('applyIPSigTauLepSF',False)
        
        self.ipSigPromptLepSF = kwargs.get('ipSigPromptLepSF',None)
        self.ipSigTauLepSF = kwargs.get('ipSigTauLepSF',None)

        self.applyWeightCP = kwargs.get('applyWeightCP',0)
        
    def getSampleName(self):
        return self.sampleName


    def ComputeMass(self,pt_1,eta_1,phi_1,mass_1,
                    pt_2,eta_2,phi_2,mass_2):
        l1 = ROOT.TLorentzVector()
        l1.SetPtEtaPhiM(pt_1,eta_1,phi_1,mass_1)
        l2 = ROOT.TLorentzVector()
        l2.SetPtEtaPhiM(pt_2,eta_2,phi_2,mass_2)
        mass = (l1+l2).M()
        return mass
        
    
    def CreateHisto(self,var,weight,cut,bins,name):

        nbins = len(bins)-1
        histname = self.sampleName+'_'+name
        hist = ROOT.TH1D(histname,"",nbins,array('d',list(bins)))
        cutstring = weight+"*("+cut+")"
        tree = self.sampleFile.Get("tree")
        if (self.additionalCut!=''):
            cutstring = weight+"*("+cut+"&&"+self.additionalCut+")"
        tree.Draw(var+">>"+histname,cutstring)
        hist.Scale(self.norm)
        return hist

        
    def GetPtEtaBinLabels(self,pt,eta):
        ptX = pt
        etaX = abs(eta)

        if ptX<self.ptMin+0.1: ptX = self.ptMin+0.1
        if ptX>self.ptMax-0.1: ptX = self.ptMax-0.1

        if etaX<self.etaMin+0.01: etaX = self.etaMin+0.01
        if etaX>self.etaMax-0.01: etaX = self.etaMax-0.01

        binPt = self.histPtBins.GetXaxis().FindBin(ptX)
        binEta = self.histEtaBins.GetXaxis().FindBin(etaX)
        binLabel = '%1i_%1i'%(binPt,binEta)

#        print(pt,eta,binPt,binEta,binLabel)
        
        return binLabel,binPt,binEta
        
    def DeclareHistos(self):

        analysisType = self.analysisType
        
        hists = {}
        # general histos ->
        for var in utils.lib_histos:
            nbins = utils.lib_histos[var][0]
            xmin  = utils.lib_histos[var][1]
            xmax  = utils.lib_histos[var][2]
            width = (xmax-xmin)/float(nbins)
            xbins = []
            for i in range(0,nbins+1):
                xb = xmin + width*float(i)
                xbins.append(xb)
                
            for sign in self.sign_labels:
                for iso in self.iso_labels:
                    for typ in self.type_labels:
                        name = '%s_%s_%s_%s'%(var,sign,iso,typ)
                        histname = self.sampleName+'_'+name
                        hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))


        if analysisType=='phiCP':
            for var in utils.lib_phiCP_histos:
                nbins = utils.lib_histos[var][0]
                xmin  = utils.lib_histos[var][1]
                xmax  = utils.lib_histos[var][2]
                width = (xmax-xmin)/float(nbins)
                xbins = []
                for i in range(0,nbins+1):
                    xb = xmin + width*float(i)
                    xbins.append(xb)
                
                for sign in self.sign_labels:
                    for iso in self.iso_labels:
                        for typ in self.type_labels:
                            name = '%s_%s_%s_%s'%(var,sign,iso,typ)
                            histname = self.sampleName+'_'+name
                            hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))

            
        # measurements of SFs related to cut on mu/e ipSig
        if analysisType=='ipSig':
            nbins = utils.lib_histos['m_vis'][0]
            xmin  = utils.lib_histos['m_vis'][1]
            xmax  = utils.lib_histos['m_vis'][2]
            width = (xmax-xmin)/float(nbins)
            xbins = []
            for i in range(0,nbins+1):
                xb = xmin + width*float(i)
                xbins.append(xb)
            for sign in self.sign_labels:
                for iso in self.iso_labels:
                    for typ in self.type_labels:
                        for binPt in range(1,self.nbinsPt+1):
                            for binEta in range(1,self.nbinsEta+1):
                                label = '%1i_%1i'%(binPt,binEta)
                                for region in ['pass','fail','incl']:
                                    name = 'm_vis_%s_%s_%s_%s_%s'%(region,label,sign,iso,typ)
                                    histname = self.sampleName+'_'+ name
                                    hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                                    name = 'm_vis_2_%s_%s_%s_%s_%s'%(region,label,sign,iso,typ)
                                    histname = self.sampleName+'_'+ name
                                    hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                                    for unc in self.scale_unc:
                                        name = 'm_vis_%s_%s_%s_%s_%s_%s'%(unc,region,label,sign,iso,typ)
                                        histname = self.sampleName+'_'+ name
                                        hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))


        if analysisType=='datacards':
            nbins_bdt = 14
            nbins_phicp = 8
            xmin_bdt = 0.3
            xmax_bdt = 1.0
            xmin_phcp = 0.
            xmax_phicp = 2.0*ROOT.TMath.Pi()

        if analysisType=='jetFakes':
            xbins = []
            xbins.append(20.)
            nbins = 0
            for i in range(1,5):
                nbins += 1
                xb = 20.0+float(i)*5.0
                xbins.append(xb)
            for i in range(1,3):
                nbins += 1
                xb = 40.0+float(i)*10.0
                xbins.append(xb)
            for i in range(1,3):
                nbins += 1
                xb = 60.0+float(i)*20.0
                xbins.append(xb)
            nbins += 1
            xbins.append(150.)
            nbins += 1
            xbins.append(400.)
            for ff_tauid in self.ff_tauids:
                for mt in self.mt_regions:
                    for njet in self.njet_labels:
                        for dm in self.dm_labels:
                            for sign in self.sign_labels:
                                for iso in self.iso_labels:
                                    for typ in self.type_labels:
                                        name = 'pt_2_%s_%s_%s_%s_%s_%s_%s'%(ff_tauid,mt,njet,dm,sign,iso,typ)
                                        histname = self.sampleName+'_'+ name
                                        hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                                        
        return hists

    def CreateHistosTuple(self):

        tree = self.sampleFile.Get('ntuple')

        # initialization
        cuts = self.analysisCuts
        ipSigPromptLepSF = self.ipSigPromptLepSF
        ipSigTauLepSF = self.ipSigTauLepSF
        applyIPSigPromptLepSF = self.applyIPSigPromptLepSF
        applyIPSigTauLepSF = self.applyIPSigTauLepSF
        applyWeightCP = self.applyWeightCP
        analysisType = self.analysisType

#        print('SetConfig ->')
#        print('ptMin=%1.0f ptMax=%1.0f'%(self.ptMin,self.ptMax))
#        print('etaMin=%1.0f etaMax=%3.1f'%(self.etaMin,self.etaMax))
        
        era = self.era
        
        RadToDeg = 180./ROOT.TMath.Pi()
        one_over_sqrt2 = 1.0/math.sqrt(2.0)
        pi_over_4 = 0.25*ROOT.TMath.Pi()
        
        channel = self.channel
        sampleName = self.sampleName
        
        massLep = utils.electron_mass
        if channel=='mt':
            massLep = utils.muon_mass
        massPi = utils.pi_mass
        massPi0 = utils.pi0_mass
        
        # creating histograms 
        hists = self.DeclareHistos()

        # run event
        run    = np.zeros(1,dtype=np.int64)
        event  = np.zeros(1,dtype=np.int64)
        
        # floats
        w_Trigger   = np.zeros(1,dtype=np.float64)
        weight      = np.zeros(1,dtype=np.float64)

        # First lepton
        
        decayMode_1 = np.zeros(1,dtype=np.int64)
        decayModePNet_1 = np.zeros(1,dtype=np.int64)
        pt_1_FastMTT = np.zeros(1,dtype=np.float64)
        
        pt_1        = np.zeros(1,dtype=np.float64)
        eta_1       = np.zeros(1,dtype=np.float64)
        phi_1       = np.zeros(1,dtype=np.float64)
        iso_1       = np.zeros(1,dtype=np.float64)        
        mt_1        = np.zeros(1,dtype=np.float64)
        mass_1      = np.zeros(1,dtype=np.float64)
        charge_1    = np.zeros(1,dtype=np.float64)
        ipsig_1     = np.zeros(1,dtype=np.float64)

        ip_x_1      = np.zeros(1,dtype=np.float64)
        ip_y_1      = np.zeros(1,dtype=np.float64)
        ip_z_1      = np.zeros(1,dtype=np.float64)
        
        pi_pt_1     = np.zeros(1,dtype=np.float64)
        pi_eta_1    = np.zeros(1,dtype=np.float64)
        pi_phi_2    = np.zeros(1,dtype=np.float64)
        pi_mass_2   = np.zeros(1,dtype=np.float64)
        pi_charge_2 = np.zeros(1,dtype=np.float64)
        
        pi0_pt_1     = np.zeros(1,dtype=np.float64)
        pi0_eta_1    = np.zeros(1,dtype=np.float64)
        pi0_phi_1    = np.zeros(1,dtype=np.float64)
        pi0_mass_1   = np.zeros(1,dtype=np.float64)

        pi2_pt_1     = np.zeros(1,dtype=np.float64)
        pi2_eta_1    = np.zeros(1,dtype=np.float64)
        pi2_phi_1    = np.zeros(1,dtype=np.float64)
        pi2_mass_1   = np.zeros(1,dtype=np.float64)
        pi2_charge_1 = np.zeros(1,dtype=np.float64)
        
        pi3_pt_1     = np.zeros(1,dtype=np.float64)
        pi3_eta_1    = np.zeros(1,dtype=np.float64)
        pi3_phi_1    = np.zeros(1,dtype=np.float64)
        pi3_mass_1   = np.zeros(1,dtype=np.float64)
        pi3_charge_1 = np.zeros(1,dtype=np.float64)

        pion_E_split_1 = np.zeros(1,dtype=np.float64)
        
        sv_x_1       = np.zeros(1,dtype=np.float64)
        sv_y_1       = np.zeros(1,dtype=np.float64)
        sv_z_1       = np.zeros(1,dtype=np.float64)

        genPart_pt_1 = np.zeros(1,dtype=np.float64)
        genPart_eta_1 = np.zeros(1,dtype=np.float64)
        genPart_phi_1 =	np.zeros(1,dtype=np.float64)

        alphaAngle_1 = np.zeros(1,dtype=np.float64)
        
        # Second lepton ->
        
        decayMode_2 = np.zeros(1,dtype=np.int64)
        decayModePNet_2 = np.zeros(1,dtype=np.int64)
        pt_2_FastMTT = np.zeros(1,dtype=np.float64)

        pt_2        = np.zeros(1,dtype=np.float64)
        eta_2       = np.zeros(1,dtype=np.float64)
        phi_2       = np.zeros(1,dtype=np.float64)
        iso_2       = np.zeros(1,dtype=np.float64)
        mt_2        = np.zeros(1,dtype=np.float64)
        mass_2      = np.zeros(1,dtype=np.float64)
        charge_2    = np.zeros(1,dtype=np.int16)
        ipsig_2     = np.zeros(1,dtype=np.float64)
        
        ip_x_2      = np.zeros(1,dtype=np.float64)
        ip_y_2      = np.zeros(1,dtype=np.float64)
        ip_z_2      = np.zeros(1,dtype=np.float64)

        pi_pt_2     = np.zeros(1,dtype=np.float64)
        pi_eta_2    = np.zeros(1,dtype=np.float64)
        pi_phi_2    = np.zeros(1,dtype=np.float64)
        pi_mass_2   = np.zeros(1,dtype=np.float64)
        pi_charge_2 = np.zeros(1,dtype=np.float64)
        
        pi0_pt_2     = np.zeros(1,dtype=np.float64)
        pi0_eta_2    = np.zeros(1,dtype=np.float64)
        pi0_phi_2    = np.zeros(1,dtype=np.float64)
        pi0_mass_2   = np.zeros(1,dtype=np.float64)

        pi2_pt_2     = np.zeros(1,dtype=np.float64)
        pi2_eta_2    = np.zeros(1,dtype=np.float64)
        pi2_phi_2    = np.zeros(1,dtype=np.float64)
        pi2_mass_2   = np.zeros(1,dtype=np.float64)
        pi2_charge_2 = np.zeros(1,dtype=np.float64)
        
        pi3_pt_2     = np.zeros(1,dtype=np.float64)
        pi3_eta_2    = np.zeros(1,dtype=np.float64)
        pi3_phi_2    = np.zeros(1,dtype=np.float64)
        pi3_mass_2   = np.zeros(1,dtype=np.float64)
        pi3_charge_2 = np.zeros(1,dtype=np.float64)

        pion_E_split_2 = np.zeros(1,dtype=np.float64)
        
        sv_x_2       = np.zeros(1,dtype=np.float64)
        sv_y_2       = np.zeros(1,dtype=np.float64)
        sv_z_2       = np.zeros(1,dtype=np.float64)

        alphaAngle_2 = np.zeros(1,dtype=np.float64)

        # general variables
        
        m_vis        = np.zeros(1,dtype=np.float64)
        met_pt       = np.zeros(1,dtype=np.float64)
        met_phi      = np.zeros(1,dtype=np.float64)
        met_covXX    = np.zeros(1,dtype=np.float64)
        met_covXY    = np.zeros(1,dtype=np.float64)
        met_covYY    = np.zeros(1,dtype=np.float64)
        aco_lep_pi   = np.zeros(1,dtype=np.float64)
        aco_lep_rho  = np.zeros(1,dtype=np.float64)
        aco_lep_a1   = np.zeros(1,dtype=np.float64)
        aco_lep_a1_FastMTT = np.zeros(1,dtype=np.float64)
        weight_CP    = np.zeros(1,dtype=np.float64)

        n_jets       = np.zeros(1,dtype=np.float64)
        n_bjets      = np.zeros(1,dtype=np.float64)
        jdeta        = np.zeros(1,dtype=np.float64)
        mjj          = np.zeros(1,dtype=np.float64)
        jpt_1        = np.zeros(1,dtype=np.float64)
        jpt_2        = np.zeros(1,dtype=np.float64)

        PVBS_x       = np.zeros(1,dtype=np.float64)
        PVBS_y       = np.zeros(1,dtype=np.float64)
        PVBS_z       = np.zeros(1,dtype=np.float64)

        aco_pi_pi    = np.zeros(1,dtype=np.float64)
        aco_pi_rho   = np.zeros(1,dtype=np.float64)
        aco_pi_a1    = np.zeros(1,dtype=np.float64)

        aco_rho_pi   = np.zeros(1,dtype=np.float64)
        aco_rho_rho  = np.zeros(1,dtype=np.float64)
        aco_rho_a1   = np.zeros(1,dtype=np.float64)

        aco_a1_pi   = np.zeros(1,dtype=np.float64)
        aco_a1_rho  = np.zeros(1,dtype=np.float64)
        aco_a1_a1   = np.zeros(1,dtype=np.float64)

        aco_pi_a1_FastMTT = np.zeros(1,dtype=np.float64)
        aco_a1_pi_FastMTT = np.zeros(1,dtype=np.float64)

        aco_rho_a1_FastMTT = np.zeros(1,dtype=np.float64)
        aco_a1_rho_FastMTT = np.zeros(1,dtype=np.float64)

        aco_a1_a1_FastMTT = np.zeros(1,dtype=np.float64)
        
        # booleans
        trg_lep     = np.zeros(1,dtype='?')
        trg_lep2    = np.zeros(1,dtype='?')
        trg_cross   = np.zeros(1,dtype='?')
        trg_doubletau = np.zeros(1,dtype='?')
        os          = np.zeros(1,dtype='?')
        
        # integers
        idDeepTau2018v2p5VSe_1   = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSmu_1  = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSjet_1 = np.zeros(1,dtype=np.int64)
        hasRefitSV_1 = np.zeros(1,dtype='?')

        idDeepTau2018v2p5VSe_2   = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSmu_2  = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSjet_2 = np.zeros(1,dtype=np.int64)
        hasRefitSV_2 = np.zeros(1,dtype='?')

        genmatch_1               = np.zeros(1,dtype=np.int64)
        genmatch_2               = np.zeros(1,dtype=np.int64)

        # branches ->

        tree.SetBranchAddress('run',run)
        tree.SetBranchAddress('event',event)
        
        # floats ->
        tree.SetBranchAddress('met_pt',met_pt)
        tree.SetBranchAddress('met_phi',met_phi)
        tree.SetBranchAddress('met_covXX',met_covXX)
        tree.SetBranchAddress('met_covXY',met_covXY)
        tree.SetBranchAddress('met_covYY',met_covYY)
        tree.SetBranchAddress('m_vis',m_vis)
        tree.SetBranchAddress('os',os)

        tree.SetBranchAddress('n_jets',n_jets)
        tree.SetBranchAddress('n_bjets',n_bjets)
        tree.SetBranchAddress('jpt_1',jpt_1)
        tree.SetBranchAddress('jpt_2',jpt_2)
        tree.SetBranchAddress('jdeta',jdeta)
        tree.SetBranchAddress('mjj',mjj)

        tree.SetBranchAddress('mt_1',mt_1)
        tree.SetBranchAddress('pt_1',pt_1)
        tree.SetBranchAddress('eta_1',eta_1)
        tree.SetBranchAddress('phi_1',phi_1)
        tree.SetBranchAddress('mass_1',mass_1)
        tree.SetBranchAddress('charge_1',charge_1)
        tree.SetBranchAddress('iso_1',iso_1)
        
        tree.SetBranchAddress('ip_x_1',ip_x_1)
        tree.SetBranchAddress('ip_y_1',ip_y_1)
        tree.SetBranchAddress('ip_z_1',ip_z_1)
        tree.SetBranchAddress('ip_LengthSig_1',ipsig_1)

        tree.SetBranchAddress('pt_2',pt_2)
        tree.SetBranchAddress('eta_2',eta_2)
        tree.SetBranchAddress('phi_2',phi_2)
        tree.SetBranchAddress('mass_2',mass_2)
        tree.SetBranchAddress('charge_2',charge_2)

        tree.SetBranchAddress('ip_x_2',ip_x_2)
        tree.SetBranchAddress('ip_y_2',ip_y_2)
        tree.SetBranchAddress('ip_z_2',ip_z_2)
        tree.SetBranchAddress('ip_LengthSig_2',ipsig_2)

        if channel=='mt' or channel=='et' or channel=='tt':
            
            tree.SetBranchAddress('PVBS_x',PVBS_x)
            tree.SetBranchAddress('PVBS_y',PVBS_y)
            tree.SetBranchAddress('PVBS_z',PVBS_z)

            if analysisType=='phiCP':
                tree.SetBranchAddress('pi0_pt_2',pi0_pt_2)
                tree.SetBranchAddress('pi0_eta_2',pi0_eta_2)
                tree.SetBranchAddress('pi0_phi_2',pi0_phi_2)

                tree.SetBranchAddress('pi_pt_2',pi_pt_2)
                tree.SetBranchAddress('pi_eta_2',pi_eta_2)
                tree.SetBranchAddress('pi_phi_2',pi_phi_2)
                tree.SetBranchAddress('pi_mass_2',pi_mass_2)
                tree.SetBranchAddress('pi_charge_2',pi_charge_2)

                tree.SetBranchAddress('pi2_pt_2',pi2_pt_2)
                tree.SetBranchAddress('pi2_eta_2',pi2_eta_2)
                tree.SetBranchAddress('pi2_phi_2',pi2_phi_2)
                tree.SetBranchAddress('pi2_mass_2',pi2_mass_2)
                tree.SetBranchAddress('pi2_charge_2',pi2_charge_2)
            
                tree.SetBranchAddress('pi3_pt_2',pi3_pt_2)
                tree.SetBranchAddress('pi3_eta_2',pi3_eta_2)
                tree.SetBranchAddress('pi3_phi_2',pi3_phi_2)
                tree.SetBranchAddress('pi3_mass_2',pi3_mass_2)
                tree.SetBranchAddress('pi3_charge_2',pi3_charge_2)

            tree.SetBranchAddress('pion_E_split_2',pion_E_split_2)
            
            tree.SetBranchAddress('sv_x_2',sv_x_2)
            tree.SetBranchAddress('sv_y_2',sv_y_2)
            tree.SetBranchAddress('sv_z_2',sv_z_2)

            tree.SetBranchAddress('hasRefitSV_2',hasRefitSV_2)

        if channel=='tt':

            if analysisType=='phiCP':
                tree.SetBranchAddress('pi0_pt_1',pi0_pt_1)
                tree.SetBranchAddress('pi0_eta_1',pi0_eta_1)
                tree.SetBranchAddress('pi0_phi_1',pi0_phi_1)

                tree.SetBranchAddress('pi_pt_1',pi_pt_1)
                tree.SetBranchAddress('pi_eta_1',pi_eta_1)
                tree.SetBranchAddress('pi_phi_1',pi_phi_1)
                tree.SetBranchAddress('pi_mass_1',pi_mass_1)
                tree.SetBranchAddress('pi_charge_1',pi_charge_1)

                tree.SetBranchAddress('pi2_pt_1',pi2_pt_1)
                tree.SetBranchAddress('pi2_eta_1',pi2_eta_1)
                tree.SetBranchAddress('pi2_phi_1',pi2_phi_1)
                tree.SetBranchAddress('pi2_mass_1',pi2_mass_1)
                tree.SetBranchAddress('pi2_charge_1',pi2_charge_1)
            
                tree.SetBranchAddress('pi3_pt_1',pi3_pt_1)
                tree.SetBranchAddress('pi3_eta_1',pi3_eta_1)
                tree.SetBranchAddress('pi3_phi_1',pi3_phi_1)
                tree.SetBranchAddress('pi3_mass_1',pi3_mass_1)
                tree.SetBranchAddress('pi3_charge_1',pi3_charge_1)

            tree.SetBranchAddress('pion_E_split_1',pion_E_split_1)
            
            tree.SetBranchAddress('sv_x_1',sv_x_1)
            tree.SetBranchAddress('sv_y_1',sv_y_1)
            tree.SetBranchAddress('sv_z_1',sv_z_1)

            tree.SetBranchAddress('hasRefitSV_1',hasRefitSV_1)

        if channel=='mt':
            tree.SetBranchAddress('alphaAngle_mu_pi_2',alphaAngle_2)
        if channel=='et':
            tree.SetBranchAddress('alphaAngle_e_pi_2',alphaAngle_2)

        if channel=='mm' or channel=='ee':
            tree.SetBranchAddress('iso_2',iso_2)

        if channel=='mt':
            tree.SetBranchAddress('aco_mu_pi',aco_lep_pi)
            tree.SetBranchAddress('aco_mu_rho',aco_lep_rho)
            tree.SetBranchAddress('aco_mu_a1',aco_lep_a1)


        if channel=='et':
            tree.SetBranchAddress('aco_e_pi',aco_lep_pi)
            tree.SetBranchAddress('aco_e_rho',aco_lep_rho)
            tree.SetBranchAddress('aco_e_a1',aco_lep_a1)

        if channel=='tt':
            
            tree.SetBranchAddress('aco_pi_pi',aco_pi_pi)
            tree.SetBranchAddress('aco_pi_rho',aco_pi_rho)
            tree.SetBranchAddress('aco_pi_a1',aco_pi_a1)
            
            tree.SetBranchAddress('aco_rho_pi',aco_rho_pi)
            tree.SetBranchAddress('aco_rho_rho',aco_rho_rho)
            tree.SetBranchAddress('aco_rho_a1',aco_rho_a1)
            
            tree.SetBranchAddress('aco_a1_pi',aco_a1_pi)
            tree.SetBranchAddress('aco_a1_rho',aco_a1_rho)
            tree.SetBranchAddress('aco_a1_a1',aco_a1_a1)

        if channel=='mt' or channel=='et' or channel=='tt':
            tree.SetBranchAddress('decayMode_2',decayMode_2)
            tree.SetBranchAddress('decayModePNet_2',decayModePNet_2)

        if channel=='tt':
            tree.SetBranchAddress('decayMode_1',decayMode_1)
            tree.SetBranchAddress('decayModePNet_1',decayModePNet_1)

        if self.ismc:
            if channel=='mm':
                tree.SetBranchAddress('w_Trigger',w_Trigger)

        tree.SetBranchAddress('weight',weight)

        if applyWeightCP==1: # SM (CP-even)
            tree.SetBranchAddress('wt_cp_sm',weight_CP)
        if applyWeightCP==2: # PS (CP-odd)
            tree.SetBranchAddress('wt_cp_ps',weight_CP)
        
        # booleans (trigger)
        if channel=='mt' or channel=='mm': tree.SetBranchAddress('trg_singlemuon',trg_lep)
        if channel=='mm': tree.SetBranchAddress('trg_singlemuon_2',trg_lep2)
        if channel=='mt': tree.SetBranchAddress('trg_mt_cross',trg_cross)
        if channel=='et' or channel=='ee': tree.SetBranchAddress('trg_singleelectron',trg_lep)
        if channel=='et': tree.SetBranchAddress('trg_et_cross',trg_cross)
        if channel=='tt': tree.SetBranchAddress('trg_doubletau',trg_doubletau)
        
        # integers
        if channel=='mt' or channel=='et' or channel=='tt':
            tree.SetBranchAddress('idDeepTau2018v2p5VSe_2',idDeepTau2018v2p5VSe_2)
            tree.SetBranchAddress('idDeepTau2018v2p5VSmu_2',idDeepTau2018v2p5VSmu_2)
            tree.SetBranchAddress('idDeepTau2018v2p5VSjet_2',idDeepTau2018v2p5VSjet_2)

        if channel=='tt':
            tree.SetBranchAddress('idDeepTau2018v2p5VSe_1',idDeepTau2018v2p5VSe_1)
            tree.SetBranchAddress('idDeepTau2018v2p5VSmu_1',idDeepTau2018v2p5VSmu_1)
            tree.SetBranchAddress('idDeepTau2018v2p5VSjet_1',idDeepTau2018v2p5VSjet_1)
            
        if self.ismc:
            tree.SetBranchAddress("genPartFlav_1",genmatch_1)
            tree.SetBranchAddress("genPartFlav_2",genmatch_2)
            
        nentries = tree.GetEntries()

        # run over entries
        for entry in range(0,nentries):

            if entry%100000==0:
                print('processed %1i out of %1i events'%(entry,nentries))
            tree.GetEntry(entry)

            
            # trigger threshold
            passTrigger = False
            if channel=='mm' or channel=='ee':
                pt1Trig = pt_1[0]
                trg1 = trg_lep[0]
                trg2 = trg_lep2[0]
                eta1Trig = abs(eta_1[0])
                pt2Trig = pt_2[0]
                eta2Trig = abs(eta_2[0])
                trigger1 = pt1Trig>cuts.ptSingleLepTrigger and eta1Trig<cuts.etaSingleLepTrigger and trg1
                trigger2 = pt2Trig>cuts.ptSingleLepTrigger and eta2Trig<cuts.etaSingleLepTrigger and trg2
                passTrigger = trigger1
            if channel=='mt' or channel=='et':
                trig_lep_acc = pt_1[0]>cuts.ptSingleLepTrigger and abs(eta_1[0])<cuts.etaSingleLepTrigger
                passSingleLepTrigger = trig_lep_acc and trg_lep[0]
                passTrigger = passSingleLepTrigger
                if cuts.useCrossTrigger:
                    trig_l_ltau = pt_1[0]>cuts.ptLepCrossTrigger and pt_1[0]<cuts.ptSingleLepTrigger and abs(eta_1[0])<cuts.etaLepCrossTrigger
                    trig_tau_ltau = pt_2[0]>cuts.ptTauCrossTrigger and abs(eta_2[0])<cuts.etaTauCrossTrigger
                    passCrossTrigger = trig_l_ltau and trig_tau_ltau and trg_cross[0]
                    passTrigger = passTrigger or passCrossTrigger

            if channel=='tt':
                trig_ditau_acc = pt_1[0]>cuts.ptDiTauTrigger and pt_2[0]>cuts.ptDiTauTrigger
                passTrigger = trig_ditau_acc and trg_doubletau[0]

            if not passTrigger: continue
            
            # kinematic cuts
            if pt_1[0]<cuts.ptLep1Cut: continue
            if abs(eta_1[0])>cuts.etaLep1Cut: continue
            if channel=='mt' or channel=='et':
                if iso_1[0]>cuts.isoLepInverseUpperCut: continue
            
            if pt_2[0]<cuts.ptLep2Cut: continue
            if abs(eta_2[0])>cuts.etaLep2Cut: continue

            if cuts.applyIPSigLep1Cut:
                if abs(ipsig_1[0])<cuts.ipsigLepCut: continue
            
            if channel=='mm' or channel=='ee':
                if iso_2[0]>cuts.isoLepCut: continue
                if cuts.applyIPSigLep2Cut:
                    if abs(ipsig_2[0])<cuts.ipsigLepCut: continue
            
            if channel=='mt' or channel=='et':
                # m_vis cut
                if m_vis[0]>cuts.mvisUpperCut: continue
                if m_vis[0]<cuts.mvisLowerCut: continue
                # tau discriminator against e and mu and jet
                if idDeepTau2018v2p5VSe_2[0]<cuts.antiE: continue
                if idDeepTau2018v2p5VSmu_2[0]<cuts.antiMu: continue
                if idDeepTau2018v2p5VSjet_2[0]<cuts.antiJet: continue

            ##### Decay mode specific selection #####
            isTauToPi_2 = decayModePNet_2[0]==0 and abs(ipsig_2[0])>cuts.ipsigTauCut        
            isTauToRho_2 = decayModePNet_2[0]==1 and decayMode_2[0]==1 and abs(pion_E_split_2[0])>cuts.tauToRhoDE
            isTauToA1_1pr_2 = decayModePNet_2[0]==2 and decayMode_2[0]==1 and abs(pion_E_split_2[0])>cuts.tauToRhoDE
            isTauToA1_3pr_2 = decayModePNet_2[0]==10 and hasRefitSV_2[0]

            passTau_2 =  isTauToPi_2 or isTauToRho_2 or isTauToA1_1pr_2 or isTauToA1_3pr_2
            if not passTau_2: continue
            
            if channel=='tt':
                isTauToPi_1 = decayModePNet_1[0]==0 and abs(ipsig_1[0])>cuts.ipsigTauCut        
                isTauToRho_1 = decayModePNet_1[0]==1 and decayMode_1[0]==1 and abs(pion_E_split_1[0])>cuts.TauToRhoDE
                isTauToA1_1pr_1 = decayModePNet_1[0]==2 and decayMode_1[0]==1 and abs(pion_E_split_1[0])>cuts.TauToRhoDE
                isTauToA1_3pr_1 = decayModePNet_1[0]==10 and hasRefitSV_1[0]

                passTau_1 =  isTauToPi_1 or isTauToRho_1 or isTauToA1_1pr_1 or isTauToA1_3pr_1 
                if not passTau_1: continue
                
            
            variables = {} 
            variables['m_vis'] = m_vis[0]
            variables['pt_1'] = pt_1[0]
            variables['eta_1'] = eta_1[0]
            variables['pt_2'] = pt_2[0]
            variables['eta_2'] = eta_2[0]
            variables['met'] = met_pt[0]
            variables['ipsig_1'] = abs(ipsig_1[0])
            variables['ipsig_2'] = abs(ipsig_2[0])

            # jet related variables
            variables['n_jets'] = n_jets[0]
            variables['n_bjets'] = n_bjets[0]
            variables['jpt_1'] = -9999.
            variables['jpt_2'] = -9999.
            variables['mjj'] = -9999.
            variables['jdeta'] = -9999.
            if n_jets[0]>0.5:
                variables['jpt_1'] = jpt_1[0]
            if n_jets[0]>1.5:
                variables['jpt_2'] = jpt_2[0]
                variables['mjj'] = mjj[0]
                variables['jdeta'] = abs(jdeta[0])
            

            sign_label = 'ss'
            lep_label = 'had'
            iso_label = 'rest'
            lep2_label = 'had'
            
            if self.ismc:
                if genmatch_1[0]==1: lep_label = 'lep'
                if genmatch_1[0]==15: lep_label = 'tau'
                if channel=='mm' or channel=='ee':
                    if genmatch_2[0]==1: lep2_label = 'lep'
                    if genmatch_2[0]==15: lep2_label = 'tau'
                else:
                    if genmatch_2[0]>=1 and genmatch_2[0]<=4: lep2_label = 'lep'
                    elif genmatch_2[0]==5: lep2_label = 'tau'
                    else: lep2_label = 'had' 
                    
                
            if os[0]:
                sign_label = 'os'
                
            directIso = iso_1[0] < cuts.isoLepCut
            inverseIso = iso_1[0] > cuts.isoLepInverseLowerCut and iso_1[0] < cuts.isoLepInverseUpperCut

            if directIso: iso_label = 'iso'
            if inverseIso: iso_label = 'antiiso'

            ##################
            ## total weight ##
            ##################
            Weight = weight[0]
            if self.ismc:
                if channel=='mm':
                    trg_weight = max(0.1,w_Trigger[0])
                    Weight /= trg_weight

            ##### CP weight (Higgs samples)
            if applyWeightCP==1 or applyWeightCP==2:
                Weight *= weight_CP[0]
            
            # prevent large weights
            if ROOT.TMath.Abs(Weight)>10.0:
                print('weight %3.1f > 10'%(Weight))
                continue

            
            ############################
            ## applying scale factors ##
            ## for IPSig cuts         ##
            ############################
            WeightSF = 1.0
            if cuts.applyIPSigLep1Cut:
                if lep_label=='lep' and applyIPSigPromptLepSF and self.ismc: 
                    WeightSF *= ipSigPromptLepSF.getSF(pt_1[0],eta_1[0])
                if lep_label=='tau' and applyIPSigTauLepSF and self.ismc:
                    WeightSF *= ipSigTauLepSF.getSF(pt_1[0],eta_1[0])

            if channel=='mm' or channel=='ee':
                if cuts.applyIPSigLep2Cut:
                    if applyIPSigPromptLepSF and self.ismc:
                        WeightSF *= ipSigPromptLepSF.getSF(pt_2[0],eta_2[0])


            if channel=='mt' or channel=='et':
                ########################################
                # Filling histogram with mt_1 variable #
                # before applying mt_1 cut             #
                ########################################
                nameAll = 'mt_1_%s_%s_all'%(sign_label,iso_label)
                name = 'mt_1_%s_%s_%s'%(sign_label,iso_label,lep_label)
                hists[nameAll].Fill(mt_1[0],Weight*WeightSF)
                hists[name].Fill(mt_1[0],Weight*WeightSF)
                # mt_1 cut
                if (mt_1[0]>cuts.mtCut): continue

            ################################
            ## Filling control histograms ##
            ################################
            for varname in variables:
                nameAll = '%s_%s_%s_all'%(varname,sign_label,iso_label)
                name = '%s_%s_%s_%s'%(varname,sign_label,iso_label,lep_label)
                hists[nameAll].Fill(variables[varname],Weight*WeightSF)
                hists[name].Fill(variables[varname],Weight*WeightSF)

#            print('genmatch_1 = %1i    genmatch_2 = %1i'%(genmatch_1[0],genmatch_2[0]))
            
            if analysisType=='ipSig':
                #################################################
                # tag-and-probe histograms for mm (ee) channels #
                # ###############################################
                if channel=='mm' or channel=='ee':
                                
                    # lep1 ->
                    bin_label,binPt,binEta = self.GetPtEtaBinLabels(pt_1[0],eta_1[0])
                    region_label = 'fail'
                    if variables['ipsig_1']>cuts.ipsigLepCut:
                        region_label = 'pass'
                    name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                    name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                
                    # lep2 ->
                    bin_label,binPt,binEta = self.GetPtEtaBinLabels(pt_2[0],eta_2[0])
                    region_label = 'fail'
                    if variables['ipsig_2']>cuts.ipsigLepCut:
                        region_label = 'pass'
                    name = 'm_vis_2_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep2_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                    name = 'm_vis_2_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                
                #################################################
                # tag-and-probe histograms for mt (et) channels #
                # ###############################################
                if channel=='mt' or channel=='et':
                    # scale uncertainties ->
                    mvis_unc = {}
                    
                    # lepton scale ->
                    pt_up = pt_1[0]*(1.0 + cuts.lepMomScale)
                    pt_down = pt_1[0]*(1.0 - cuts.lepMomScale)
                    mvis_unc['lepUp'] = self.ComputeMass(pt_up,eta_1[0],phi_1[0],mass_1[0],
                                                         pt_2[0],eta_2[0],phi_2[0],mass_2[0])
                    mvis_unc['lepDown'] = self.ComputeMass(pt_down,eta_1[0],phi_1[0],mass_1[0],
                                                           pt_2[0],eta_2[0],phi_2[0],mass_2[0])
                    
                    # tau scale ->
                    if lep2_label=='tau':
                        pt_up = pt_2[0]*(1.0 + cuts.tauMomScale)
                        pt_down = pt_2[0]*(1.0 - cuts.tauMomScale)
                        mvis_unc['tauUp'] = self.ComputeMass(pt_1[0],eta_1[0],phi_1[0],mass_1[0],
                                                             pt_up,eta_2[0],phi_2[0],mass_2[0])
                        mvis_unc['tauDown'] = self.ComputeMass(pt_1[0],eta_1[0],phi_1[0],mass_1[0],
                                                               pt_down,eta_2[0],phi_2[0],mass_2[0])
                        mvis_unc['lfakeUp'] = variables['m_vis']
                        mvis_unc['lfakeDown'] = variables['m_vis']
                    elif lep2_label=='lep':
                        pt_up = pt_2[0]*(1.0 + cuts.lepTauMomScale)
                        pt_down = pt_2[0]*(1.0 - cuts.lepTauMomScale)
                        mvis_unc['tauUp'] = variables['m_vis']
                        mvis_unc['tauDown'] = variables['m_vis']
                        mvis_unc['lfakeUp'] = self.ComputeMass(pt_1[0],eta_1[0],phi_1[0],mass_1[0],
                                                               pt_up,eta_2[0],phi_2[0],mass_2[0])
                        mvis_unc['lfakeDown'] = self.ComputeMass(pt_1[0],eta_1[0],phi_1[0],mass_1[0],
                                                                 pt_down,eta_2[0],phi_2[0],mass_2[0])
                    else:
                        mvis_unc['tauUp'] = variables['m_vis']
                        mvis_unc['tauDown'] = variables['m_vis']
                        mvis_unc['lfakeUp'] = variables['m_vis']
                        mvis_unc['lfakeDown'] = variables['m_vis']
                        
                    # lep1 ->
                    bin_label,binPt,binEta = self.GetPtEtaBinLabels(pt_1[0],eta_1[0])
                    # inclusive selection
                    region_label = 'incl'
                    name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                    name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                    for unc in mvis_unc:
                        name = 'm_vis_%s_%s_%s_%s_%s_%s'%(unc,region_label,bin_label,sign_label,iso_label,lep_label)
                        hists[name].Fill(mvis_unc[unc],Weight)
                        name = 'm_vis_%s_%s_%s_%s_%s_all'%(unc,region_label,bin_label,sign_label,iso_label)
                        hists[name].Fill(mvis_unc[unc],Weight)
                
                    if variables['ipsig_1']>cuts.ipsigLepCut:
                        # passing probes
                        WeightSF = 1.0
                        if lep_label=='lep' and self.ismc and ipSigPromptLepSF!=None:
                            WeightSF *= ipSigPromptLepSF.getSF(pt_1[0],eta_1[0])
                            region_label = 'pass'
                            name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                            hists[name].Fill(variables['m_vis'],Weight*WeightSF)
                            name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                            hists[name].Fill(variables['m_vis'],Weight*WeightSF)
                            for unc in mvis_unc:
                                name = 'm_vis_%s_%s_%s_%s_%s_%s'%(unc,region_label,bin_label,sign_label,iso_label,lep_label)
                                hists[name].Fill(mvis_unc[unc],Weight*WeightSF)
                                name = 'm_vis_%s_%s_%s_%s_%s_all'%(unc,region_label,bin_label,sign_label,iso_label)
                                hists[name].Fill(mvis_unc[unc],Weight*WeightSF)
                        

                    else:
                        # failing probes
                        region_label = 'fail'
                        name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                        hists[name].Fill(variables['m_vis'],Weight)
                        name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                        hists[name].Fill(variables['m_vis'],Weight)
                        for unc in mvis_unc:
                            name = 'm_vis_%s_%s_%s_%s_%s_%s'%(unc,region_label,bin_label,sign_label,iso_label,lep_label)
                            hists[name].Fill(mvis_unc[unc],Weight)
                            name = 'm_vis_%s_%s_%s_%s_%s_all'%(unc,region_label,bin_label,sign_label,iso_label)
                            hists[name].Fill(mvis_unc[unc],Weight)

            ####################################
            ## phi(CP) studies with signal MC ##
            ####################################
            if analysisType=='phiCP':
                if channel=='mt' or channel=='et':
                    variablesCP['aco_lep_pi_plus'] = -9999.
                    variablesCP['aco_lep_pi_minus'] = -9999.
                    variablesCP['aco_lep_pi'] = -9999.
                    variablesCP['aco_lep_piIP'] = -9999.
                    variablesCP['aco_lepIP_pi'] = -9999.
                    variablesCP['aco_lepIP_piIP'] = -9999.
            
                    variablesCP['aco_lep_rho_plus'] = -9999.
                    variablesCP['aco_lep_rho_minus'] = -9999.
                    variablesCP['aco_lep_rho'] = -9999.
                    variablesCP['aco_lepIP_rho'] = -9999.
                    variablesCP['aco_lep_rhoECut'] = -9999.
                    variablesCP['aco_lep_rhoGen'] = -9999.
                    variablesCP['aco_lep_rhoReco'] = -9999.
                    variablesCP['aco_lep_rhoCollinear'] = -9999.
                    variablesCP['aco_lep_rhoRecoGen'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenECut'] = -9999.
                    
                    variablesCP['aco_lep_rhoRecoIP1p0'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenIP1p0'] = -9999.
                    variablesCP['aco_lep_rhoRecoIP1p0ECut'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenIP1p0ECut'] = -9999.
                    
                    variablesCP['aco_lep_rhoRecoIP1p2'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenIP1p2'] = -9999.
                    variablesCP['aco_lep_rhoRecoIP1p2ECut'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenIP1p2ECut'] = -9999.
                    
                    variablesCP['aco_lep_rhoRecoIP1p5'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenIP1p5'] = -9999.
                    variablesCP['aco_lep_rhoRecoIP1p5ECut'] = -9999.
                    variablesCP['aco_lep_rhoRecoGenIP1p5ECut'] = -9999.
                
                    variablesCP['aco_lep_a1_plus'] = -9999.
                    variablesCP['aco_lep_a1_minus'] = -9999.
                    variablesCP['aco_lep_a1'] = -9999.
                    variablesCP['aco_lepIP_a1'] = -9999.
                    variablesCP['aco_lep_a1DP'] = -9999.
                    variablesCP['aco_lep_a1PVGen'] = -9999.
                    variablesCP['aco_lepIP_a1PVGen'] = -9999.
                    variablesCP['aco_lep_a1PVDESY'] = -9999.
                    variablesCP['aco_lepIP_a1PVDESY'] = -9999.
                    variablesCP['aco_lep_a1PVIC'] = -9999.
                    
                    variablesCP['alpha_lep_pi'] = -9999.
                    variablesCP['alpha_lep_rho'] = -9999.
                    variablesCP['alpha_lep_a1'] = -9999.
                    
                    alpha = -9999.
                    cosa = pv_utils.CosAlpha(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi,
                                         ip_x_2[0],ip_y_2[0],ip_z_2[0],False)
                    alpha_IP = ROOT.TMath.ACos(cosa)
                    aco = -9999.
                    pv = ROOT.TLorentzVector()
                    aco_a1 = -9999.
                
                    if isTauToPi_2:
                        # tau->pi+v decay 
                        alpha = alpha_IP
                        variablesCP['alpha_lep_pi'] = RadToDeg * alpha
                        variablesCP['aco_lep_pi'] = RadToDeg*aco_lep_pi[0]
                        if abs(ipsig_2[0])>cuts.ipsigTauCut:
                            variablesCP['aco_lep_piIP'] = RadToDeg*aco_lep_pi[0]
                        if abs(ipsig_1[0])>cuts.ipsigLepCut:
                            variablesCP['aco_lepIP_pi'] = RadToDeg*aco_lep_pi[0]
                        if abs(ipsig_1[0])>cuts.ipsigLepCut and abs(ipsig_2[0])>cuts.ipsigLepCut:
                            variablesCP['aco_lepIP_piIP'] = RadToDeg*aco_lep_pi[0]
                        if alpha > pi_over_4:
                            variablesCP['aco_lep_pi_plus'] = RadToDeg*aco_lep_pi[0]
                        else:
                            variablesCP['aco_lep_pi_minus'] = RadToDeg*aco_lep_pi[0]

                    if isTauToRho_2:
                        # tau->rho(p+pi-)+v decay
                        Pi0 = ROOT.TLorentzVector()
                        Pi0.SetPtEtaPhiM(pi0_pt_2[0],pi0_eta_2[0],pi0_phi_2[0],massPi0)
                        Pi = ROOT.TLorentzVector()
                        Pi.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                        magPi0 = Pi0.Vect().Mag()
                        magPi = Pi.Vect().Mag()
                        cosa = pv_utils.CosAlpha(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi,Pi0.Px(),Pi0.Py(),Pi0.Pz(),True)
                        alpha = ROOT.TMath.ACos(cosa)
                        variablesCP['alpha_lep_rho'] = RadToDeg * alpha
                        P1 = ROOT.TLorentzVector()
                        P1.SetPtEtaPhiM(pt_1_FastMTT[0],eta_1[0],phi_1[0],massLep)
                        R1 = ROOT.TLorentzVector()
                        R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0.0)
                        ip_vec_2 = ROOT.TVector3(ip_x_2[0],ip_y_2[0],ip_z_2[0])
                        PGen = ROOT.TLorentzVector()
                        PGen.SetPtEtaPhiM(genPart_pt_2[0],genPart_eta_2[0],genPart_phi_2[0],utils.tau_mass)

                        deltaE = abs((Pi.E()-Pi0.E())/(Pi.E()+Pi0.E()))
                    
                        # Gen
                        P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'gen')
                        firstNeg = charge_1[0] < 0.
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                        if math.isnan(aco): aco = -9999.
                        variablesCP['aco_lep_rhoGen'] = RadToDeg*aco

                        # Collinear
                        P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'collinear')
                        firstNeg = charge_1[0] < 0.
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                        if math.isnan(aco): aco = -9999.
                        variablesCP['aco_lep_rhoCollinear'] = RadToDeg*aco

                        # Reco
                        P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'reco')
                        firstNeg = charge_1[0] < 0.
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                        if math.isnan(aco): aco = -9999.
                        variablesCP['aco_lep_rhoReco'] = RadToDeg*aco
                    
                        # Reco-Gen
                        P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'reco_gen')
                        firstNeg = charge_1[0] < 0.
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                        if math.isnan(aco): aco = -9999.
                        variablesCP['aco_lep_rhoRecoGen'] = RadToDeg*aco
                        variablesCP['aco_lep_rho'] = RadToDeg*aco_lep_rho[0]
                        if alpha > pi_over_4:
                            variablesCP['aco_lep_rho_plus'] = RadToDeg*aco_lep_rho[0]
                        else:
                            variablesCP['aco_lep_rho_minus'] = RadToDeg*aco_lep_rho[0]
                    if isTauToA1_3pr_2:
                        # tau->a1(3pr)+v decay
                        PGen = ROOT.TLorentzVector()
                        PGen.SetPtEtaPhiM(genPart_pt_2[0],genPart_eta_2[0],genPart_phi_2[0],utils.tau_mass)
                    
                        sv_pv = ROOT.TVector3(sv_x_2[0]-PVBS_x[0],sv_y_2[0]-PVBS_y[0],sv_z_2[0]-PVBS_z[0])
                        Pi1 = ROOT.TLorentzVector()
                        Pi1.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                        Pi2 = ROOT.TLorentzVector()
                        Pi2.SetPtEtaPhiM(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],massPi)
                        Pi3 = ROOT.TLorentzVector()
                        Pi3.SetPtEtaPhiM(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],massPi)
                        lv_a1 = Pi1 + Pi2 + Pi3
                        cosa = pv_utils.CosAlpha(lv_a1.Pt(),lv_a1.Eta(),lv_a1.Phi(),lv_a1.M(),sv_pv.X(),sv_pv.Y(),sv_pv.Z(),False)
                        alpha = ROOT.TMath.ACos(cosa)
                        variablesCP['alpha_lep_a1'] = RadToDeg * alpha

                        PV_Mag = math.sqrt(PVBS_x[0]*PVBS_x[0]+PVBS_y[0]*PVBS_y[0]+PVBS_z[0]*PVBS_z[0])
                        SV_Mag = math.sqrt(sv_x_2[0]*sv_x_2[0]+sv_y_2[0]*sv_y_2[0]+sv_z_2[0]*sv_z_2[0])
                        P1 = ROOT.TLorentzVector()
                        P1.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],massLep)
                        R1 = ROOT.TLorentzVector()
                        R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0.0)
                        
                        PV = ROOT.TLorentzVector()
                        PV.SetXYZT(PVBS_x[0],PVBS_y[0],PVBS_z[0],PV_Mag)
                        SV = ROOT.TLorentzVector()
                        SV.SetXYZT(sv_x_2[0],sv_y_2[0],sv_z_2[0],SV_Mag)

                        P_os,P_ss1,P_ss2 = pv_utils.sortA1(Pi1,Pi2,Pi3,pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])
                        firstNeg = charge_1[0] < 0.

                        # DESY method
                        P2,R2  = pv_utils.PolVectA1(PV,SV,
                                                    pt_2_FastMTT[0],P_os,P_ss1,P_ss2,PGen,
                                                    charge_2[0],'recoDESY')
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                        if math.isnan(aco): aco = -9999.
                        variablesCP['aco_lep_a1PVDESY'] = RadToDeg * aco
                        aco_DESY = aco
                        
                        variablesCP['aco_lep_a1DP'] = RadToDeg*aco_lep_a1[0]
                        variablesCP['aco_lep_a1'] = RadToDeg*aco_lep_a1_FastMTT[0]
                        if abs(ipsig_1[0])>cuts.ipsigLepCut:
                            variablesCP['aco_lepIP_a1'] = RadToDeg*aco_lep_a1_FastMTT[0]
                        if alpha > pi_over_4:
                            variablesCP['aco_lep_a1_plus'] = RadToDeg*aco_lep_a1_FastMTT[0]
                        else:
                            variablesCP['aco_lep_a1_minus'] = RadToDeg*aco_lep_a1_FastMTT[0]

                
        for hist in hists:
            hists[hist].Scale(self.norm)
        return hists



            
