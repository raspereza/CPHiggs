import ROOT 
import math
from array import array
import numpy as np
import os
import CPHiggs.IP.utils as utils
from CPHiggs.IP.ScaleFactor import ScaleFactor
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1
import CPHiggs.IP.pv_utils as pv_utils

# Cuts for Z->tau+tau and Z->ll selection
class AnalysisCuts:
    def __init__(self,**kwargs):
        self.mtCut = kwargs.get('mtCut',999999.)
        self.mvisCut = kwargs.get('mvisCut',999999.)
        
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
        
        self.antiMu  = kwargs.get('antiMu',4)
        self.antiE   = kwargs.get('antiE',6)
        self.antiJet = kwargs.get('antiJet',5)
        self.useCrossTrigger = kwargs.get('useCrossTrigger',False)

        self.tauToRhoDE = kwargs.get('TauToRhoDE',0.2)
        
        self.isoLepCut = kwargs.get('isoLepCut',0.15)

        self.ipsigLepCut = kwargs.get('ipsigLepCut',1.0)

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
        print("mvisCut",self.mvisCut)
        
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
def RunSamplesTuple(samples,name):

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
        self.sampleName = samplename + '_' + era
        self.sampleFile = ROOT.TFile(filename,"READ")
        self.channel = channel
        self.era = era
        self.norm = norm
        self.ismc = True
        self.isdata = isdata        
        if isdata:
            self.norm = 1.0
            self.ismc = False

        self.sign_labels = ['os','ss']
        self.iso_labels = ['iso','antiiso','rest']
        self.type_labels = ['tau','lep','had','all']
        self.scale_unc = ['lepUp','lepDown','tauUp','tauDown','lfakeUp','lfakeDown']
        
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

        if ptX<self.ptMin: ptX = self.ptMin+0.01
        if ptX>self.ptMax: ptX = self.ptMax-0.01

        if etaX<self.etaMin: etaX = self.etaMin+0.01
        if etaX>self.etaMax: etaX = self.etaMax-0.01

        binPt = self.histPtBins.FindBin(ptX)
        binEta = self.histEtaBins.FindBin(etaX)
        binLabel = '%1i_%1i'%(binPt,binEta)

        #        print(pt,eta,binPt,binEta,binLabel)
        
        return binLabel,binPt,binEta
        
    def DeclareHistos(self):
        
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

        # tag-and-probe histos
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
                                for unc in self.scale_unc:
                                    name = 'm_vis_%s_%s_%s_%s_%s_%s'%(unc,region,label,sign,iso,typ)
                                    histname = self.sampleName+'_'+ name
                                    hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))

#        for hist in hists:
#            print(hist)
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

        print('applyWeightCP = %1i'%(applyWeightCP))
        
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
        
        sv_x_2       = np.zeros(1,dtype=np.float64)
        sv_y_2       = np.zeros(1,dtype=np.float64)
        sv_z_2       = np.zeros(1,dtype=np.float64)

        genPart_pt_2 = np.zeros(1,dtype=np.float64)
        genPart_eta_2 = np.zeros(1,dtype=np.float64)
        genPart_phi_2 =	np.zeros(1,dtype=np.float64)

        alphaAngle_2 = np.zeros(1,dtype=np.float64)

        # general variables
        
        m_vis        = np.zeros(1,dtype=np.float64)
        met_pt       = np.zeros(1,dtype=np.float64)
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

        tree.SetBranchAddress('genPart_pt_2',genPart_pt_2)
        tree.SetBranchAddress('genPart_eta_2',genPart_eta_2)
        tree.SetBranchAddress('genPart_phi_2',genPart_phi_2)
        
        if channel=='mt' or channel=='et' or channel=='tt':
            
            tree.SetBranchAddress('FastMTT_pt_1_constraint',pt_1_FastMTT)
            tree.SetBranchAddress('FastMTT_pt_2_constraint',pt_2_FastMTT)

            tree.SetBranchAddress('PVBS_x',PVBS_x)
            tree.SetBranchAddress('PVBS_y',PVBS_y)
            tree.SetBranchAddress('PVBS_z',PVBS_z)

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

            tree.SetBranchAddress('sv_x_2',sv_x_2)
            tree.SetBranchAddress('sv_y_2',sv_y_2)
            tree.SetBranchAddress('sv_z_2',sv_z_2)

            tree.SetBranchAddress('hasRefitSV_2',hasRefitSV_2)

        if channel=='tt':

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
#            tree.SetBranchAddress('aco_mu_a1_FASTMTT_NoMassConstraint',aco_lep_a1_FastMTT)
            tree.SetBranchAddress('aco_mu_a1_FASTMTT_MassConstraint',aco_lep_a1_FastMTT)


        if channel=='et':
            tree.SetBranchAddress('aco_e_pi',aco_lep_pi)
            tree.SetBranchAddress('aco_e_rho',aco_lep_rho)
            tree.SetBranchAddress('aco_e_a1',aco_lep_a1)
#            tree.SetBranchAddress('aco_e_a1_FASTMTT_NoMassConstraint',aco_lep_a1_FastMTT)
            tree.SetBranchAddress('aco_e_a1_FASTMTT_MassConstraint',aco_lep_a1_FastMTT)

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
            
            tree.SetBranchAddress('aco_pi_a1_FASTMTT_MassConstraint',aco_pi_a1_FastMTT)
            tree.SetBranchAddress('aco_a1_pi,FASTMTT_MassConstraint',aco_a1_pi_FastMTT)
            tree.SetBranchAddress('aco_rho_a1_FASTMTT_MassConstraint',aco_rho_a1_FastMTT)
            tree.SetBranchAddress('aco_a1_rho_FASTMTT_MassConstraint',aco_a1_rho_FastMTT)
            tree.SetBranchAddress('aco_a1_a1_FASTMTT_MassConstraint',aco_a1_a1_FastMTT)

            
        if channel=='mt' or channel=='et' or channel=='tt':
            tree.SetBranchAddress('decayMode_2',decayMode_2)
            tree.SetBranchAddress('decayModePNet_2',decayModePNet_2)

        if channel=='tt':
            tree.SetBranchAddress('decayMode_1',decayMode_1)
            tree.SetBranchAddress('decayModePNet_1',decayModePNet_1)
            
        tree.SetBranchAddress('weight',weight)

        if applyWeightCP==1: # SM (CP-even)
            tree.SetBranchAddress('wt_cp_sm',weight_CP)
        if applyWeightCP==2: # PS (CP-odd)
            tree.SetBranchAddress('wt_cp_ps',weight_CP)
        
        # booleans (trigger)
        if channel=='mt' or channel=='mm': tree.SetBranchAddress('trg_singlemuon',trg_lep)
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
                ptTrig = pt_1[0]
                etaTrig = abs(eta_1[0])
                if pt_2[0]>pt_1[0]:
                    ptTrig = pt_2[0]
                    etaTrig = abs(eta_2[0])

                passTrigger = ptTrig>cuts.ptSingleLepTrigger and etaTrig<cuts.etaSingleLepTrigger and trg_lep[0]
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
                # mT cut
                if mt_1[0]>cuts.mtCut: continue
                if m_vis[0]>cuts.mvisCut: continue
                # tau discriminator against e and mu and jet
                if idDeepTau2018v2p5VSe_2[0]<cuts.antiE: continue
                if idDeepTau2018v2p5VSmu_2[0]<cuts.antiMu: continue
                if idDeepTau2018v2p5VSjet_2[0]<cuts.antiJet: continue

            variables = {} 
            variables['m_vis'] = m_vis[0]
            variables['pt_1'] = pt_1[0]
            variables['eta_1'] = eta_1[0]
            variables['mt_1'] = mt_1[0]
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
            

            if channel=='mt' or channel=='et':
                variables['aco_lep_pi_plus'] = -9999.
                variables['aco_lep_pi_minus'] = -9999.
                variables['aco_lep_pi'] = -9999.
                variables['aco_lep_piIP'] = -9999.
            
                variables['aco_lep_rho_plus'] = -9999.
                variables['aco_lep_rho_minus'] = -9999.
                variables['aco_lep_rho'] = -9999.
                variables['aco_lep_rhoECut'] = -9999.
                variables['aco_lep_rhoGen'] = -9999.
                variables['aco_lep_rhoReco'] = -9999.
                variables['aco_lep_rhoCollinear'] = -9999.
                variables['aco_lep_rhoRecoGen'] = -9999.
                variables['aco_lep_rhoRecoGenECut'] = -9999.
                
                variables['aco_lep_rhoRecoIP1p0'] = -9999.
                variables['aco_lep_rhoRecoGenIP1p0'] = -9999.
                variables['aco_lep_rhoRecoIP1p0ECut'] = -9999.
                variables['aco_lep_rhoRecoGenIP1p0ECut'] = -9999.
                
                variables['aco_lep_rhoRecoIP1p2'] = -9999.
                variables['aco_lep_rhoRecoGenIP1p2'] = -9999.
                variables['aco_lep_rhoRecoIP1p2ECut'] = -9999.
                variables['aco_lep_rhoRecoGenIP1p2ECut'] = -9999.
                
                variables['aco_lep_rhoRecoIP1p5'] = -9999.
                variables['aco_lep_rhoRecoGenIP1p5'] = -9999.
                variables['aco_lep_rhoRecoIP1p5ECut'] = -9999.
                variables['aco_lep_rhoRecoGenIP1p5ECut'] = -9999.
                
                variables['aco_lep_a1_plus'] = -9999.
                variables['aco_lep_a1_minus'] = -9999.
                variables['aco_lep_a1'] = -9999.
                variables['aco_lep_a1DP'] = -9999.
                variables['aco_lep_a1PVGen'] = -9999.
                variables['aco_lep_a1PVDESY'] = -9999.
                variables['aco_lep_a1PVIC'] = -9999.
                
                variables['alpha_lep_pi'] = -9999.
                variables['alpha_lep_rho'] = -9999.
                variables['alpha_lep_a1'] = -9999.

                alpha = -9999.
                cosa = pv_utils.CosAlpha(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi,ip_x_2[0],ip_y_2[0],ip_z_2[0],False)
                alpha_IP = ROOT.TMath.ACos(cosa)
                aco = -9999.
                pv = ROOT.TLorentzVector()
                aco_a1 = -9999.
                
                if decayModePNet_2[0]==0: 
                    alpha = alpha_IP
                    variables['alpha_lep_pi'] = RadToDeg * alpha
                    """
                    P1 = ROOT.TLorentzVector()
                    P1.SetPtEtaPhiM(10.,eta_1[0],phi_1[0],massLep)
                    P2 = ROOT.TLorentzVector()
                    P2.SetPtEtaPhiM(10.,pi_eta_2[0],pi_phi_2[0],massPi)
                    R1 = ROOT.TLorentzVector()
                    R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0)
                    R2 = ROOT.TLorentzVector()
                    R2.SetXYZT(ip_x_2[0],ip_y_2[0],ip_z_2[0],0)
                    firstNeg = charge_1[0] < 0.
                    aco = pv.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','Impact-Parameter')
                    P1.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],massLep)
                    P2.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                    aco_old = pv.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','Impact-Parameter')
                    if math.isnan(aco): aco = -9999.
                    """
                    variables['aco_lep_pi'] = RadToDeg*aco_lep_pi[0]
                    if abs(ipsig_2[0])>cuts.ipsigLepCut:
                        variables['aco_lep_piIP'] = RadToDeg*aco_lep_pi[0]
                    if alpha > pi_over_4:
                        variables['aco_lep_pi_plus'] = RadToDeg*aco_lep_pi[0]
                    else:
                        variables['aco_lep_pi_minus'] = RadToDeg*aco_lep_pi[0]
#                    print('PNetDM %2i -> alpha = %5.3f : %5.3f -- phi(CP) = %5.3f : %5.3f'%(decayModePNet_2[0],alpha,alphaAngle_2[0],aco_old,aco))
#                    print('')
#                elif decayModePNet_2[0]==1 and decayMode_2[0]==1 and pi0_pt_2[0]>1 and pt_2_FastMTT[0]>10. and abs(ipsig_2[0])>cuts.ipsigLepCut:
                elif decayModePNet_2[0]==1 and decayMode_2[0]==1 and pi0_pt_2[0]>1 and pt_2_FastMTT[0]>10.:
                    Pi0 = ROOT.TLorentzVector()
                    Pi0.SetPtEtaPhiM(pi0_pt_2[0],pi0_eta_2[0],pi0_phi_2[0],massPi0)
                    Pi = ROOT.TLorentzVector()
                    Pi.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                    magPi0 = Pi0.Vect().Mag()
                    magPi = Pi.Vect().Mag()
                    cosa = pv_utils.CosAlpha(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi,Pi0.Px(),Pi0.Py(),Pi0.Pz(),True)
                    alpha = ROOT.TMath.ACos(cosa)
                    variables['alpha_lep_rho'] = RadToDeg * alpha
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
                    variables['aco_lep_rhoGen'] = RadToDeg*aco

                    # Collinear
                    P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'collinear')
                    firstNeg = charge_1[0] < 0.
                    aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_lep_rhoCollinear'] = RadToDeg*aco

                    # Reco
                    P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'reco')
                    firstNeg = charge_1[0] < 0.
                    aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_lep_rhoReco'] = RadToDeg*aco
                    if deltaE>cuts.tauToRhoDE:
                        variables['aco_lep_rhoRecoECut'] = RadToDeg*aco
                    if ipsig_2[0]>1.0:
                        variables['aco_lep_rhoRecoIP1p0'] = RadToDeg*aco
                        if deltaE>cuts.tauToRhoDE:
                            variables['aco_lep_rhoRecoIP1p0ECut'] = RadToDeg*aco

                    if ipsig_2[0]>1.2:
                        variables['aco_lep_rhoRecoIP1p2'] = RadToDeg*aco
                        if deltaE>cuts.tauToRhoDE:
                            variables['aco_lep_rhoRecoIP1p2ECut'] = RadToDeg*aco
                    
                    if ipsig_2[0]>1.5:
                        variables['aco_lep_rhoRecoIP1p5']
                        if deltaE>cuts.tauToRhoDE:
                            variables['aco_lep_rhoRecoIP1p5ECut']
                    
                    # Reco-Gen
                    P2,R2 = pv_utils.PolVectRho(pt_2_FastMTT[0],Pi,Pi0,ip_vec_2,PGen,'reco_gen')
                    firstNeg = charge_1[0] < 0.
                    aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_lep_rhoRecoGen'] = RadToDeg*aco
                    if deltaE>cuts.tauToRhoDE:
                        variables['aco_lep_rhoRecoGenECut'] = RadToDeg*aco
                    if ipsig_2[0]>1.0:
                        variables['aco_lep_rhoRecoGenIP1p0'] = RadToDeg*aco
                        if deltaE>cuts.tauToRhoDE:
                            variables['aco_lep_rhoRecoGenIP1p0ECut'] = RadToDeg*aco

                    if ipsig_2[0]>1.2:
                        variables['aco_lep_rhoRecoGenIP1p2'] = RadToDeg*aco
                        if deltaE>cuts.tauToRhoDE:
                            variables['aco_lep_rhoRecoGenIP1p2ECut'] = RadToDeg*aco
                    
                    if ipsig_2[0]>1.5:
                        variables['aco_lep_rhoRecoGenIP1p5']
                        if deltaE>cuts.tauToRhoDE:
                            variables['aco_lep_rhoRecoGenIP1p5ECut']
                    

                    
                    # Reco-Gen
                        
                    variables['aco_lep_rho'] = RadToDeg*aco_lep_rho[0]
                    if deltaE>cuts.tauToRhoDE:
                        variables['aco_lep_rhoECut'] = RadToDeg*aco_lep_rho[0]
                        
                    if alpha > pi_over_4:
                        variables['aco_lep_rho_plus'] = RadToDeg*aco_lep_rho[0]
                    else:
                        variables['aco_lep_rho_minus'] = RadToDeg*aco_lep_rho[0]
#                    print('PNetDM %2i -> alpha = %5.3f : %5.3f -- phi(CP) = %5.3f : %5.3f'%(decayModePNet_2[0],alpha_IP,alphaAngle_2[0],aco_lep_rho[0],aco))
#                    print('')
                elif decayModePNet_2[0]==10 and hasRefitSV_2[0] and pt_2_FastMTT[0]>10.:
#                    print('run = %1i   event = %1i'%(run[0],event[0]))
#                    print('PNetDM = %2i  hasRefitSV_2 = %1i'%(decayModePNet_2[0],hasRefitSV_2[0]))
#                    print('muon   (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pt_1[0],eta_1[0],phi_1[0],massLep,charge_1[0]))
#                    print('muon impact par. (x,y,z) = (%8.6f,%8.6f,%8.6f)'%(ip_x_1[0],ip_y_1[0],ip_z_1[0]))
#                    print('FastMTT_pt_2_constraint  = %5.1f'%(pt_2_FastMTT[0]))
#                    print('pi1    (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],pi_mass_2[0],int(pi_charge_2[0])))
#                    print('pi2    (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],pi2_mass_2[0],int(pi2_charge_2[0])))
#                    print('pi3    (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],pi3_mass_2[0],int(pi3_charge_2[0])))

                    
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
                    variables['alpha_lep_a1'] = RadToDeg * alpha

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

#                    print('PV (x,y,z) = (%6.5f,%6.5f,%6.5f)'%(PV.X(),PV.Y(),PV.Z()))
#                    print('SV (x,y,z) = (%6.5f,%6.5f,%6.5f)'%(SV.X(),SV.Y(),SV.Z()))
                    
                    P_os,P_ss1,P_ss2 = pv_utils.sortA1(Pi1,Pi2,Pi3,pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])
                    firstNeg = charge_1[0] < 0.

                    # Gen
                    P2,R2 = pv_utils.PolVectA1(PV,SV,
                                               pt_2_FastMTT[0],P_os,P_ss1,P_ss2,PGen,
                                               charge_2[0],'gen')
                    aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_lep_a1PVGen'] = RadToDeg * aco

                    # IC method
                    P2,R2 = pv_utils.PolVectA1(PV,SV,
                                               pt_2_FastMTT[0],P_os,P_ss1,P_ss2,PGen,
                                               charge_2[0],'recoIC')
                    aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_lep_a1PVIC'] = RadToDeg * aco

                    # DESY method
                    P2,R2  = pv_utils.PolVectA1(PV,SV,
                                                pt_2_FastMTT[0],P_os,P_ss1,P_ss2,PGen,
                                                charge_2[0],'recoDESY')
                    aco = pv_utils.acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_lep_a1PVIC'] = RadToDeg * aco

                    variables['aco_lep_a1DP'] = RadToDeg*aco_lep_a1[0]
                    variables['aco_lep_a1'] = RadToDeg*aco_lep_a1_FastMTT[0]
                    if alpha > pi_over_4:
                        variables['aco_lep_a1_plus'] = RadToDeg*aco_lep_a1_FastMTT[0]
                    else:
                        variables['aco_lep_a1_minus'] = RadToDeg*aco_lep_a1_FastMTT[0]

#                    print('aco_mu_a1                        = %5.3f (tuple) : %5.3f (my code)'%(aco_lep_a1[0],aco_a1))
#                    print('aco_mu_a1_FASTMTT_MassConstraint = %5.3f (tuple) : %5.3f (IC code) : %5.3f (DESY code)'%(aco_lep_a1_FastMTT[0],aco_IC,aco))
#                    print('')

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
                    
            ################################
            ## Filling control histograms ##
            ################################
            for varname in variables:
                nameAll = '%s_%s_%s_all'%(varname,sign_label,iso_label)
                name = '%s_%s_%s_%s'%(varname,sign_label,iso_label,lep_label)
                hists[nameAll].Fill(variables[varname],Weight*WeightSF)
                hists[name].Fill(variables[varname],Weight*WeightSF)

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
                name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep2_label)
                hists[name].Fill(variables['m_vis'],Weight)
                name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
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
                    mvis_unc['tauUp'] =	variables['m_vis']
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
                
        for hist in hists:
            hists[hist].Scale(self.norm)
        return hists
