import ROOT 
import math
from array import array
import numpy as np
import os
import CPHiggs.Analysis.utils as utils
from CPHiggs.Analysis.ScaleFactor import ScaleFactor
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1
import CPHiggs.Analysis.pv_utils as pv_utils

def cosDeltaPhi(x1,x2):
    px1 = x1.Px()
    px2 = x2.Px()
    py1 = x1.Py()
    py2 = x2.Py()
    pt1 = x1.Pt()
    pt2 = x2.Pt()
    num = px1*px2+py1*py2
    den = pt1*pt2
    if den<0.001: den=0.001
    cosinus = num/den
    return cosinus

# Cuts for Z->tau+tau and Z->ll selection
class AnalysisCuts:
    def __init__(self,**kwargs):
        
        self.mtCut = kwargs.get('mtCut',65.)
        self.mvisUpperCut = kwargs.get('mvisUpperCut',999999.)
        self.mvisLowerCut = kwargs.get('mvisLowerCut',0.)
        
        self.etaLep1Cut = kwargs.get('etaLep1Cut',2.4)
        self.etaLep2Cut = kwargs.get('etaLep2Cut',2.5)
        
        self.ptLep1Cut = kwargs.get('ptLep1Cut',21.)
        self.ptLep2Cut = kwargs.get('ptLep2Cut',20.)

        self.ptSingleLepTrigger  = kwargs.get('ptSingleLepTrigger',26.) # 31 for single-e
        self.etaSingleLepTrigger = kwargs.get('etaSingleLepTrigger',2.4) # 2.1 for single-e

        self.ptLepCrossTrigger  = kwargs.get('ptLepCrossTrigger',21.)  # 25 for e-tau
        self.etaLepCrossTrigger = kwargs.get('etaLepCrossTrigger',2.1) # 2.1 for e-tau
        self.ptTauCrossTrigger  = kwargs.get('ptTauCrossTrigger',32.)  # 35 for e-tau
        self.etaTauCrossTrigger = kwargs.get('etaTauCrossTrigger',2.1) # 2.1 for e-tau

        self.ptDiTauTrigger = kwargs.get('ptDitauTrigger',40.) # 
        
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
        self.antiJetInverted = kwargs.get('antiJetInverted',4)
        self.useCrossTrigger = kwargs.get('useCrossTrigger',False)

        self.tauToRhoDE = kwargs.get('TauToRhoDE',0.2)
        
        self.isoLepCut = kwargs.get('isoLepCut',0.15)

        self.ipsigLepCut = kwargs.get('ipsigLepCut',1.0)
        self.ipsigTauCut = kwargs.get('ipsigTauCut',1.25)
        
        self.isoLepInverseLowerCut = kwargs.get('isoLepInverseLowerCut',0.05)
        self.isoLepInverseUpperCut = kwargs.get('isoLepInverseUpperCut',1.00)

        self.applyBVeto = kwargs.get('applyBVeto',False)
        
        self.applyIPSigLep1Cut = kwargs.get('applyIPSigLep1Cut',False)
        self.applyIPSigLep2Cut = kwargs.get('applyIPSigLep2Cut',False)

        self.lepMomScale = kwargs.get('lepMomScale',0.002) # lepton mom scale unc. 0.002/0.025 for muons/electrons
        self.tauMomScale = kwargs.get('tauMomScale',0.02)  # tau mom scale unc. 
        self.lepTauMomScale = kwargs.get('lepTauMomScale',0.04) # lep->tau fake mom scale unc. 0.02/0.04 for muons/electrons

        self.signalCat = kwargs.get('SignalCategory',1) # index of signal category in BDT
        self.signalRegion = kwargs.get('SignalRegion',False) # select to signal category
        
        print('')
        print("Setting cuts ->")

        self.mtUpperCut = 65.
        self.mtLowerCut = 70.
        self.bdtWCut = 0.2
        
        print("mtCut",self.mtCut)
        print("mtUpperCut",self.mtUpperCut)
        print("mtLowerCut",self.mtLowerCut)
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
        print("antiJetInverted",self.antiJetInverted)

        print("useCrossTrigger",self.useCrossTrigger)
        print("BVeto",self.applyBVeto)

        print("ipsigLepCut",self.ipsigLepCut)
        print("ipsigTauCut",self.ipsigTauCut)
        print("applyIPSigLep1Cut",self.applyIPSigLep1Cut)
        print("applyIPSigLep2Cut",self.applyIPSigLep2Cut)
        print("TauToRhoDE",self.tauToRhoDE)
        
        print("lepMomScale",self.lepMomScale)
        print("tauMomScale",self.tauMomScale)
        print("lepTauMomScale",self.lepTauMomScale)
        print("signalCat",self.signalCat)

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
        if self.analysisType not in ['baseline','ipSig','phiCP','jetFakes','datacardsPhiCP']:
            print('')
            print('Unknown analysis type : %s'%(self.analysisType))
            print('')
            exit()
        if isdata:
            self.norm = 1.0
            self.ismc = False

        
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

        self.applyIPSigCorr  = kwargs.get('applyIPSigCorr',False)
        self.applyFakeFactor = kwargs.get('applyFakeFactor',False)
        
        self.ipSigCorr  = kwargs.get('ipSigCorr',None)
        self.fakeFactor = kwargs.get('fakeFactor',None)
        
        self.applyWeightCP = kwargs.get('applyWeightCP',False)
        self.cppKinFit = kwargs.get('cppKinFit',False)
        self.phiScan = kwargs.get('phiScan',False)
        
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
        ########################################
        ######### baseline analysis ############
        ########################################
        if analysisType=='baseline':
            for var in utils.lib_histos:
            
                nbins = utils.lib_histos[var][0]
                xmin  = utils.lib_histos[var][1]
                xmax  = utils.lib_histos[var][2]
                xbins = utils.createBins(nbins,xmin,xmax)

                for region in utils.region_labels:
                    for typ in utils.type_labels:
                        name = '%s_%s_%s'%(var,region,typ)
                        histname = self.sampleName+'_'+name
                        hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                        for ff in utils.ff_labels:
                            name = '%s_%s_%s_%s'%(var,region,ff,typ)
                            histname = self.sampleName+'_'+name
                            hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                """
                        for njets in utils.njets_labels:
                            name = '%s_%s_%s_%s'%(var,region,njets,typ)
                            histname = self.sampleName+'_'+name
                            hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                            for ff in utils.ff_labels:
                                name = '%s_%s_%s_%s_%s'%(var,region,njets,ff,typ)
                                histname = self.sampleName+'_'+name
                                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                """             

        ######################################
        ########## phi(CP) studies ###########
        ######################################
        if analysisType=='phiCP':
            for var in utils.lib_phiCP_histos:

                nbins = utils.lib_phiCP_histos[var][0]
                xmin  = utils.lib_phiCP_histos[var][1]
                xmax  = utils.lib_phiCP_histos[var][2]
                xbins = utils.createBins(nbins,xmin,xmax)
                
                name = '%s'%(var)
                histname = self.sampleName+'_'+name
                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                if self.applyWeightCP:
                    for hypothesis in utils.cp_hypotheses:
                        name = '%s_%s'%(var,hypothesis)
                        histname = self.sampleName+'_'+name
                        hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))


        ######################################################
        ## measurements of SFs related to cut on mu/e ipSig ##
        ######################################################
        if analysisType=='ipSig':
            
            nbins = utils.lib_histos['m_vis'][0]
            xmin  = utils.lib_histos['m_vis'][1]
            xmax  = utils.lib_histos['m_vis'][2]
            xbins = utils.createBins(nbins,xmin,xmax)

            for sign in utils.sign_labels:
                for iso in utils.iso_labels:
                    for typ in utils.type_labels:
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


        ##########################################
        #######  datacard production #############
        ##########################################
        if analysisType=='datacardsPhiCP':

            # BDT variable : binning
            nbins_bdt = 100
            xmin_bdt = 0.0
            xmax_bdt = 1.0
            xbins_bdt = utils.createBins(nbins_bdt,xmin_bdt,xmax_bdt)
            
            # phi(CP) variable binning
            nbins_phicp = 40
            xmin_phicp = 0.
            xmax_phicp = 360.
            xbins_phicp =  utils.createBins(nbins_phicp,xmin_phicp,xmax_phicp)
            
            for var in utils.lib_datacards_1Dhistos:
                for sign in utils.sign_labels:
                    for iso in utils.iso_labels:
                        for typelep in utils.type_labels:
                            for typetau in utils.type_labels:
                                name = '%s_%s_%s_%s_%s'%(var,sign,iso,typelep,typetau)
                                histname = self.sampleName+'_'+name
                                hists[name] = ROOT.TH1D(histname,"",nbins_bdt,array('d',list(xbins_bdt)))
                                if self.applyWeightCP:
                                    for hypothesis in utils.cp_hypotheses:
                                        name = '%s_%s_%s_%s_%s_%s'%(var,sign,iso,typelep,typetau,hypothesis)
                                        histname = self.sampleName+'_'+name
                                        hists[name] = ROOT.TH1D(histname,"",nbins_bdt,array('d',list(xbins_bdt)))
                                
            for dm in utils.dm_labels:
                decay_label = 'lep_%s'%(dm)
                var = 'phicp_vs_bdt_%s'%(decay_label)
                for sign in utils.sign_labels:
                    for iso in utils.iso_labels:
                        for typelep in utils.type_labels:
                            for typetau in utils.type_labels:
                                name = '%s_%s_%s_%s_%s'%(var,sign,iso,typelep,typetau)
                                histname = self.sampleName+'_'+name
                                hists[name] = ROOT.TH2D(histname,"",
                                                        nbins_bdt,array('d',list(xbins_bdt)),
                                                        nbins_phicp,array('d',list(xbins_phicp)))
                                if self.applyWeightCP:
                                    for hypothesis in utils.cp_hypotheses:
                                        name = '%s_%s_%s_%s_%s_%s'%(var,sign,iso,typelep,typetau,hypothesis)
                                        histname = self.sampleName+'_'+name
                                        hists[name] = ROOT.TH2D(histname,"",
                                                                nbins_bdt,array('d',list(xbins_bdt)),
                                                                nbins_phicp,array('d',list(xbins_phicp)))
                                            

        ####################################
        ####### jetFakes modeling ##########
        ####################################
        if analysisType=='jetFakes':
            nbins = utils.lib_histos['pt_2'][0]
            xmin  = utils.lib_histos['pt_2'][1]
            xmax  = utils.lib_histos['pt_2'][2]
            xbins = utils.createBins(nbins,xmin,xmax)
            for dm in utils.dm_labels:
                for njets in utils.njets_labels:
                    for eta in utils.eta_labels:
                        for tauid in utils.tauid_labels:
                            for typ in utils.type_labels:
                                # QCD : ss, antiiso, lowmt
                                name = 'pt_2_qcd_%s_%s_%s_%s_%s'%(dm,njets,eta,tauid,typ)
                                histname = self.sampleName+'_'+ name
                                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                                # WJets : ss, iso, highmt
                                name = 'pt_2_wj_%s_%s_%s_%s_%s'%(dm,njets,eta,tauid,typ)
                                histname = self.sampleName+'_'+ name
                                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                                name = 'pt_2_top_%s_%s_%s_%s_%s'%(dm,njets,eta,tauid,typ)
                                histname = self.sampleName+'_'+ name
                                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
                                name = 'pt_2_qcd_looseIso_%s_%s_%s_%s_%s'%(dm,njets,eta,tauid,typ)
                                histname = self.sampleName+'_'+ name
                                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))
#        for hist in hists:
#            print(hist)
                                        
        return hists

    def CreateHistosTuple(self):

        tree = self.sampleFile.Get('ntuple')

        # initialization
        cuts = self.analysisCuts
        
        applyIPSigCorr = self.applyIPSigCorr
        ipSigCorr = self.ipSigCorr
        applyWeightCP = self.applyWeightCP
        analysisType = self.analysisType
        fakeFactor = self.fakeFactor
        applyFakeFactor = self.applyFakeFactor
        cppKinFit = self.cppKinFit
        phiScan = self.phiScan

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
        if channel=='mt' or channel=='mm':
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
        w_IP_Significance = np.zeros(1,dtype=np.float64)

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
        pi_phi_1    = np.zeros(1,dtype=np.float64)
        pi_mass_1   = np.zeros(1,dtype=np.float64)
        pi_charge_1 = np.zeros(1,dtype=np.float64)
        
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

        sv_cov00_1   = np.zeros(1,dtype=np.float64)
        sv_cov10_1   = np.zeros(1,dtype=np.float64)
        sv_cov11_1   = np.zeros(1,dtype=np.float64)
        sv_cov20_1   = np.zeros(1,dtype=np.float64)
        sv_cov21_1   = np.zeros(1,dtype=np.float64)
        sv_cov22_1   = np.zeros(1,dtype=np.float64)

        PVBS_cov00   = np.zeros(1,dtype=np.float64)
        PVBS_cov10   = np.zeros(1,dtype=np.float64)
        PVBS_cov20   = np.zeros(1,dtype=np.float64)
        PVBS_cov11   = np.zeros(1,dtype=np.float64)
        PVBS_cov21   = np.zeros(1,dtype=np.float64)
        PVBS_cov22   = np.zeros(1,dtype=np.float64)
        
        genPart_pt_1 = np.zeros(1,dtype=np.float64)
        genPart_eta_1 = np.zeros(1,dtype=np.float64)
        genPart_phi_1 =	np.zeros(1,dtype=np.float64)

        genPart_pt_2 = np.zeros(1,dtype=np.float64)
        genPart_eta_2 = np.zeros(1,dtype=np.float64)
        genPart_phi_2 =	np.zeros(1,dtype=np.float64)

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

        sv_cov00_2   = np.zeros(1,dtype=np.float64)
        sv_cov10_2   = np.zeros(1,dtype=np.float64)
        sv_cov11_2   = np.zeros(1,dtype=np.float64)
        sv_cov20_2   = np.zeros(1,dtype=np.float64)
        sv_cov21_2   = np.zeros(1,dtype=np.float64)
        sv_cov22_2   = np.zeros(1,dtype=np.float64)
        
        
        alphaAngle_2 = np.zeros(1,dtype=np.float64)

        # general variables
        
        m_vis        = np.zeros(1,dtype=np.float64)
        FastMTT_mass = np.zeros(1,dtype=np.float64)
        FastMTT_pt   = np.zeros(1,dtype=np.float64)
        met_pt       = np.zeros(1,dtype=np.float64)
        met_phi      = np.zeros(1,dtype=np.float64)
        met_covXX    = np.zeros(1,dtype=np.float64)
        met_covXY    = np.zeros(1,dtype=np.float64)
        met_covYY    = np.zeros(1,dtype=np.float64)
        aco_lep_pi   = np.zeros(1,dtype=np.float64)
        aco_lep_rho  = np.zeros(1,dtype=np.float64)
        aco_lep_a1   = np.zeros(1,dtype=np.float64)
        weight_cp_sm = np.zeros(1,dtype=np.float64)
        weight_cp_ps = np.zeros(1,dtype=np.float64)
        weight_cp_mm = np.zeros(1,dtype=np.float64)
        
        n_jets       = np.zeros(1,dtype=np.float64)
        n_bjets      = np.zeros(1,dtype=np.float64)
        jdeta        = np.zeros(1,dtype=np.float64)
        mjj          = np.zeros(1,dtype=np.float64)
        jpt_1        = np.zeros(1,dtype=np.float64)
        jpt_2        = np.zeros(1,dtype=np.float64)
        jeta_1       = np.zeros(1,dtype=np.float64)
        jeta_2       = np.zeros(1,dtype=np.float64)

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


        BDT_W_score = np.zeros(1,dtype=np.float64)
        dR = np.zeros(1,dtype=np.float64)
        
        # booleans
        trg_lep       = np.zeros(1,dtype='?')
        trg_lep2      = np.zeros(1,dtype='?')
        trg_cross     = np.zeros(1,dtype='?')
        trg_doubletau = np.zeros(1,dtype='?')
        os            = np.zeros(1,dtype='?')
        
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
        
        # BDTs
        bdt_pred = np.zeros(1,dtype=np.float64)
        class_pred = np.zeros(1,dtype=np.float64)
        
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
        tree.SetBranchAddress('FastMTT_mass',FastMTT_mass)
        tree.SetBranchAddress('FastMTT_pt',FastMTT_pt)
        if not analysisType=='phiCP':
            if self.ismc: tree.SetBranchAddress('w_IP_Significance',w_IP_Significance)

        tree.SetBranchAddress('n_jets',n_jets)
        tree.SetBranchAddress('n_bjets',n_bjets)
        tree.SetBranchAddress('jpt_1',jpt_1)
        tree.SetBranchAddress('jpt_2',jpt_2)
        tree.SetBranchAddress('jeta_1',jeta_1)
        tree.SetBranchAddress('jeta_2',jeta_2)
        tree.SetBranchAddress('jdeta',jdeta)
        tree.SetBranchAddress('mjj',mjj)

        tree.SetBranchAddress('mt_1',mt_1)
        tree.SetBranchAddress('pt_1',pt_1)
        tree.SetBranchAddress('eta_1',eta_1)
        tree.SetBranchAddress('phi_1',phi_1)
        tree.SetBranchAddress('mass_1',mass_1)
        tree.SetBranchAddress('charge_1',charge_1)
        if channel=='mm' or channel=='ee' or channel=='mt' or channel=='et':
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

                tree.SetBranchAddress('genPart_pt_1',genPart_pt_1)
                tree.SetBranchAddress('genPart_eta_1',genPart_eta_1)
                tree.SetBranchAddress('genPart_phi_1',genPart_phi_1)
                
                tree.SetBranchAddress('genPart_pt_2',genPart_pt_2)
                tree.SetBranchAddress('genPart_eta_2',genPart_eta_2)
                tree.SetBranchAddress('genPart_phi_2',genPart_phi_2)

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

                tree.SetBranchAddress('PVBS_cov00',PVBS_cov00)
                tree.SetBranchAddress('PVBS_cov10',PVBS_cov10)
                tree.SetBranchAddress('PVBS_cov20',PVBS_cov20)
                tree.SetBranchAddress('PVBS_cov11',PVBS_cov11)
                tree.SetBranchAddress('PVBS_cov21',PVBS_cov21)
                tree.SetBranchAddress('PVBS_cov22',PVBS_cov22)
                
                tree.SetBranchAddress('sv_cov00_2',sv_cov00_2)
                tree.SetBranchAddress('sv_cov10_2',sv_cov10_2)
                tree.SetBranchAddress('sv_cov11_2',sv_cov11_2)
                tree.SetBranchAddress('sv_cov20_2',sv_cov20_2)
                tree.SetBranchAddress('sv_cov21_2',sv_cov21_2)
                tree.SetBranchAddress('sv_cov22_2',sv_cov22_2)
                
            tree.SetBranchAddress('FastMTT_pt_1_constraint',pt_1_FastMTT)
            tree.SetBranchAddress('FastMTT_pt_2_constraint',pt_2_FastMTT)

            tree.SetBranchAddress('pion_E_split_2',pion_E_split_2)
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

                tree.SetBranchAddress('sv_cov00_1',sv_cov00_1)
                tree.SetBranchAddress('sv_cov10_1',sv_cov10_1)
                tree.SetBranchAddress('sv_cov11_1',sv_cov11_1)
                tree.SetBranchAddress('sv_cov20_1',sv_cov20_1)
                tree.SetBranchAddress('sv_cov21_1',sv_cov21_1)
                tree.SetBranchAddress('sv_cov22_1',sv_cov22_1)


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
            tree.SetBranchAddress('aco_mu_a1_FASTMTT_MassConstraint',aco_lep_a1)


        if channel=='et':
            tree.SetBranchAddress('aco_e_pi',aco_lep_pi)
            tree.SetBranchAddress('aco_e_rho',aco_lep_rho)
            tree.SetBranchAddress('aco_e_a1_FASTMTT_MassConstraint',aco_lep_a1)

        if channel=='tt':
            tree.SetBranchAddress('aco_pi_pi',aco_pi_pi)
            tree.SetBranchAddress('aco_rho_rho',aco_rho_rho)
            tree.SetBranchAddress('aco_a1_a1',aco_a1_a1)

            tree.SetBranchAddress('aco_pi_rho',aco_pi_rho)
            tree.SetBranchAddress('aco_rho_pi',aco_rho_pi)
            
            tree.SetBranchAddress('aco_pi_a1_FASTMTT_MassConstraint',aco_pi_a1)
            tree.SetBranchAddress('aco_a1_pi_FASTMTT_MassConstraint',aco_a1_pi)
            
            tree.SetBranchAddress('aco_rho_a1_FASTMTT_MassConstraint',aco_rho_a1)            
            tree.SetBranchAddress('aco_a1_rho_FASTMTT_MassConstraint',aco_a1_rho)

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

        if applyWeightCP:
            tree.SetBranchAddress('wt_cp_sm',weight_cp_sm)
            tree.SetBranchAddress('wt_cp_ps',weight_cp_ps)
            tree.SetBranchAddress('wt_cp_mm',weight_cp_mm)
            
        if analysisType=='datacardsPhiCP' or analysisType=='baseline':
            tree.SetBranchAddress('BDT_pred_score',bdt_pred)
            tree.SetBranchAddress('BDT_pred_class',class_pred)

        if analysisType=='baseline':
            tree.SetBranchAddress('BDT_W_score',BDT_W_score)
            
        tree.SetBranchAddress('dR',dR)
        
        # booleans (trigger)
        # muons
        if channel=='mt' or channel=='mm': tree.SetBranchAddress('trg_singlemuon',trg_lep)
        if channel=='mm': tree.SetBranchAddress('trg_singlemuon_2',trg_lep2)
        if channel=='mt': tree.SetBranchAddress('trg_mt_cross',trg_cross)
        # electrons
        if channel=='et' or channel=='ee': tree.SetBranchAddress('trg_singleelectron',trg_lep)
        if channel=='ee': tree.SetBranchAddress('trg_singleelectron_2',trg_lep2)
        if channel=='et': tree.SetBranchAddress('trg_et_cross',trg_cross)
        # di-taus
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

            if entry%10000==0 and entry>0:
                print('processed %1i out of %1i events'%(entry,nentries))
            tree.GetEntry(entry)
            
            # trigger threshold
            passTrigger = False
            if channel=='mm' or channel=='ee':
                trg1 = trg_lep[0]
                trg2 = trg_lep2[0]
                pt1Trig = pt_1[0]
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

            if cuts.applyBVeto:
                if n_bjets[0]>0: continue
            
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

            isTauToPi_2 = False
            isTauToRho_2 = False
            isTauToA1_1pr_2 = False
            isTauToA1_3pr_2 = False
            
            if channel=='mt' or channel=='et':
                # m_vis cut
                if m_vis[0]>cuts.mvisUpperCut: continue
                if m_vis[0]<cuts.mvisLowerCut: continue
                # tau discriminator against e and mu and jet
                if idDeepTau2018v2p5VSe_2[0]<cuts.antiE: continue
                if idDeepTau2018v2p5VSmu_2[0]<cuts.antiMu: continue
                if idDeepTau2018v2p5VSjet_2[0]<1: continue

                ##### Decay mode specific selection #####
                isTauToPi_2 = decayModePNet_2[0]==0 and abs(ipsig_2[0])>cuts.ipsigTauCut        
                isTauToRho_2 = decayModePNet_2[0]==1 and decayMode_2[0]==1 and abs(pion_E_split_2[0])>cuts.tauToRhoDE
                isTauToA1_1pr_2 = decayModePNet_2[0]==2 and decayMode_2[0]==1 and abs(pion_E_split_2[0])>cuts.tauToRhoDE
                isTauToA1_3pr_2 = decayModePNet_2[0]==10 and hasRefitSV_2[0]

                passTau_2 =  isTauToPi_2 or isTauToRho_2 or isTauToA1_1pr_2 or isTauToA1_3pr_2
                if not passTau_2:
                    continue
            
            isTauToPi_1 = False
            isTauToRho_1 = False
            isTauToA1_1pr_1 = False
            isTauToA1_3pr_1 = False
            if channel=='tt':
                # m_vis cut
                if m_vis[0]>cuts.mvisUpperCut: continue
                if m_vis[0]<cuts.mvisLowerCut: continue

                # tau discriminator against e and mu and jet
                if idDeepTau2018v2p5VSe_1[0]<cuts.antiE: continue
                if idDeepTau2018v2p5VSmu_1[0]<cuts.antiMu: continue
                if idDeepTau2018v2p5VSjet_1[0]<cuts.antiJet: continue

                # decay modes of the first tau
                isTauToPi_1 = decayModePNet_1[0]==0 and abs(ipsig_1[0])>cuts.ipsigTauCut        
                isTauToRho_1 = decayModePNet_1[0]==1 and decayMode_1[0]==1 and abs(pion_E_split_1[0])>cuts.tauToRhoDE
                isTauToA1_1pr_1 = decayModePNet_1[0]==2 and decayMode_1[0]==1 and abs(pion_E_split_1[0])>cuts.tauToRhoDE
                isTauToA1_3pr_1 = decayModePNet_1[0]==10 and hasRefitSV_1[0]

                passTau_1 =  isTauToPi_1 or isTauToRho_1 or isTauToA1_1pr_1 or isTauToA1_3pr_1 
                if not passTau_1: continue
                
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
                if not passTau_2:
                    continue


            #############################################
            #### event labeling #####
            #############################################
            sign_label = 'ss'
            lep_label = 'had'
            iso_label = 'antiiso'
            lep2_label = 'had'
            mt_label = 'lowmt'
            deeptau_label = 'rest'
            njets_label = 'njets0'
            dm_label = 'rest'
            eta_label = 'barrel'
            
            if mt_1[0]>cuts.mtCut: mt_label = 'highmt'

            absEta = ROOT.TMath.Abs(eta_2[0])
            if absEta>1.48:
                eta_label='endcap'
                
            if n_jets[0]==0: njets_label='njets0'
            elif n_jets[0]==1: njets_label='njets1'
            else: njets_label='njets2'

            DecayMode = -1
            if channel=='mt' or channel=='et':
                if isTauToPi_2:
                    dm_label='pi'
                    DecayMode = 0
                if isTauToRho_2:
                    dm_label='rho'
                    DecayMode = 1
                if isTauToA1_1pr_2:
                    dm_label='a1_1pr'
                    DecayMode = 2
                if isTauToA1_3pr_2:
                    dm_label='a1_3pr'
                    DecayMode = 10
            
                if idDeepTau2018v2p5VSjet_2[0]<cuts.antiJet:
                    deeptau_label = 'inverted'
                if idDeepTau2018v2p5VSjet_2[0]>=cuts.antiJet:
                    deeptau_label = 'nominal'

                if dm_label=='rest': continue
                if deeptau_label=='rest': continue
                
                
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
            invertedIso = iso_1[0] > cuts.isoLepInverseLowerCut and iso_1[0] < cuts.isoLepInverseUpperCut

            iso_loose_ff = iso_1[0] > 0.05 
            iso_ff = iso_1[0] > 0.05 and iso_1[0] < 0.2
            iso_closure_ff = iso_1[0] > 0.2
            
            if directIso: iso_label='iso'

            #####################
            # defining regions ->
            #####################
            region_flags={}
            for reg in utils.region_labels:
                region_flags[reg] = False

            acceptEvent = False
                
            if mt_label=='lowmt' and sign_label=='os'and iso_label=='iso':
                region_flags['lowmt_os_iso'] = True
                acceptEvent = True
            if mt_label=='lowmt' and sign_label=='ss'and iso_label=='iso':
                region_flags['lowmt_ss_iso'] = True
                acceptEvent = True

            if mt_label=='lowmt' and sign_label=='ss'and iso_closure_ff:
                region_flags['qcd_closure_ff'] = True
                acceptEvent = True

            qcd_ff_region = mt_1[0]<cuts.mtUpperCut and sign_label=='ss' and iso_ff 
            if qcd_ff_region:
                region_flags['qcd_ff'] = True
                acceptEvent = True

            qcd_ff_looseIso_region = mt_1[0]<cuts.mtUpperCut and sign_label=='ss' and iso_loose_ff
            if qcd_ff_looseIso_region:
                region_flags['qcd_looseIso_ff'] = True
                acceptEvent = True

            wj_ff_region = mt_1[0]>cuts.mtLowerCut and sign_label=='os' and directIso and n_bjets[0]==0 and BDT_W_score[0]>cuts.bdtWCut
            if wj_ff_region:
                region_flags['wj_ff'] = True
                acceptEvent = True
            
            if (channel=='mt' or channel=='et'):
                if not acceptEvent: continue
            
            aco_lep = {
                'pi'  : aco_lep_pi[0],
                'rho' : aco_lep_rho[0],
                'a1_1pr' : aco_lep_rho[0],
                'a1_3pr' : aco_lep_a1[0],
            }


            #################################
            #### Definition of variables ####
            #################################
            
            lep1LV = ROOT.TLorentzVector()
            lep2LV = ROOT.TLorentzVector()
            metLV = ROOT.TLorentzVector()
            lep1LV.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],massLep)
            if channel=='mm' or channel=='ee':
                lep2LV.SetPtEtaPhiM(pt_2[0],eta_2[0],phi_2[0],massLep)
            else:
                lep2LV.SetPtEtaPhiM(pt_2[0],eta_2[0],phi_2[0],mass_2[0])
            metLV.SetPtEtaPhiM(met_pt[0]*ROOT.TMath.Cos(met_phi[0]),met_pt[0]*ROOT.TMath.Sin(met_phi[0]),0.,0.)

            ptLV = lep1LV+lep2LV+metLV
            pt_tt = ptLV.Pt()

            cosDeltaPhiQCD = cosDeltaPhi(metLV,lep2LV)
            CMetQCD = metLV.Pt()*cosDeltaPhiQCD/lep2LV.Pt()
            metLepLV = metLV+lep1LV
            cosDeltaPhiW = cosDeltaPhi(metLepLV,lep2LV)
            CMetW = metLepLV.Pt()*cosDeltaPhiW/lep2LV.Pt()
            
            variables = {} 
            variables['m_vis'] = m_vis[0]
            variables['m_FastMTT'] = FastMTT_mass[0]
            variables['pt_FastMTT'] = FastMTT_pt[0]
            variables['pt_tt'] = pt_tt 
            variables['pt_1'] = pt_1[0]
            variables['eta_1'] = eta_1[0]
            variables['pt_2'] = pt_2[0]
            variables['eta_2'] = eta_2[0]
            variables['mt_1'] = mt_1[0]
            variables['met'] = met_pt[0]
            variables['CMetQCD'] = CMetQCD
            variables['CMetW'] = CMetW
#            variables['ipsig_1'] = abs(ipsig_1[0])
#            variables['ipsig_2'] = abs(ipsig_2[0])
            variables['dm_2'] = DecayMode
            variables['dR'] = dR[0]

#            print('dR = %5.3f'%(dR[0]))

            # jet related variables
            variables['n_jets'] = n_jets[0]
            variables['n_bjets'] = n_bjets[0]
            variables['jpt_1'] = -9999.
            variables['jpt_2'] = -9999.
            variables['jeta_1'] = -9999.
            variables['jeta_2'] = -9999.
            variables['mjj'] = -9999.
            variables['jdeta'] = -9999.
            if n_jets[0]>0.5:
                variables['jpt_1'] = jpt_1[0]
                variables['jeta_1'] = jeta_1[0]
            if n_jets[0]>1.5:
                variables['jpt_2'] = jpt_2[0]
                variables['jeta_2'] = jeta_2[0]
                variables['mjj'] = mjj[0]
                variables['jdeta'] = abs(jdeta[0])

            if channel=='mt' or channel=='et':
                cat = int(class_pred[0]+0.1)
                bdt_name = 'bdt_ditau'
                if cat==1:
                    bdt_name = 'bdt_signal'
                if cat==2:
                    bdt_name = 'bdt_fakes'
                variables[bdt_name] = bdt_pred[0]

            variables['aco_lep_pi'] = -99999.
            variables['aco_lep_rho'] = -99999.
            variables['aco_lep_a1_1pr'] = -99999.
            variables['aco_lep_a1_3pr'] = -99999. 
            if isTauToPi_2:
                variables['aco_lep_pi'] = RadToDeg*aco_lep_pi[0]
            if isTauToRho_2:
                variables['aco_lep_rho'] = RadToDeg*aco_lep_rho[0]
            if isTauToA1_1pr_2:
                variables['aco_lep_a1_1pr'] = RadToDeg*aco_lep_rho[0]
            if isTauToA1_3pr_2:
                variables['aco_lep_a1_3pr'] = RadToDeg*aco_lep_a1[0]
            
            ##################
            ## total weight ##
            ##################
            Weight = weight[0]
            # prevent large weights
            if ROOT.TMath.Abs(Weight)>10.0:
                print('weight %3.1f > 10'%(Weight))
                continue

            ############################
            ## applying scale factors ##
            ## for IPSig cuts         ##
            ############################
#            if not cuts.applyIPSigLep1Cut and self.ismc:
#                Weight = Weight/w_IP_Significance[0]
                
            ####################################
            # Applying jet->tau fake factors ###
            ####################################
            FF_Weights = {}
            for ff_label in utils.ff_labels:
                FF_Weights[ff_label] = 1.0
                
            if applyFakeFactor and fakeFactor!=None:
                for FF in FF_Weights:
                    if FF not in ['ar']:
#                        print(pt_2[0],n_jets[0],DecayMode)
                        pt_a = [pt_2[0]]
                        njets_a = [n_jets[0]]
                        dm_a = [DecayMode]
                        FF_Weights[FF] = self.fakeFactor.evaluate(pt_a,
                                                                  njets_a,
                                                                  dm_a,
                                                                  FF,
                                                                  'nom')[0]
#            for FF in FF_Weights:            
#                print('%6s -> %6.4f'%(FF,FF_Weights[FF]))
                    
            #####################################
            ####### JetFakes measurements #######
            #####################################
            if analysisType=='jetFakes':
                if channel=='mm' or channel=='ee': continue
                # histogram for measurement of jet->tau fake
                if region_flags['qcd_ff']:
                    name = 'pt_2_qcd_%s_%s_%s_%s_%s'%(dm_label,njets_label,eta_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_qcd_%s_%s_%s_%s_all'%(dm_label,njets_label,eta_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_qcd_%s_%s_all_%s_%s'%(dm_label,njets_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_qcd_%s_%s_all_%s_all'%(dm_label,njets_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                if region_flags['wj_ff']:
                    name = 'pt_2_wj_%s_%s_%s_%s_%s'%(dm_label,njets_label,eta_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_wj_%s_%s_%s_%s_all'%(dm_label,njets_label,eta_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_wj_%s_%s_all_%s_%s'%(dm_label,njets_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_wj_%s_%s_all_%s_all'%(dm_label,njets_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                if region_flags['lowmt_os_iso']:
                    name = 'pt_2_top_%s_%s_%s_%s_%s'%(dm_label,njets_label,eta_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_top_%s_%s_%s_%s_all'%(dm_label,njets_label,eta_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_top_%s_%s_all_%s_%s'%(dm_label,njets_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_top_%s_%s_all_%s_all'%(dm_label,njets_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                if region_flags['qcd_looseIso_ff']:
                    name = 'pt_2_qcd_looseIso_%s_%s_%s_%s_%s'%(dm_label,njets_label,eta_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_qcd_looseIso_%s_%s_%s_%s_all'%(dm_label,njets_label,eta_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_qcd_looseIso_%s_%s_all_%s_%s'%(dm_label,njets_label,deeptau_label,lep2_label)
                    hists[name].Fill(variables['pt_2'],Weight)
                    name = 'pt_2_qcd_looseIso_%s_%s_all_%s_all'%(dm_label,njets_label,deeptau_label)
                    hists[name].Fill(variables['pt_2'],Weight)

            #####################################
            ######### datacardsPhiCP ############
            #####################################
            if analysisType=='datacardsPhiCP':
                if channel=='ee' or channel=='mm': continue
                
                cat = int(class_pred[0]+0.1)
                bdt_name = 'bdt_ditau'
                if cat==1:
                    bdt_name = 'bdt_signal'
                if cat==2:
                    bdt_name = 'bdt_fakes'
                
                CPWeights = {}
                if applyWeightCP:
                    CPWeights['sm'] = weight_cp_sm[0]
                    CPWeights['ps'] = weight_cp_ps[0]
                    CPWeights['mm'] = weight_cp_mm[0]
                
                if self.ismc:
                    name = '%s_%s_%s_%s_%s'%(bdt_name,sign_label,iso_label,lep_label,lep2_label)
                    hists[name].Fill(bdt_pred[0],Weight)
                    name = '%s_%s_%s_all_%s'%(bdt_name,sign_label,iso_label,lep2_label)
                    hists[name].Fill(bdt_pred[0],Weight)
                    name = '%s_%s_%s_all_all'%(bdt_name,sign_label,iso_label)
                    hists[name].Fill(bdt_pred[0],Weight)
                    if applyWeightCP:
                        for CPWeight in CPWeights:
                            name = '%s_%s_%s_%s_%s_%s'%(bdt_name,sign_label,iso_label,lep_label,lep2_label,CPWeight)
                            hists[name].Fill(bdt_pred[0],Weight*CPWeights[CPWeight])
                            name = '%s_%s_%s_all_%s_%s'%(bdt_name,sign_label,iso_label,lep2_label,CPWeight)
                            hists[name].Fill(bdt_pred[0],Weight*CPWeights[CPWeight])
                            name = '%s_%s_%s_all_all_%s'%(bdt_name,sign_label,iso_label,CPWeight)
                            hists[name].Fill(bdt_pred[0],Weight*CPWeights[CPWeight])
                else:
                    name = '%s_%s_%s_all_all'%(bdt_name,sign_label,iso_label)
                    hists[name].Fill(bdt_pred[0],1.0)
                    
                # the signal category : filling 2D distributions
                if cat==cuts.signalCat:
                    decay_label='phicp_vs_bdt_lep_%s'%(dm_label)
                    phiCP = RadToDeg*aco_lep[dm_label]
                    if self.ismc:
                        name = '%s_%s_%s_%s_%s'%(decay_label,sign_label,iso_label,lep_label,lep2_label)
                        hists[name].Fill(bdt_pred[0],phiCP,Weight)
                        name = '%s_%s_%s_all_%s'%(decay_label,sign_label,iso_label,lep2_label)
                        hists[name].Fill(bdt_pred[0],phiCP,Weight)
                        name = '%s_%s_%s_all_all'%(decay_label,sign_label,iso_label)
                        hists[name].Fill(bdt_pred[0],phiCP,Weight)
                        if applyWeightCP:
                            for CPWeight in CPWeights:
                                name = '%s_%s_%s_%s_%s_%s'%(decay_label,sign_label,iso_label,lep_label,lep2_label,CPWeight)
                                hists[name].Fill(bdt_pred[0],phiCP,Weight*CPWeights[CPWeight])
                                name = '%s_%s_%s_all_%s_%s'%(decay_label,sign_label,iso_label,lep2_label,CPWeight)
                                hists[name].Fill(bdt_pred[0],phiCP,Weight*CPWeights[CPWeight])
                                name = '%s_%s_%s_all_all_%s'%(decay_label,sign_label,iso_label,CPWeight)
                                hists[name].Fill(bdt_pred[0],phiCP,Weight*CPWeights[CPWeight])
                    else:
                        name = '%s_%s_%s_all_all'%(decay_label,sign_label,iso_label)
                        hists[name].Fill(bdt_pred[0],phiCP,1.0)
                        
            ################################
            ## Filling control histograms ##
            ##    baseline analysis       ##
            ################################
            if analysisType=='baseline':
                if channel=='mm' or channel=='ee':
                    for varname in variables:
                        for reg in ['incl_os_iso','incl_ss_iso']:
                            if region_flags[reg]:
                                name = '%s_incl_%s_iso_all'%(varname,sign_label)
                                hists[name].Fill(variables[varname],Weight)
                else:
                    for reg in utils.region_labels:
                        if region_flags[reg]:
                            for varname in variables:
                                if deeptau_label=='nominal':
                                    name = '%s_%s_all'%(varname,reg)
                                    hists[name].Fill(variables[varname],Weight)
                                    name = '%s_%s_%s'%(varname,reg,lep2_label)
                                    hists[name].Fill(variables[varname],Weight)
                                    # number of jets
                                    #name = '%s_%s_%s_all'%(varname,reg,njets_label)
                                    #hists[name].Fill(variables[varname],Weight)
                                    #name = '%s_%s_%s_%s'%(varname,reg,njets_label,lep2_label)
                                    #hists[name].Fill(variables[varname],Weight)
                                else:
                                    for FF in FF_Weights:
                                        name = '%s_%s_%s_all'%(varname,reg,FF)
                                        hists[name].Fill(variables[varname],Weight*FF_Weights[FF])
                                        name = '%s_%s_%s_%s'%(varname,reg,FF,lep2_label)
                                        hists[name].Fill(variables[varname],Weight*FF_Weights[FF])
                                        # number of jets 
                                        #name = '%s_%s_%s_%s_all'%(varname,reg,njets_label,FF)
                                        #hists[name].Fill(variables[varname],Weight*FF_Weight[FF])
                                        #name = '%s_%s_%s_%s_%s'%(varname,reg,njets_label,FF,lep2_label)
                                        #hists[name].Fill(variables[varname],Weight*FF_Weight[FF])

            #######################################
            #### IP significance scale factors ####
            #######################################
            if analysisType=='ipSig':
                # tag-and-probe histograms for mm (ee) channels
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
                    name = 'm_vis_2_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                    name = 'm_vis_2_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
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
                
                # tag-and-probe histograms for mt (et) channels
                if channel=='mt' or channel=='et':
                    # deepTau Vs Jets->
                    if deeptau_label!='nominal': continue

                    # lep1 ->
                    bin_label,binPt,binEta = self.GetPtEtaBinLabels(pt_1[0],eta_1[0])
                    # inclusive selection
                    region_label = 'incl'
                    name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                    name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                    hists[name].Fill(variables['m_vis'],Weight)
                
                    if variables['ipsig_1']>cuts.ipsigLepCut:
                        # print('passed')
                        # passing probes
                        WeightSF = 1.0
                        if lep_label=='lep' and self.ismc and ipSigPromptLepSF!=None:
                            WeightSF *= ipSigPromptLepSF.getSF(pt_1[0],eta_1[0])
                        region_label = 'pass'
                        name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                        hists[name].Fill(variables['m_vis'],Weight)
                        name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                        hists[name].Fill(variables['m_vis'],Weight)
                    else:
                        # failing probes
                        # print('failed')
                        region_label = 'fail'
                        name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                        hists[name].Fill(variables['m_vis'],Weight)
                        name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                        hists[name].Fill(variables['m_vis'],Weight)
                        

            ####################################
            ## phi(CP) studies with signal MC ##
            ####################################
            if analysisType=='phiCP':
                variablesCP = {}
                
                P1 = ROOT.TLorentzVector()
                R1 = ROOT.TLorentzVector()
                P2 = ROOT.TLorentzVector()
                R2 = ROOT.TLorentzVector()

                Pi1_1 = ROOT.TLorentzVector()
                Pi2_1 = ROOT.TLorentzVector()
                Pi3_1 = ROOT.TLorentzVector()
                
                Pi1_2 = ROOT.TLorentzVector()
                Pi2_2 = ROOT.TLorentzVector()
                Pi3_2 = ROOT.TLorentzVector()

                P1_vis = ROOT.TLorentzVector()
                P2_vis = ROOT.TLorentzVector()
                
                P1Gen = ROOT.TLorentzVector()                        
                P1Gen.SetPtEtaPhiM(genPart_pt_1[0],genPart_eta_1[0],genPart_phi_1[0],utils.tau_mass)
                P2Gen = ROOT.TLorentzVector()
                P2Gen.SetPtEtaPhiM(genPart_pt_2[0],genPart_eta_2[0],genPart_phi_2[0],utils.tau_mass)
                
                R1Gen = ROOT.TLorentzVector()
                R2Gen = ROOT.TLorentzVector()
                
                MET = ROOT.TVector3(met_pt[0]*math.cos(met_phi[0]),
                                    met_pt[0]*math.sin(met_phi[0]),
                                    0.0)

                variablesCP['aco_lep_pi'] = -9999.                
                variablesCP['aco_lep_rho'] = -9999.

                variablesCP['aco_lep_a1'] = -9999.
                variablesCP['aco_lep_a1_chi5'] = -9999.
                variablesCP['aco_lep_a1_chi10'] = -9999.
                variablesCP['aco_lep_a1_chi20'] = -9999.
                
                variablesCP['aco_lep_a1_KinFit'] = -9999.
                variablesCP['aco_lep_a1_chi5_KinFit'] = -9999.
                variablesCP['aco_lep_a1_chi10_KinFit'] = -9999.
                variablesCP['aco_lep_a1_chi20_KinFit'] = -9999.
                
                variablesCP['aco_lep_a1_GenPV'] = -9999.
                variablesCP['aco_lep_a1_GenTheta'] = -9999.

                variablesCP['theta_GJ'] = -9999.
                variablesCP['theta_GJ_max'] = -9999.
                
                variablesCP['dAlpha'] = -9999.
                variablesCP['dAlpha_KinFit'] = -9999.
                variablesCP['dAlpha_KinFitSV'] = -9999.
                
                variablesCP['dTheta'] = -9999.
                variablesCP['dTheta_KinFit'] = -9999.
                variablesCP['dTheta_KinFitSV'] = -9999.
                
                variablesCP['dPt'] = -9999.
                variablesCP['dPt_KinFit'] = -9999.
                
                variablesCP['aco_pi_a1'] = -9999.
                variablesCP['aco_rho_a1'] = -9999.
                variablesCP['aco_a1_a1'] = -9999.
                
                variablesCP['aco_pi_a1_chi5'] = -9999.
                variablesCP['aco_rho_a1_chi5'] = -9999.
                variablesCP['aco_a1_a1_chi5'] = -9999.

                variablesCP['aco_pi_a1_chi10'] = -9999.
                variablesCP['aco_rho_a1_chi10'] = -9999.
                variablesCP['aco_a1_a1_chi10'] = -9999.
                
                variablesCP['aco_pi_a1_chi20'] = -9999.
                variablesCP['aco_rho_a1_chi20'] = -9999.
                variablesCP['aco_a1_a1_chi20'] = -9999.
                
                variablesCP['aco_pi_a1_KinFit'] = -9999.
                variablesCP['aco_rho_a1_KinFit'] = -9999.
                variablesCP['aco_a1_a1_KinFit'] = -9999.

                variablesCP['aco_pi_a1_chi5_KinFit'] = -9999.
                variablesCP['aco_rho_a1_chi5_KinFit'] = -9999.
                variablesCP['aco_a1_a1_chi5_KinFit'] = -9999.

                variablesCP['aco_pi_a1_chi10_KinFit'] = -9999.
                variablesCP['aco_rho_a1_chi10_KinFit'] = -9999.
                variablesCP['aco_a1_a1_chi10_KinFit'] = -9999.

                variablesCP['aco_pi_a1_chi20_KinFit'] = -9999.
                variablesCP['aco_rho_a1_chi20_KinFit'] = -9999.
                variablesCP['aco_a1_a1_chi20_KinFit'] = -9999.

                variablesCP['aco_pi_a1_GenPV'] = -9999.
                variablesCP['aco_rho_a1_GenPV'] = -9999.
                variablesCP['aco_a1_a1_GenPV'] = -9999.

                variablesCP['aco_pi_a1_GenTheta'] = -9999.
                variablesCP['aco_rho_a1_GenTheta'] = -9999.
                variablesCP['aco_a1_a1_GenTheta'] = -9999.
                
                variablesCP['aco_pi_pi'] = -9999.
                variablesCP['aco_pi_rho'] = -9999.
                variablesCP['aco_rho_rho'] = -9999.

                variablesCP['mass'] = -9999.
                variablesCP['mass_KinFit'] = -9999.
                
                variablesCP['mass_a1_a1'] = -9999.
                variablesCP['mass_a1_a1_KinFit'] = -9999.

                variablesCP['chi2_pi_a1'] = -99999.
                variablesCP['chi2_rho_a1'] = -99999.
                variablesCP['chi2_a1_a1'] = -99999.
                variablesCP['chi2_lep_a1'] = -99999.
                
                mass1_TT = -100.
                mass2_TT = -100.
                mass12_TT = -99999.
                dummyMass = -100.
                
                firstNegative = charge_1[0]<0

                if channel=='tt':
                    
                    solution1 = False
                    solution2 = False
                    chi2_1 = 10000.
                    chi2_2 = 10000.
                    label1 = 'Impact-Parameter'
                    label2 = 'Impact-Parameter'
                    
                    if isTauToPi_1:
                        label1 = 'Impact-Parameter'
                        P1.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],utils.pi_mass) 
                        R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0.)
                    if isTauToRho_1:
                        label1 = 'Decay-Plane'
                        P1.SetPtEtaPhiM(pi_pt_1[0],pi_eta_1[0],pi_phi_1[0],massPi)
                        R1.SetPtEtaPhiM(pi0_pt_1[0],pi0_eta_1[0],pi0_phi_1[0],massPi0)
                    if isTauToA1_3pr_1:
                        label1 = 'PV'
                        sv = ROOT.TVector3(sv_x_1[0]-PVBS_x[0],sv_y_1[0]-PVBS_y[0],sv_z_1[0]-PVBS_z[0])
                        Pi1_1.SetPtEtaPhiM(pi_pt_1[0],pi_eta_1[0],pi_phi_1[0],massPi)
                        Pi2_1.SetPtEtaPhiM(pi2_pt_1[0],pi2_eta_1[0],pi2_phi_1[0],massPi)
                        Pi3_1.SetPtEtaPhiM(pi3_pt_1[0],pi3_eta_1[0],pi3_phi_1[0],massPi)                        
                        P1_vis = Pi1_1 + Pi2_1 + Pi3_1
                        P2_vis.SetPtEtaPhiM(pt_2[0],eta_2[0],phi_2[0],mass_2[0])
                        P_os,P_ss1,P_ss2 = pv_utils.sortA1(Pi1_1,Pi2_1,Pi3_1,pi_charge_1[0],pi2_charge_1[0],pi3_charge_1[0])
                        if cppKinFit:
                            #Ptau,Pother,chi2 = pv_utils.FastMTT_3pr_cpp(P1_vis,P2_vis,
                            #                                            met_pt,met_phi,met_covXX,met_covXY,met_covYY,
                            #                                            sv,
                            #                                            sv_cov00_1,sv_cov10_1,sv_cov11_1,
                            #                                            sv_cov20_1,sv_cov21_1,sv_cov22_1,
                            #                                            dummyMass)
                            #if chi2<1.0e+10: mass1_TT = (Ptau+Pother).M()
                            Ptau,Pother,chi2 = pv_utils.FastMTT_3pr_cpp(P1_vis,P2_vis,
                                                                        met_pt,met_phi,
                                                                        met_covXX,met_covXY,met_covYY,
                                                                        sv,
                                                                        sv_cov00_1,sv_cov10_1,sv_cov11_1,
                                                                        sv_cov20_1,sv_cov21_1,sv_cov22_1,
                                                                        phiScan,
                                                                        utils.massH)
                        else:
                            Pother,Ptau,chi2 = pv_utils.FastMTT_3pr(P2_vis,P1_vis,
                                                                    MET,
                                                                    met_covXX,met_covXY,met_covYY,
                                                                    sv,
                                                                    sv_cov00_1,sv_cov10_1,sv_cov11_1,
                                                                    sv_cov20_1,sv_cov21_1,sv_cov22_1)
                        
                        if chi2<1.0e+10:
                            solution1 = True
                            chi2_1 = chi2
                            results = pv_utils.AnalysisGenerator(pt_1_FastMTT[0],Ptau,P1_vis,P1Gen,sv)
                            variablesCP['dPt'] = results['dPt']
                            variablesCP['dPt_KinFit'] = results['dPt_KinFit']
                            variablesCP['dAlpha'] = results['dAlpha']
                            variablesCP['dAlpha_KinFit'] = results['dAlpha_KinFit']
                            variablesCP['dAlpha_KinFitSV'] = results['dAlpha_KinFitSV']
                            variablesCP['dTheta'] = results['dTheta']
                            variablesCP['dTheta_KinFit'] = results['dTheta_KinFit']
                            variablesCP['dTheta_KinFitSV'] = results['dTheta_KinFitSV']
                            
                        P1 = Ptau
                        R1 = pv_utils.PolVectA1(Ptau,P_os,P_ss1,P_ss2,charge_1[0])
                        R1Gen = pv_utils.PolVectA1(P1Gen,P_os,P_ss1,P_ss2,charge_1[0])
                        P1ThetaGen = pv_utils.GeneratorVector(P1Gen,P1_vis,sv)                        
                        R1ThetaGen = pv_utils.PolVectA1(P1ThetaGen,P_os,P_ss1,P_ss2,charge_1[0])
                        
                    if isTauToPi_2:
                        label2 = 'Impact-Parameter'
                        P2.SetPtEtaPhiM(pt_2[0],eta_2[0],phi_2[0],utils.pi_mass) 
                        R2.SetXYZT(ip_x_2[0],ip_y_2[0],ip_z_2[0],0.)
                    if isTauToRho_2:
                        label2 = 'Decay-Plane'
                        P2.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                        R2.SetPtEtaPhiM(pi0_pt_2[0],pi0_eta_2[0],pi0_phi_2[0],massPi0)
                    if isTauToA1_3pr_2:
                        label2 = 'PV'
                        sv = ROOT.TVector3(sv_x_2[0]-PVBS_x[0],sv_y_2[0]-PVBS_y[0],sv_z_2[0]-PVBS_z[0])
                        
                        Pi1_2.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                        Pi2_2.SetPtEtaPhiM(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],massPi)
                        Pi3_2.SetPtEtaPhiM(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],massPi)                        
                        P2_vis = Pi1_2 + Pi2_2 + Pi3_2
                        P1_vis.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],mass_1[0])
                        P_os,P_ss1,P_ss2 = pv_utils.sortA1(Pi1_2,Pi2_2,Pi3_2,pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])
                        if cppKinFit:
                            #Ptau,Pother,chi2 = pv_utils.FastMTT_3pr_cpp(P2_vis,P1_vis,
                            #                                            met_pt,met_phi,
                            #                                            met_covXX,met_covXY,met_covYY,
                            #                                            sv,
                            #                                            sv_cov00_2,sv_cov10_2,sv_cov11_2,
                            #                                            sv_cov20_2,sv_cov21_2,sv_cov22_2,
                            #                                            dummyMass)
                            #if chi2<1.0e+10: mass2_TT = (Ptau+Pother).M()                        
                            Ptau,Pother,chi2 = pv_utils.FastMTT_3pr_cpp(P2_vis,P1_vis,
                                                                        met_pt,met_phi,
                                                                        met_covXX,met_covXY,met_covYY,
                                                                        sv,
                                                                        sv_cov00_2,sv_cov10_2,sv_cov11_2,
                                                                        sv_cov20_2,sv_cov21_2,sv_cov22_2,
                                                                        phiScan,
                                                                        utils.massH)
                        else:
                            Pother,Ptau,chi2 = pv_utils.FastMTT_3pr(P1_vis,P2_vis,
                                                                    MET,
                                                                    met_covXX,met_covXY,met_covYY,
                                                                    sv,
                                                                    sv_cov00_2,sv_cov10_2,sv_cov11_2,
                                                                    sv_cov20_2,sv_cov21_2,sv_cov22_2)
                        
                        if chi2<1.0e+10:
                            solution2 = True
                            chi2_2 = chi2
                            results = pv_utils.AnalysisGenerator(pt_2_FastMTT[0],Ptau,P2_vis,P2Gen,sv)
                            variablesCP['dPt'] = results['dPt']
                            variablesCP['dPt_KinFit'] = results['dPt_KinFit']
                            variablesCP['dAlpha'] = results['dAlpha']
                            variablesCP['dAlpha_KinFit'] = results['dAlpha_KinFit']
                            variablesCP['dAlpha_KinFitSV'] = results['dAlpha_KinFitSV']
                            variablesCP['dTheta'] = results['dTheta']
                            variablesCP['dTheta_KinFit'] = results['dTheta_KinFit']
                            variablesCP['dTheta_KinFitSV'] = results['dTheta_KinFitSV']
                            
                        P2 = Ptau
                        R2 = pv_utils.PolVectA1(Ptau,P_os,P_ss1,P_ss2,charge_2[0])
                        R2Gen = pv_utils.PolVectA1(P2Gen,P_os,P_ss1,P_ss2,charge_2[0])
                        P2ThetaGen = pv_utils.GeneratorVector(P2Gen,P2_vis,sv)                        
                        R2ThetaGen = pv_utils.PolVectA1(P2ThetaGen,P_os,P_ss1,P_ss2,charge_1[0])
                    
                    # computing acoplanarities
                    if isTauToA1_3pr_1 and (isTauToPi_2 or isTauToRho_2):
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNegative,label1,label2)
                        aco_gen = pv_utils.acoCP(P1Gen,P2,R1Gen,R2,firstNegative,label1,label2)
                        aco_genTheta = pv_utils.acoCP(P1ThetaGen,P2,R1ThetaGen,R2,firstNegative,label1,label2)
                        if isTauToPi_2:
                            variablesCP['aco_pi_a1_GenPV'] = RadToDeg*aco_gen
                            variablesCP['aco_pi_a1_GenTheta'] = RadToDeg*aco_genTheta
                            variablesCP['mass'] = FastMTT_mass[0]
                            variablesCP['mass_KinFit'] = mass1_TT
                            if solution1:
                                variablesCP['chi2_pi_a1'] = chi2_1
                                variablesCP['aco_pi_a1_KinFit'] = RadToDeg*aco
                                if chi2_1<5:
                                    variablesCP['aco_pi_a1_chi5_KinFit'] = RadToDeg*aco
                                if chi2_1<10:
                                    variablesCP['aco_pi_a1_chi10_KinFit'] = RadToDeg*aco
                                if chi2_1<20:
                                    variablesCP['aco_pi_a1_chi20_KinFit'] = RadToDeg*aco
                        if isTauToRho_2:
                            variablesCP['aco_rho_a1_GenPV'] = RadToDeg*aco_gen
                            variablesCP['aco_rho_a1_GenTheta'] = RadToDeg*aco_genTheta
                            variablesCP['mass'] = FastMTT_mass[0]
                            variablesCP['mass_KinFit'] = mass1_TT
                            if solution1:
                                variablesCP['chi2_rho_a1'] = chi2_1
                                variablesCP['aco_rho_a1_KinFit'] = RadToDeg*aco
                                if chi2_1<5:
                                    variablesCP['aco_rho_a1_chi5_KinFit'] = RadToDeg*aco
                                if chi2_1<10:
                                    variablesCP['aco_rho_a1_chi10_KinFit'] = RadToDeg*aco
                                if chi2_1<20:
                                    variablesCP['aco_rho_a1_chi20_KinFit'] = RadToDeg*aco
                            
                    if isTauToA1_3pr_2 and (isTauToPi_1 or isTauToRho_1):
                        aco = pv_utils.acoCP(P1,P2,R1,R2,firstNegative,label1,label2)
                        aco_gen = pv_utils.acoCP(P1,P2Gen,R1,R2Gen,firstNegative,label1,label2)
                        aco_genTheta = pv_utils.acoCP(P1,P2ThetaGen,R1,R2ThetaGen,firstNegative,label1,label2)
                        if isTauToPi_1:
                            variablesCP['aco_pi_a1_GenPV'] = RadToDeg*aco_gen
                            variablesCP['aco_pi_a1_GenTheta'] = RadToDeg*aco_genTheta
                            variablesCP['mass'] = FastMTT_mass[0]
                            variablesCP['mass_KinFit'] = mass2_TT
                            if solution2:
                                variablesCP['chi2_pi_a1'] = chi2_2
                                variablesCP['aco_pi_a1_KinFit'] = RadToDeg*aco
                                if chi2_2<5:
                                    variablesCP['aco_pi_a1_chi5_KinFit'] = RadToDeg*aco
                                if chi2_2<10:
                                    variablesCP['aco_pi_a1_chi10_KinFit'] = RadToDeg*aco
                                if chi2_2<20:
                                    variablesCP['aco_pi_a1_chi20_KinFit'] = RadToDeg*aco
                        if isTauToRho_1:
                            variablesCP['aco_rho_a1_GenPV'] = RadToDeg*aco_gen
                            variablesCP['aco_rho_a1_GenTheta'] = RadToDeg*aco_genTheta
                            variablesCP['mass'] = FastMTT_mass[0]
                            variablesCP['mass_KinFit'] = mass2_TT
                            if solution2:
                                variablesCP['chi2_rho_a1'] = chi2_2
                                variablesCP['aco_rho_a1_KinFit'] = RadToDeg*aco
                                if chi2_2<5:
                                    variablesCP['aco_rho_a1_chi5_KinFit'] = RadToDeg*aco
                                if chi2_2<10:
                                    variablesCP['aco_rho_a1_chi10_KinFit'] = RadToDeg*aco
                                if chi2_2<20:
                                    variablesCP['aco_rho_a1_chi20_KinFit'] = RadToDeg*aco
                                
                    if isTauToA1_3pr_1 and isTauToA1_3pr_2:
                        
                        sv1 = ROOT.TVector3(sv_x_1[0]-PVBS_x[0],sv_y_1[0]-PVBS_y[0],sv_z_1[0]-PVBS_z[0])
                        Pi1_1.SetPtEtaPhiM(pi_pt_1[0],pi_eta_1[0],pi_phi_1[0],massPi)
                        Pi2_1.SetPtEtaPhiM(pi2_pt_1[0],pi2_eta_1[0],pi2_phi_1[0],massPi)
                        Pi3_1.SetPtEtaPhiM(pi3_pt_1[0],pi3_eta_1[0],pi3_phi_1[0],massPi)                        
                        P1_vis = Pi1_1 + Pi2_1 + Pi3_1
                        
                        sv2 = ROOT.TVector3(sv_x_2[0]-PVBS_x[0],sv_y_2[0]-PVBS_y[0],sv_z_2[0]-PVBS_z[0])
                        Pi1_2.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                        Pi2_2.SetPtEtaPhiM(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],massPi)
                        Pi3_2.SetPtEtaPhiM(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],massPi)                        
                        P2_vis = Pi1_2 + Pi2_2 + Pi3_2

                        P1_os,P1_ss1,P1_ss2 = pv_utils.sortA1(Pi1_1,Pi2_1,Pi3_1,pi_charge_1[0],pi2_charge_1[0],pi3_charge_1[0])
                        P2_os,P2_ss1,P2_ss2 = pv_utils.sortA1(Pi1_2,Pi2_2,Pi3_2,pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])

                        if cppKinFit:
                            #P1_tau,P2_tau,chi2 = pv_utils.FastMTT_3pr_3pr_cpp(P1_vis,P2_vis,
                            #                                                  met_pt, met_phi,
                            #                                                  met_covXX,met_covXY,met_covYY,
                            #                                                  sv1,
                            #                                                  sv_cov00_1,sv_cov10_1,sv_cov11_1,
                            #                                                  sv_cov20_1,sv_cov21_1,sv_cov22_1,
                            #                                                  sv2,
                            #                                                  sv_cov00_2,sv_cov10_2,sv_cov11_2,
                            #                                                  sv_cov20_2,sv_cov21_2,sv_cov22_2,
                            #                                                  dummyMass)
                            #if chi2<1.0e+10:
                            #    variablesCP['mass_a1_a1_KinFit'] = (P1_tau+P2_tau).M()
                            #    variablesCP['mass_a1_a1'] = FastMTT_mass[0]
                            P1_tau,P2_tau,chi2 = pv_utils.FastMTT_3pr_3pr_cpp(P1_vis,P2_vis,
                                                                              met_pt, met_phi,
                                                                              met_covXX,met_covXY,met_covYY,
                                                                              sv1,
                                                                              sv_cov00_1,sv_cov10_1,sv_cov11_1,
                                                                              sv_cov20_1,sv_cov21_1,sv_cov22_1,
                                                                              sv2,
                                                                              sv_cov00_2,sv_cov10_2,sv_cov11_2,
                                                                              sv_cov20_2,sv_cov21_2,sv_cov22_2,
                                                                              utils.massH)
                        else:
                            P1_tau,P2_tau,chi2 = pv_utils.FastMTT_3pr_3pr(P1_vis,P2_vis,
                                                                          MET,
                                                                          met_covXX,met_covXY,met_covYY,
                                                                          sv1,
                                                                          sv_cov00_1,sv_cov10_1,sv_cov11_1,
                                                                          sv_cov20_1,sv_cov21_1,sv_cov22_1,
                                                                          sv2,
                                                                          sv_cov00_2,sv_cov10_2,sv_cov11_2,
                                                                          sv_cov20_2,sv_cov21_2,sv_cov22_2,
                                                                          )      
                        
                        if chi2<1.0e+10:
                            R1_tau = pv_utils.PolVectA1(P1_tau,P1_os,P1_ss1,P1_ss2,charge_1[0])
                            R2_tau = pv_utils.PolVectA1(P2_tau,P2_os,P2_ss1,P2_ss2,charge_2[0])
                            aco_a1_a1_comb = pv_utils.acoCP(P1_tau,P2_tau,R1_tau,R2_tau,firstNegative,'PV','PV')
                            aco_a1_a1_comb += math.pi
                            if aco_a1_a1_comb > 2.0*math.pi: aco_a1_a1_comb -= 2.0*math.pi
                            variablesCP['aco_a1_a1_KinFit'] = RadToDeg*aco_a1_a1_comb
                            variablesCP['chi2_a1_a1'] = chi2
                            if chi2<20:
                                variablesCP['aco_a1_a1_chi20_KinFit'] = RadToDeg*aco_a1_a1_comb
                            if chi2<10:
                                variablesCP['aco_a1_a1_chi10_KinFit'] = RadToDeg*aco_a1_a1_comb
                            if chi2<5:
                                variablesCP['aco_a1_a1_chi5_KinFit'] = RadToDeg*aco_a1_a1_comb

                        aco_a1_a1_gen = pv_utils.acoCP(P1Gen,P2Gen,R1Gen,R2Gen,firstNegative,'PV','PV')
                        aco_a1_a1_gen += math.pi
                        if aco_a1_a1_gen>2.0*math.pi: aco_a1_a1_gen -= 2.0*math.pi
                        
                        P1ThetaGen = pv_utils.GeneratorVector(P1Gen,P1_vis,sv1)
                        R1ThetaGen = pv_utils.PolVectA1(P1ThetaGen,P1_os,P1_ss1,P1_ss2,charge_1[0])
                        P2ThetaGen = pv_utils.GeneratorVector(P2Gen,P2_vis,sv2)
                        R2ThetaGen = pv_utils.PolVectA1(P2ThetaGen,P2_os,P2_ss1,P2_ss2,charge_1[0])
                        aco_a1_a1_gen_Theta = pv_utils.acoCP(P1ThetaGen,P2ThetaGen,R1ThetaGen,R2ThetaGen,firstNegative,'PV','PV')
                        aco_a1_a1_gen_Theta += math.pi
                        if aco_a1_a1_gen_Theta>2.0*math.pi: aco_a1_a1_gen_Theta -= 2.0*math.pi                        
                        variablesCP['aco_a1_a1_GenPV'] = RadToDeg*aco_a1_a1_gen
                        variablesCP['aco_a1_a1_GenTheta'] = RadToDeg*aco_a1_a1_gen_Theta
                        
                    # acoplanarities from tuple
                    if isTauToPi_1 and isTauToPi_2:
                        variablesCP['aco_pi_pi'] = RadToDeg*aco_pi_pi[0]
                    if isTauToRho_1 and isTauToRho_2:
                        variablesCP['aco_rho_rho'] = RadToDeg*aco_rho_rho[0]
                    if isTauToA1_3pr_1 and isTauToA1_3pr_2:
                        variablesCP['aco_a1_a1'] = RadToDeg*aco_a1_a1[0]
                        
                    if isTauToPi_1 and isTauToRho_2:
                        variablesCP['aco_pi_rho'] = RadToDeg*aco_pi_rho[0]
                    if isTauToRho_1 and isTauToPi_2:
                        variablesCP['aco_pi_rho'] = RadToDeg*aco_rho_pi[0]
                        
                    if isTauToPi_1 and isTauToA1_3pr_2:
                        variablesCP['aco_pi_a1'] = RadToDeg*aco_pi_a1[0]                        
                    if isTauToA1_3pr_1 and isTauToPi_2:
                        variablesCP['aco_pi_a1'] = RadToDeg*aco_a1_pi[0]

                    if isTauToRho_1 and isTauToA1_3pr_2:
                        variablesCP['aco_rho_a1'] = RadToDeg*aco_rho_a1[0]                        
                    if isTauToA1_3pr_1 and isTauToRho_2:
                        variablesCP['aco_rho_a1'] = RadToDeg*aco_a1_rho[0]

                ####################################
                ###### Leptonic channels ###########
                ####################################
                if channel=='mt' or channel=='et':
                    accept = sign_label=='os' and iso_label=='iso' and deeptau_label=='nominal' and mt_label=='lowmt'
                    if not accept: continue 

                    solution = False
                    leptonDecay = True
                    
                    firstNegative = charge_1[0] < 0
                    P1 = ROOT.TLorentzVector()
                    P1.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],massLep)
                    P1_vis = ROOT.TLorentzVector()
                    P1_vis.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],massLep)
                    R1 = ROOT.TLorentzVector()
                    R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0.0)
                    if isTauToPi_2:
                        # tau->pi+v decay 
                        variablesCP['aco_lep_pi'] = RadToDeg*aco_lep_pi[0]
                    if isTauToRho_2:
                        # tau->rho(p+pi-)+v decay
                        variablesCP['aco_lep_rho'] = RadToDeg * aco_lep_rho[0]                        
                    if isTauToA1_3pr_2:
                        # tau->a1(3pr)+v decay                        
                        sv = ROOT.TVector3(sv_x_2[0]-PVBS_x[0],sv_y_2[0]-PVBS_y[0],sv_z_2[0]-PVBS_z[0])
                        Pi1_2.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                        Pi2_2.SetPtEtaPhiM(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],massPi)
                        Pi3_2.SetPtEtaPhiM(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],massPi)                        
                        P2_vis = Pi1_2 + Pi2_2 + Pi3_2
                        variablesCP['aco_lep_a1'] = RadToDeg*aco_lep_a1[0]
                        P1Gen.SetPtEtaPhiM(genPart_pt_1[0],genPart_eta_1[0],genPart_phi_1[0],utils.tau_mass)
                        P2Gen.SetPtEtaPhiM(genPart_pt_2[0],genPart_eta_2[0],genPart_phi_2[0],utils.tau_mass)
                        P_os,P_ss1,P_ss2 = pv_utils.sortA1(Pi1_2,Pi2_2,Pi3_2,pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])
                        if cppKinFit:
                            Ptau,Pother,chi2 = pv_utils.FastMTT_3pr_cpp(P2_vis,P1_vis,
                                                                        met_pt,met_phi,
                                                                        met_covXX,met_covXY,met_covYY,
	                                                                sv,
	                                                                sv_cov00_2,sv_cov10_2,sv_cov11_2,
                                                                        sv_cov20_2,sv_cov21_2,sv_cov22_2,
                                                                        phiScan,
                                                                        utils.massH)
                        else:
                            Pother, Ptau, chi2 = pv_utils.FastMTT_3pr(P1_vis,P2_vis,
                                                                      MET,
                                                                      met_covXX,met_covXY,met_covYY,
                                                                      sv,
                                                                      sv_cov00_2,sv_cov10_2,sv_cov11_2,
                                                                      sv_cov20_2,sv_cov21_2,sv_cov22_2)
                        
                        if chi2<1.0e+10:
                            results = pv_utils.AnalysisGenerator(pt_2_FastMTT[0],Ptau,P2_vis,P2Gen,sv)
                            variablesCP['dPt'] = results['dPt']
                            variablesCP['dPt_KinFit'] = results['dPt_KinFit']
                            variablesCP['dAlpha'] = results['dAlpha']
                            variablesCP['dAlpha_KinFit'] = results['dAlpha_KinFit']
                            variablesCP['dTheta'] = results['dTheta']
                            variablesCP['dTheta_KinFit'] = results['dTheta_KinFit']
                            variablesCP['dTheta_KinFitSV'] = results['dTheta_KinFitSV']
                            R2 = pv_utils.PolVectA1(Ptau,P_os,P_ss1,P_ss2,charge_2[0])
                            P2 = Ptau
                            aco = pv_utils.acoCP(P1,P2,R1,R2,firstNegative,'Impact-Parameter','PV')
                            variablesCP['aco_lep_a1_KinFit'] = RadToDeg*aco
                            variablesCP['chi2_lep_a1'] = chi2
                            if chi2<5:
                                variablesCP['aco_lep_a1_chi5_KinFit'] = RadToDeg*aco
                            if chi2<10:
                                variablesCP['aco_lep_a1_chi10_KinFit'] = RadToDeg*aco
                            if chi2<20:
                                variablesCP['aco_lep_a1_chi20_KinFit'] = RadToDeg*aco
                        R2Gen = pv_utils.PolVectA1(P2Gen,P_os,P_ss1,P_ss2,charge_2[0])                        
                        aco_Gen = pv_utils.acoCP(P1,P2Gen,R1,R2Gen,firstNegative,'Impact-Parameter','PV')
                        
                        P2_GenTheta = pv_utils.GeneratorVector(P2Gen,P2_vis,sv)
                        R2_GenTheta = pv_utils.PolVectA1(P2_GenTheta,P_os,P_ss1,P_ss2,charge_2[0])
                        aco_GenTheta = pv_utils.acoCP(P1,P2_GenTheta,R1,R2_GenTheta,firstNegative,'Impact-Parameter','PV')
                        variablesCP['aco_lep_a1_GenPV'] = RadToDeg*aco_Gen
                        variablesCP['aco_lep_a1_GenTheta'] = RadToDeg*aco_GenTheta

                for varname in variablesCP:
                    cp_weigths = {
                        'sm' : weight_cp_sm[0],
                        'ps' : weight_cp_ps[0],
                        'mm' : weight_cp_mm[0],
                    }
                    hists[varname].Fill(variablesCP[varname],Weight)
                    if applyWeightCP:
                        for cp_weight in cp_weigths: 
                            name = '%s_%s'%(varname,cp_weight)
                            hists[name].Fill(variablesCP[varname],Weight*cp_weigths[cp_weight])


            
        for hist in hists:
            hists[hist].Scale(self.norm)
        return hists



            
