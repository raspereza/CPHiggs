import ROOT 
import math
from array import array
import numpy as np
import os
import CPHiggs.IP.utils as utils

# Cuts for Z->tau+tau and Z->ll selection
class AnalysisCuts:
    def __init__(self,**kwargs):
        self.mtCut = kwargs.get('mtCut',999999.)
        
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

        self.isoLepCut = kwargs.get('isoLepCut',0.15)

        self.ipsigLepCut = kwargs.get('ipsigLepCut',1.0)

        self.isoLepInverseLowerCut = kwargs.get('isoLepInverseLowerCut',0.20)
        self.isoLepInverseUpperCut = kwargs.get('isoLepInverseUpperCut',0.50)
        
        print('')
        print("Setting cuts ->")

        print("mtCut",self.mtCut)

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
        self.norm = norm
        self.isdata = isdata
        if isdata:
            self.norm = 1.0

        self.sign_labels = ['os','ss']
        self.iso_labels = ['iso','antiiso','rest']
        self.type_labels = ['tau','lep','had','all']
            
        print('%s : %s : %s : norm = %7.3f'%(era,self.channel,self.sampleName,self.norm))
        
    def getSampleName(self):
        return self.sampleName
        
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

    def SetConfig(self,analysisCuts,histPtBins,histEtaBins):
        self.analysisCuts = analysisCuts

        self.histEtaBins = histEtaBins
        self.histPtBins = histPtBins

        self.nbinsPt = histPtBins.GetNbinsX()
        self.nbinsEta = histEtaBins.GetNbinsX()

        self.ptMin = histPtBins.GetBinLowEdge(1)
        self.ptMax = histPtBins.GetBinLowEdge(self.nbinsPt+1)

        self.etaMin = histEtaBins.GetBinLowEdge(1)
        self.etaMax = histEtaBins.GetBinLowEdge(self.nbinsEta+1)

        
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

        print(pt,eta,binPt,binEta,binLabel)
        
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
                        nbins = utils.lib_histos[var][0]
                        xmin = utils.lib_histos[var][1]
                        xmax = utils.lib_histos[var][1]
                        hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))


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
                            for region in ['pass','fail']:
                                name = 'm_vis_%s_%s_%s_%s_%s'%(region,label,sign,iso,typ)
                                histname = self.sampleName+'_'+ name
                                hists[name] = ROOT.TH1D(histname,"",nbins,array('d',list(xbins)))

#        for hist in hists:
#            print(hist)
        return hists

    def CreateHistosTuple(self):

        tree = self.sampleFile.Get('ntuple')

        # initialization
        cuts = self.analysisCuts
        channel = self.channel
        sampleName = self.sampleName
        
        # creating histograms 
        hists = self.DeclareHistos()

        # floats
        weight      = np.zeros(1,dtype=np.float64)

        pt_1        = np.zeros(1,dtype=np.float64)
        eta_1       = np.zeros(1,dtype=np.float64)
        phi_1       = np.zeros(1,dtype=np.float64)
        iso_1       = np.zeros(1,dtype=np.float64)        
        mt_1        = np.zeros(1,dtype=np.float64)
        ipsig_1     = np.zeros(1,dtype=np.float64)
        
        pt_2        = np.zeros(1,dtype=np.float64)
        eta_2       = np.zeros(1,dtype=np.float64)
        phi_2       = np.zeros(1,dtype=np.float64)
        iso_2       = np.zeros(1,dtype=np.float64)
        mt_2        = np.zeros(1,dtype=np.float64)
        ipsig_2     = np.zeros(1,dtype=np.float64)
        
        m_vis        = np.zeros(1,dtype=np.float64)
        met_pt       = np.zeros(1,dtype=np.float64)
        
        # booleans
        trg_lep     = np.zeros(1,dtype='?')
        trg_cross   = np.zeros(1,dtype='?')
        os          = np.zeros(1,dtype='?')
        
        # integers
        idDeepTau2018v2p5VSe_2   = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSmu_2  = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSjet_2 = np.zeros(1,dtype=np.int64)

        genmatch_1               = np.zeros(1,dtype=np.int64)
        genmatch_2               = np.zeros(1,dtype=np.int64)

        # branches ->

        # floats ->
        tree.SetBranchAddress('met_pt',met_pt)
        tree.SetBranchAddress('m_vis',m_vis)
        tree.SetBranchAddress('os',os)
        
        tree.SetBranchAddress('mt_1',mt_1)
        tree.SetBranchAddress('pt_1',pt_1)
        tree.SetBranchAddress('eta_1',eta_1)
        tree.SetBranchAddress('phi_1',phi_1)
        tree.SetBranchAddress('iso_1',iso_1)
        tree.SetBranchAddress('ip_LengthSig_1',ipsig_1)
        
        tree.SetBranchAddress('pt_2',pt_2)
        tree.SetBranchAddress('eta_2',eta_2)
        tree.SetBranchAddress('phi_2',phi_2)
        tree.SetBranchAddress('mt_2',mt_2)
        tree.SetBranchAddress('ip_LengthSig_2',ipsig_2)
              
        if channel=='mm' or channel=='ee':
            tree.SetBranchAddress('iso_2',iso_2)
        
        tree.SetBranchAddress('weight',weight)
        
        # booleans (trigger)
        if channel=='mt' or channel=='mm': tree.SetBranchAddress('trg_singlemuon',trg_lep)
        if channel=='mt': tree.SetBranchAddress('trg_mt_cross',trg_cross)
        if channel=='et' or channel=='ee': tree.SetBranchAddress('trg_singleelectron',trg_lep)
        if channel=='et': tree.SetBranchAddress('trg_et_cross',trg_cross)
            
        # integers
        if channel=='mt' or channel=='et':
            tree.SetBranchAddress('idDeepTau2018v2p5VSe_2',idDeepTau2018v2p5VSe_2)
            tree.SetBranchAddress('idDeepTau2018v2p5VSmu_2',idDeepTau2018v2p5VSmu_2)
            tree.SetBranchAddress('idDeepTau2018v2p5VSjet_2',idDeepTau2018v2p5VSjet_2)

        if not self.isdata:
            tree.SetBranchAddress("genPartFlav_1",genmatch_1)
            tree.SetBranchAddress("genPartFlav_2",genmatch_2)
            
        nentries = tree.GetEntries()

        # run over entries
        for entry in range(0,nentries):

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


            if not passTrigger: continue
            
            # kinematic cuts
            if pt_1[0]<cuts.ptLep1Cut: continue
            if abs(eta_1[0])>cuts.etaLep1Cut: continue
            if iso_1[0]>cuts.isoLepInverseUpperCut: continue
            
            if pt_2[0]<cuts.ptLep2Cut: continue
            if abs(eta_2[0])>cuts.etaLep2Cut: continue

            if channel=='mm' or channel=='ee':
                if iso_2[0]>cuts.isoLepCut: continue
            
            if channel=='mt' or channel=='et':
                # mT cut
                if mt_1[0]>cuts.mtCut: continue
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
            
            lep_label = 'had'
            sign_label = 'os'
            iso_label = 'rest'
            
            if not self.isdata:
                if genmatch_1[0]==1: lep_label = 'lep'
                if genmatch_1[0]==15: lep_label = 'tau'

            if not os[0]:
                sign_label = 'ss'
                
            directIso = iso_1[0] < cuts.isoLepCut
            inverseIso = iso_1[0] > cuts.isoLepInverseLowerCut and iso_1[0] < cuts.isoLepInverseUpperCut

            if directIso: iso_label = 'iso'
            if inverseIso: iso_label = 'antiiso'

            Weight = weight[0]
            
            for varname in variables:
                name = '%s_%s_%s_all'%(varname,sign_label,iso_label)
                hists[name].Fill(variables[varname],Weight)
                name = '%s_%s_%s_%s'%(varname,sign_label,iso_label,lep_label)
                hists[name].Fill(variables[varname],Weight)

           
            if channel=='mm' or channel=='ee':
                # lep1
                bin_label,binPt,binEta = self.GetPtEtaBinLabels(pt_1[0],eta_1[0])
                region_label = 'fail'
                if variables['ipsig_1']>cuts.ipsigLepCut:
                    region_label = 'pass'
                name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                hists[name].Fill(variables['m_vis'],Weight)
                name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                hists[name].Fill(variables['m_vis'],Weight)
                # lep2
                bin_label,binPt,binEta = self.GetPtEtaBinLabels(pt_2[0],eta_2[0])
                region_label = 'fail'
                if variables['ipsig_2']>cuts.ipsigLepCut:
                    region_label = 'pass'
                name = 'm_vis_%s_%s_%s_%s_%s'%(region_label,bin_label,sign_label,iso_label,lep_label)
                hists[name].Fill(variables['m_vis'],Weight)
                name = 'm_vis_%s_%s_%s_%s_all'%(region_label,bin_label,sign_label,iso_label)
                hists[name].Fill(variables['m_vis'],Weight)
                
                
        for hist in hists:
            hists[hist].Scale(self.norm)
        return hists
