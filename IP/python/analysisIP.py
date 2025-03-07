import ROOT 
import math
from array import array
import numpy as np
import os
import CPHiggs.IP.utils as utils
from CPHiggs.IP.ScaleFactor import ScaleFactor
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1

# three LorentzVectors (pions) and three doubles (charges) 
# sorts three charged pions
def sortA1(pi1_input,pi2_input,pi3_input,q1,q2,q3):

    pi1 = ROOT.TLorentzVector()
    pi1.SetXYZT(pi1_input.X(),pi1_input.Y(),pi1_input.Z(),pi1_input.T())
    pi2 = ROOT.TLorentzVector()
    pi2.SetXYZT(pi2_input.X(),pi2_input.Y(),pi2_input.Z(),pi2_input.T())    
    pi3 = ROOT.TLorentzVector()
    pi3.SetXYZT(pi3_input.X(),pi3_input.Y(),pi3_input.Z(),pi3_input.T())

    condition1 = (q1*q2)<0. and (q1*q3)<0.
    condition2 = (q2*q1)<0. and (q2*q3)<0.

    if condition1:
        p1 = pi1
        p2 = pi2
        p3 = pi3
    elif condition2:
        p1 = pi2
        p2 = pi1
        p3 = pi3
    else:
        p1 = pi3
        p2 = pi1
        p3 = pi2


    dM12 = abs((p1+p2).M()-utils.rho_mass)
    dM13 = abs((p1+p3).M()-utils.rho_mass)

    if dM13<dM12:
        temp = p2
        p2 = p3
        p3 = temp

    """
    print('a1 sorted -> ')
    print('pi_os  (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p1.Pt(),p1.Eta(),p1.Phi(),p1.M()))
    print('pi_ss1 (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p2.Pt(),p2.Eta(),p2.Phi(),p2.M()))
    print('pi_ss2 (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p3.Pt(),p3.Eta(),p3.Phi(),p3.M()))
    print('M(os_ss1) = %5.4f -- M(os_ss2) = %5.4f -- m(rho) = %5.4f'%((p1+p2).M(),(p1+p3).M(),utils.rho_mass))
    """
    
    return p1,p2,p3


def rotateToGJ(visTau,tau):

    vistau_P = visTau.Vect()
    tau_P = tau.Vect()
    mass_vis_tau = visTau.M()
    
    tau_mom = tau.P()
    vistau_mom = visTau.P()

    tau_dir = vistau_P.Unit()

    cos_theta_GJX_1 = -0.5*(utils.tau_mass*utils.tau_mass+mass_vis_tau*mass_vis_tau+2.0*tau.E()*visTau.E())/(tau_mom*vistau_mom)
    cos_theta_GJX_2 = -0.5*(utils.tau_mass*utils.tau_mass+mass_vis_tau*mass_vis_tau-2.0*tau.E()*visTau.E())/(tau_mom*vistau_mom)
    cos_theta_GJ_1 = max(-1.0,min(1.0,cos_theta_GJX_1))
    cos_theta_GJ_2 = max(-1.0,min(1.0,cos_theta_GJX_2))

    sin_theta_GJ_1 = ROOT.TMath.Sqrt(1.0-cos_theta_GJ_1*cos_theta_GJ_1)
    sin_theta_GJ_2 = ROOT.TMath.Sqrt(1.0-cos_theta_GJ_2*cos_theta_GJ_2)

    #    print('E(tau) = %5.3f  E(vis) = %5.3f  m(tau) = %5.3f  m(vis) = %5.3f'%(tau.E(),visTau.E(),utils.tau_mass,mass_vis_tau))
    #    print('cos_theta_GJ_1 = %10.8f   cos_theta_GJ_2 = %10.8f'%(cos_theta_GJX_1,cos_theta_GJX_2)) 
    
    # orthogonal basis
    
    n1 = vistau_P.Unit()
    n3 = n1.Cross(tau_dir).Unit()
    n2 = n3.Cross(n1).Unit()

    P_perp  = tau_P.Dot(n3)
    P_paral = ROOT.TMath.Sqrt(tau_mom*tau_mom-P_perp*P_perp) 
    p_paral = []
    p_paral.append(n1*cos_theta_GJ_1+n2*sin_theta_GJ_1)
    p_paral.append(n1*cos_theta_GJ_1-n2*sin_theta_GJ_1)
    p_paral.append(n1*cos_theta_GJ_2+n2*sin_theta_GJ_2)
    p_paral.append(n1*cos_theta_GJ_2-n2*sin_theta_GJ_2)
    n_paral = p_paral[0]
    max_cos = n_paral.Dot(tau_dir)
    for i in range(1,4):
        x = p_paral[i].Dot(tau_dir)
        if x>max_cos:
            max_cos = x
            n_paral = p_paral[i]
    

    new_tau_P = n3*P_perp+n_paral*P_paral
    new_tau = ROOT.TLorentzVector()
    new_tau.SetXYZT(new_tau_P.X(),new_tau_P.Y(),new_tau_P.Z(),tau.E())
    return new_tau

def rotateToGJMax(vistau,tau):
    vistau_P = vistau.Vect()
    tau_P = tau.Vect()
    vistau_dir = vistau_P.Unit()
    tau_dir = tau_P.Unit()
    vistau_mom = vistau.P()
    tau_mom = tau.P()
    
    mass_vis_tau = vistau.M()

    theta_GJ = math.acos(max(-1.0,min(1.0,tau_dir.Dot(vistau_dir))))
    sin_GJ_max = 0.5*(utils.tau_mass*utils.tau_mass-mass_vis_tau*mass_vis_tau)/(utils.tau_mass*vistau_mom)
    theta_GJ_max = math.asin(max(-1.0,min(1.0,sin_GJ_max)))

    new_tau = ROOT.TLorentzVector()
    new_tau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())
    if theta_GJ>theta_GJ_max:
#        print('Here we are')
        n1 = vistau_dir
        n3 = n1.Cross(tau_dir).Unit()
        n2 = n3.Cross(n1).Unit()

        P_perp  = tau_P.Dot(n3)
        P_paral = math.sqrt(tau_mom*tau_mom-P_perp*P_perp) 

        n_paral_1 = (n1*math.cos(theta_GJ_max)-n2*math.sin(theta_GJ_max)).Unit()
        n_paral_2 = (n1*math.cos(theta_GJ_max)+n2*math.sin(theta_GJ_max)).Unit()
        
        new_dir_1 = (n3*P_perp+n_paral_1*P_paral).Unit()
        new_dir_2 = (n3*P_perp+n_paral_2*P_paral).Unit()
        condition = new_dir_1.Dot(tau_dir) < new_dir_2.Dot(tau_dir)
        
        new_dir = new_dir_1
        if condition:
            new_dir = new_dir_2

        new_P = tau_mom*new_dir
        new_tau = ROOT.TLorentzVector()
        new_tau.SetXYZT(new_P.X(),new_P.Y(),new_P.Z(),tau.E())
#        print('Rotated tau (x,y,z,t) = (%5.3f,%5.3f,%5.3f,%5.3f)'%(new_tau.X(),new_tau.Y(),new_tau.Z(),new_tau.T()))
    return new_tau
    
# returns full tau 4-momentum and polarimetric vector
def PolVectRho(Pt,Pi,Pi0):

    visTau = Pi+Pi0
    tau = ROOT.TLorentzVector()
    tau.SetPtEtaPhiM(Pt,visTau.Eta(),visTau.Phi(),utils.tau_mass)

    tempTau = ROOT.TLorentzVector()
    tempTau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())

    tempPi = ROOT.TLorentzVector()
    tempPi.SetXYZT(Pi.X(),Pi.Y(),Pi.Z(),Pi.T())

    tempPi0 = ROOT.TLorentzVector()
    tempPi0.SetXYZT(Pi0.X(),Pi0.Y(),Pi0.Z(),Pi0.T())

    boost = -tau.BoostVector()
    
    tempTau.Boost(boost)
    tempPi.Boost(boost)
    tempPi0.Boost(boost)

    tempQ = tempPi - tempPi0
    tempN = tempTau - tempPi - tempPi0
    tempN.SetPtEtaPhiM(tempN.Pt(),tempN.Eta(),tempN.Phi(),0.)
    
    tempPV = 2*(tempQ*tempN)*tempQ#-tempQ.M2()*tempN
    pv = ROOT.TLorentzVector()
    pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),tempPV.T())
    pv.Boost(-boost)
    
    return tau,pv


# returns corrected full tau 4-momentum and polarimetric vector in 
def PolVectA1(PV,SV,
              Pt,P1,P2,P3,
              q,q1,q2,q3):    

#    P1,P2,P3 = sortA1(Pi1,Pi2,Pi3,q1,q2,q3)    
#    print('PolVectA1 : q1 = %3.1f  q2 = %3.1f  q3 = %3.1f   q = %3.1f'%(q1,q2,q3,q))
#    print('PolVectA1 : PV(x) = %6.4f  PV(y) = %6.4f  PV(z) = %6.4f'%(PV.X(),PV.Y(),PV.Z()))
#    print('PolVectA1 : SV(x) = %6.4f  SV(y) = %6.4f  SV(z) = %6.4f'%(SV.X(),SV.Y(),SV.Z()))

    n = SV - PV
    
    Ptau = ROOT.TLorentzVector()
    Ptau.SetPtEtaPhiM(Pt,n.Eta(),n.Phi(),utils.tau_mass)
    visTau = P1+P2+P3
    tau = rotateToGJMax(visTau,Ptau)
    tempTau = ROOT.TLorentzVector()
    tempTau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())
    tempP1 = ROOT.TLorentzVector()
    tempP1.SetXYZT(P1.X(),P1.Y(),P1.Z(),P1.T())
    tempP2 = ROOT.TLorentzVector()
    tempP2.SetXYZT(P2.X(),P2.Y(),P2.Z(),P2.T())
    tempP3 = ROOT.TLorentzVector()
    tempP3.SetXYZT(P3.X(),P3.Y(),P3.Z(),P3.T())
    
    boost = -tau.BoostVector()

    tempTau.Boost(boost)
    tempP1.Boost(boost)
    tempP2.Boost(boost)
    tempP3.Boost(boost)
    
    PA1 = PolarimetricA1(tempTau,
                         tempP1,
                         tempP2,
                         tempP3,
                         q)

    tempPV = -PA1.PVC().Vect()
    pv = ROOT.TLorentzVector()
    pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),0.)
    pv.Boost(-boost)
    return tau,pv
    
# P1, P2, R1, R2 - TLorentzVectors
# firstNeg = True if muon is negative, otherwise False
# method_1, method_2 - str
# available options for method_{1,2} :
#   - Impact-Parameter
#   - Decay-Plane
#   - Decay-Plane-a1
#   - PV 
def acoCP(P1_input, P2_input, R1_input, R2_input,
          firstNeg, method_1, method_2):



    P1 = ROOT.TLorentzVector()
    P1.SetXYZT(P1_input.X(),P1_input.Y(),P1_input.Z(),P1_input.T())
    P2 = ROOT.TLorentzVector()
    P2.SetXYZT(P2_input.X(),P2_input.Y(),P2_input.Z(),P2_input.T())
    R1 = ROOT.TLorentzVector()
    R1.SetXYZT(R1_input.X(),R1_input.Y(),R1_input.Z(),R1_input.T())
    R2 = ROOT.TLorentzVector()
    R2.SetXYZT(R2_input.X(),R2_input.Y(),R2_input.Z(),R2_input.T())

    Prongsum = P1 + P2
    boost = -Prongsum.BoostVector()
    
    if method_1=='Impact-Parameter':
        temp = R1.Vect().Unit()
        R1.SetXYZT(temp.X(),temp.Y(),temp.Z(),0.)
    if method_2=='Impact-Parameter':
        temp = R2.Vect().Unit()
        R2.SetXYZT(temp.X(),temp.Y(),temp.Z(),0.)
    
    """
    if method_1=='Impact-Parameter' and method_2=='Impact-Parameter':
        print('%s -- %s'%(method_1,method_2))
        print('before boost in P1+P2 frame -> ')
        print('P1 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(P1.X(),P1.Y(),P1.Z(),P1.T()))
        print('P2 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(P2.X(),P2.Y(),P2.Z(),P2.T()))
        print('R1 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(R1.X(),R1.Y(),R1.Z(),R1.T()))
        print('R2 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(R2.X(),R2.Y(),R2.Z(),R2.T()))
    """
    
    P1.Boost(boost)
    P2.Boost(boost)
    R1.Boost(boost)
    R2.Boost(boost)
    
    y1 = 1
    y2 = 1
    
    if method_1=='Decay-Plane': # rho(+/-) -> pi(+/-) + pi0
        y1 = (P1.E()-R1.E())/(P1.E()+R1.E())
    elif method_1=='Decay-Plane-a1': # a1(+/-) -> rho0 + pi(+/-)
        y1 = (R1.E()-P1.E())/(P1.E()+R1.E())
        
    if method_2=='Decay-Plane': # neutral pion
        y2 = (P2.E()-R2.E())/(P2.E()+R2.E())
    elif method_2=='Decay-Plane-a1': # a1(+/-) -> rho0 + pi(+/-)
        y2 = (R2.E()-P2.E())/(P2.E()+R2.E())

    y = y1*y2
    
    vecP1 = P1.Vect().Unit()
    vecP2 = P2.Vect().Unit()
    vecR1 = R1.Vect()
    vecR2 = R2.Vect()

    R1transv = vecR1 - vecP1*(vecP1.Dot(vecR1))
    R2transv = vecR2 - vecP2*(vecP2.Dot(vecR2))

    n1 = R1transv.Unit()
    n2 = R2transv.Unit()

    cos_phi = n1.Dot(n2)
    acop = math.acos(cos_phi)

    """
    if method_1=='Impact-Parameter' and method_2=='Impact-Parameter':
        print('after boost in P1+P2 frame -> ')
        print('P1 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(P1.X(),P1.Y(),P1.Z(),P1.T()))
        print('P2 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(P2.X(),P2.Y(),P2.Z(),P2.T()))
        print('R1 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(R1.X(),R1.Y(),R1.Z(),R1.T()))
        print('R2 (x,y,z,t)     = (%5.3f,%5.3f,%5.3f,%5.3f)'%(R2.X(),R2.Y(),R2.Z(),R2.T()))
        print('R1transv (x,y,z) = (%5.3f,%5.3f,%5.3f)'%(R1transv.X(),R1transv.Y(),R1transv.Z()))
        print('R2transv (x,y,z) = (%5.3f,%5.3f,%5.3f)'%(R2transv.X(),R2transv.Y(),R2transv.Z()))
        print('n1 (x,y,z) = (%5.3f,%5.3f,%5.3f) len = %5.3f'%(n1.X(),n1.Y(),n1.Z(),n1.Mag()))
        print('n2 (x,y,z) = (%5.3f,%5.3f,%5.3f) len = %5.3f'%(n2.X(),n2.Y(),n2.Z(),n1.Mag()))
        print('cos(phi) = %5.4f'%(cos_phi))
    """

    sign = vecP2.Dot(n1.Cross(n2))
    
    if firstNeg:
        sign = vecP1.Dot(n2.Cross(n1))

    if sign<0:
        acop = 2.0*ROOT.TMath.Pi() - acop

    if y<0:
        acop = acop + ROOT.TMath.Pi()
        if acop>2.*ROOT.TMath.Pi():
            acop = acop - 2.*ROOT.TMath.Pi()

    return acop

            
def CosAlpha(pt,eta,phi,mass,Rx,Ry,Rz,prn):

    P = ROOT.TLorentzVector()
    P.SetPtEtaPhiM(pt,eta,phi,mass)
    ez = ROOT.TVector3(0.,0.,1.)
    p  = ROOT.TVector3(P.Px(),P.Py(),P.Pz())
    n  = ROOT.TVector3(Rx,Ry,Rz)

    ez = ez.Unit();
    p = p.Unit();
    n = n.Unit();

    denom = ez.Cross(p).Mag() * n.Cross(p).Mag()

    if denom > 0.: 
        cosa = abs( ez.Cross(p).Dot(n.Cross(p)) / denom )
    else:
        if prn:
            print('px = %3.1f  py = %3.1f  pz = %3.1f'%(p.X(),p.Y(),p.Z()))
            print('nx = %3.1f  ny = %3.1f  nz = %3.1f'%(Rx,Ry,Rz))
            print('')
        cosa = 1.0/math.sqrt(2.0)
    
    return cosa


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
        
        pt_2        = np.zeros(1,dtype=np.float64)
        eta_2       = np.zeros(1,dtype=np.float64)
        phi_2       = np.zeros(1,dtype=np.float64)
        iso_2       = np.zeros(1,dtype=np.float64)
        mass_2      = np.zeros(1,dtype=np.float64)
        charge_2    = np.zeros(1,dtype=np.int16)

        mt_2        = np.zeros(1,dtype=np.float64)
        ipsig_2     = np.zeros(1,dtype=np.float64)
        
        decayMode_2 = np.zeros(1,dtype=np.int64)
        decayModePNet_2 = np.zeros(1,dtype=np.int64)

        pt_1_FastMTT = np.zeros(1,dtype=np.float64)
        pt_2_FastMTT = np.zeros(1,dtype=np.float64)
        
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

        alphaAngle_2 = np.zeros(1,dtype=np.float64)

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

        sv_x_2       = np.zeros(1,dtype=np.float64)
        sv_y_2       = np.zeros(1,dtype=np.float64)
        sv_z_2       = np.zeros(1,dtype=np.float64)
        
        # booleans
        trg_lep     = np.zeros(1,dtype='?')
        trg_cross   = np.zeros(1,dtype='?')
        os          = np.zeros(1,dtype='?')
        hasRefitSV_2 = np.zeros(1,dtype='?')
        
        # integers
        idDeepTau2018v2p5VSe_2   = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSmu_2  = np.zeros(1,dtype=np.int64)
        idDeepTau2018v2p5VSjet_2 = np.zeros(1,dtype=np.int64)

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

        if channel=='mt' or channel=='et':
            
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

            tree.SetBranchAddress('hasRefitSV_2',hasRefitSV_2)
            tree.SetBranchAddress('FastMTT_pt_1_constraint',pt_1_FastMTT)
            tree.SetBranchAddress('FastMTT_pt_2_constraint',pt_2_FastMTT)

            tree.SetBranchAddress('sv_x_2',sv_x_2)
            tree.SetBranchAddress('sv_y_2',sv_y_2)
            tree.SetBranchAddress('sv_z_2',sv_z_2)

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
            
        if channel=='mt' or channel=='et':
            tree.SetBranchAddress('decayMode_2',decayMode_2)
            tree.SetBranchAddress('decayModePNet_2',decayModePNet_2)

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
            
        # integers
        if channel=='mt' or channel=='et':
            tree.SetBranchAddress('idDeepTau2018v2p5VSe_2',idDeepTau2018v2p5VSe_2)
            tree.SetBranchAddress('idDeepTau2018v2p5VSmu_2',idDeepTau2018v2p5VSmu_2)
            tree.SetBranchAddress('idDeepTau2018v2p5VSjet_2',idDeepTau2018v2p5VSjet_2)

        if self.ismc:
            tree.SetBranchAddress("genPartFlav_1",genmatch_1)
            tree.SetBranchAddress("genPartFlav_2",genmatch_2)
            
        nentries = tree.GetEntries()

        # run over entries
        for entry in range(0,nentries):

            if entry%100000==0:
                print('processed %1i out of %1i events'%(entry,nentries))
            tree.GetEntry(entry)

            # prevent large weights
            if weight[0]>100: continue
            
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
                variables['aco_DM0_plus'] = -9999.
                variables['aco_DM0_minus'] = -9999.
                variables['aco_DM0'] = -9999.
            
                variables['aco_DM1_plus'] = -9999.
                variables['aco_DM1_minus'] = -9999.
                variables['aco_DM1'] = -9999.
                variables['aco_DM1_PV'] = -9999.

                variables['aco_DM10_plus'] = -9999.
                variables['aco_DM10_minus'] = -9999.
                variables['aco_DM10'] = -9999.
                variables['aco_DM10_PV'] = -9999.
                variables['aco_DM10_DP'] = -9999.
            
                variables['alpha_DM0'] = -9999.
                variables['alpha_DM1'] = -9999.
                variables['alpha_DM10'] = -9999.

                alpha = -9999.
                cosa = CosAlpha(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi,ip_x_2[0],ip_y_2[0],ip_z_2[0],False)
                alpha_IP = ROOT.TMath.ACos(cosa)
                aco = -9999.
                pv = ROOT.TLorentzVector()
                aco_a1 = -9999.
                
                if decayModePNet_2[0]==0 and abs(ipsig_2[0])>cuts.ipsigLepCut:
                    alpha = alpha_IP
                    variables['alpha_DM0'] = RadToDeg * alpha
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
                    aco = acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','Impact-Parameter')
                    P1.SetPtEtaPhiM(pt_1[0],eta_1[0],phi_1[0],massLep)
                    P2.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                    aco_old = acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','Impact-Parameter')
                    if math.isnan(aco): aco = -9999.
                    """
                    variables['aco_DM0'] = RadToDeg*aco_lep_pi[0]
                    if alpha > pi_over_4:
                        variables['aco_DM0_plus'] = RadToDeg*aco_lep_pi[0]
                    else:
                        variables['aco_DM0_minus'] = RadToDeg*aco_lep_pi[0]
#                    print('PNetDM %2i -> alpha = %5.3f : %5.3f -- phi(CP) = %5.3f : %5.3f'%(decayModePNet_2[0],alpha,alphaAngle_2[0],aco_old,aco))
#                    print('')
                elif decayModePNet_2[0]==1 and decayMode_2[0]==1 and pi0_pt_2[0]>1:
                    Pi0 = ROOT.TLorentzVector()
                    Pi0.SetPtEtaPhiM(pi0_pt_2[0],pi0_eta_2[0],pi0_phi_2[0],massPi0)
                    Pi = ROOT.TLorentzVector()
                    Pi.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                    magPi0 = Pi0.Vect().Mag()
                    magPi = Pi.Vect().Mag()
                    cosa = CosAlpha(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi,Pi0.Px(),Pi0.Py(),Pi0.Pz(),True)
                    alpha = ROOT.TMath.ACos(cosa)
                    variables['alpha_DM1'] = RadToDeg * alpha
                    P1 = ROOT.TLorentzVector()
                    P1.SetPtEtaPhiM(pt_1_FastMTT[0],eta_1[0],phi_1[0],massLep)
                    R1 = ROOT.TLorentzVector()
                    R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0.0)
                    P2,R2 = PolVectRho(pt_2_FastMTT[0],Pi,Pi0)
                    firstNeg = charge_1[0] < 0.
                    aco = acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')
                    if math.isnan(aco): aco = -9999.
                    variables['aco_DM1_PV'] = RadToDeg*aco
                    variables['aco_DM1'] = RadToDeg*aco_lep_rho[0]
                    if alpha > pi_over_4:
                        variables['aco_DM1_plus'] = RadToDeg*aco_lep_rho[0]
                    else:
                        variables['aco_DM1_minus'] = RadToDeg*aco_lep_rho[0]
#                    print('PNetDM %2i -> alpha = %5.3f : %5.3f -- phi(CP) = %5.3f : %5.3f'%(decayModePNet_2[0],alpha_IP,alphaAngle_2[0],aco_lep_rho[0],aco))
#                    print('')
                elif decayModePNet_2[0]==10 and hasRefitSV_2[0] and pt_2_FastMTT[0]>10.:
                    """
                    print('run = %1i   event = %1i'%(run[0],event[0]))
                    print('PNetDM = %2i  hasRefitSV_2 = %1i'%(decayModePNet_2[0],hasRefitSV_2[0]))
                    print('muon   (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pt_1[0],eta_1[0],phi_1[0],massLep,charge_2[0]))
                    print('muon impact par. (x,y,z) = (%8.6f,%8.6f,%8.6f)'%(ip_x_1[0],ip_y_1[0],ip_z_1[0]))
                    print('FastMTT_pt_2_constraint  = %5.1f'%(pt_2_FastMTT[0]))
                    print('pi1    (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],pi_mass_2[0],int(pi_charge_2[0])))
                    print('pi2    (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],pi2_mass_2[0],int(pi2_charge_2[0])))
                    print('pi3    (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f) q = %2i'%(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],pi3_mass_2[0],int(pi3_charge_2[0])))
                    """
                    sv_pv = ROOT.TVector3(sv_x_2[0]-PVBS_x[0],sv_y_2[0]-PVBS_y[0],sv_z_2[0]-PVBS_z[0])
                    Pi1 = ROOT.TLorentzVector()
                    Pi1.SetPtEtaPhiM(pi_pt_2[0],pi_eta_2[0],pi_phi_2[0],massPi)
                    Pi2 = ROOT.TLorentzVector()
                    Pi2.SetPtEtaPhiM(pi2_pt_2[0],pi2_eta_2[0],pi2_phi_2[0],massPi)
                    Pi3 = ROOT.TLorentzVector()
                    Pi3.SetPtEtaPhiM(pi3_pt_2[0],pi3_eta_2[0],pi3_phi_2[0],massPi)
                    lv_a1 = Pi1 + Pi2 + Pi3
                    cosa = CosAlpha(lv_a1.Pt(),lv_a1.Eta(),lv_a1.Phi(),lv_a1.M(),sv_pv.X(),sv_pv.Y(),sv_pv.Z(),False)
                    alpha = ROOT.TMath.ACos(cosa)
                    variables['alpha_DM10'] = RadToDeg * alpha

                    PV_Mag = math.sqrt(PVBS_x[0]*PVBS_x[0]+PVBS_y[0]*PVBS_y[0]+PVBS_z[0]*PVBS_z[0])
                    SV_Mag = math.sqrt(sv_x_2[0]*sv_x_2[0]+sv_y_2[0]*sv_y_2[0]+sv_z_2[0]*sv_z_2[0])
                    P1 = ROOT.TLorentzVector()
                    P1.SetPtEtaPhiM(pt_1_FastMTT[0],eta_1[0],phi_1[0],massLep)
                    R1 = ROOT.TLorentzVector()
                    R1.SetXYZT(ip_x_1[0],ip_y_1[0],ip_z_1[0],0.0)
                    
                    PV = ROOT.TLorentzVector()
                    PV.SetXYZT(PVBS_x[0],PVBS_y[0],PVBS_z[0],PV_Mag)
                    SV = ROOT.TLorentzVector()
                    SV.SetXYZT(sv_x_2[0],sv_y_2[0],sv_z_2[0],SV_Mag)
                    P_os,P_ss1,P_ss2 = sortA1(Pi1,Pi2,Pi3,pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])
                    P2,R2 = PolVectA1(PV,SV,
                                      pt_2_FastMTT[0],P_os,P_ss1,P_ss2,
                                      charge_2[0],pi_charge_2[0],pi2_charge_2[0],pi3_charge_2[0])
                    firstNeg = charge_1[0] < 0.
                    aco = acoCP(P1,P2,R1,R2,firstNeg,'Impact-Parameter','PV')

                    if math.isnan(aco): aco = -9999.
                    variables['aco_DM10_PV'] = RadToDeg*aco
                    variables['aco_DM10_DP'] = RadToDeg*aco_lep_a1[0]
                    variables['aco_DM10'] = RadToDeg*aco_lep_a1_FastMTT[0]
                    if alpha > pi_over_4:
                        variables['aco_DM10_plus'] = RadToDeg*aco_lep_a1_FastMTT[0]
                    else:
                        variables['aco_DM10_minus'] = RadToDeg*aco_lep_a1_FastMTT[0]

#                    print('aco_mu_a1                        = %5.3f (tuple) : %5.3f (my code)'%(aco_lep_a1[0],aco_a1))
#                    print('aco_mu_a1_FASTMTT_MassConstraint = %5.3f (tuple) : %5.3f (my code)'%(aco_lep_a1_FastMTT[0],aco))
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
