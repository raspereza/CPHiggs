import ROOT
from array import array
import os

class ScaleFactor:

    def __init__(self,**kwargs):
        self.filename = kwargs.get('filename','None')
        self.labelSF = kwargs.get('label','promptSF')
        if os.path.isfile(self.filename):
            print('Loading file with scale factors for IPSig cut : %s'%(self.filename))
        else:
            print('No file %s is found'%(self.filename))
            print('Quitting')
            exit()
        self.sfFile = ROOT.TFile(self.filename,'read')
        self.eff_data = self.sfFile.Get('effData')
        self.eff_mc = self.sfFile.Get('effMC')
        self.nbinsPt = self.eff_data.GetNbinsX()
        self.nbinsEta = self.eff_data.GetNbinsY()
        self.minPt = self.eff_data.GetXaxis().GetBinLowEdge(1)
        self.maxPt = self.eff_data.GetXaxis().GetBinLowEdge(self.nbinsPt+1)
        self.minEta = self.eff_data.GetYaxis().GetBinLowEdge(1)
        self.maxEta = self.eff_data.GetYaxis().GetBinLowEdge(self.nbinsEta+1)
        self.binsEta = []
        self.binsEta.append(self.minEta)
        for ib in range(1,self.nbinsEta+1):
            self.binsEta.append(self.eff_data.GetYaxis().GetBinCenter(ib))
        self.binsEta.append(self.maxEta)
        self.nbinsExtrapEta = self.nbinsEta + 1
        self.EtaBinning = ROOT.TH1D(self.labelSF,"",self.nbinsExtrapEta,array('d',list(self.binsEta)))
        print('nbinsPt  = %1i'%(self.nbinsPt))
        print('nbinsEta = %1i'%(self.nbinsEta))
        print('minPt    = %2.0f'%(self.minPt))
        print('maxPt    = %2.0f'%(self.maxPt))
        print('minEta   = %4.2f'%(self.minEta))
        print('maxEta   = %4.2f'%(self.maxEta))
        self.hfits = {}
        for i in range(1,self.nbinsEta+1):
            istr = '%1i'%(i)
            self.hfits[i] = self.sfFile.Get('hfit_binEta'+istr)
        
        
    def getSF(self,pt,eta,**kwargs):
        absEta = abs(eta)
        Eta = ROOT.TMath.Min(absEta,self.maxEta-0.001)
        Pt = ROOT.TMath.Max(self.minPt+0.01,ROOT.TMath.Min(pt,self.maxPt-0.01))
        binEta = self.eff_data.GetYaxis().FindBin(Eta)
        bin_id = self.EtaBinning.GetXaxis().FindBin(Eta)
        sf = 1.0
        if bin_id==1:
            bin_pt = self.hfits[1].GetXaxis().FindBin(Pt)
            sf = self.hfits[1].GetBinContent(bin_pt)
        elif bin_id==self.nbinsExtrapEta:
            bin_pt = self.hfits[self.nbinsEta].GetXaxis().FindBin(Pt)
            sf = self.hfits[self.nbinsEta].GetBinContent(bin_pt)
        else:
            bin1 = bin_id-1
            bin2 = bin_id
            x1 = self.EtaBinning.GetXaxis().GetBinLowEdge(bin_id)
            x2 = self.EtaBinning.GetXaxis().GetBinLowEdge(bin_id+1)
            dx = x2 - x1
            bin_pt = self.hfits[bin1].GetXaxis().FindBin(Pt)
            sf1 = self.hfits[bin1].GetBinContent(bin_pt)
            sf2 = self.hfits[bin2].GetBinContent(bin_pt)
            dsf = sf2 - sf1
            dsf_dx = dsf/dx
            sf = sf1 + dsf_dx*(Eta-x1)
            
        return sf
        
    def getEffData(self,pt,eta):
        ptX = pt
        etaX = abs(eta)
        if ptX<self.minPt: ptX = self.minPt+0.01
        if ptX>self.maxPt: ptX = self.maxPt-0.01
        if etaX<self.minEta: etaX = self.minEta+0.01
        if etaX>self.maxEta: etaX = self.maxEta-0.01
        eff = self.eff_data.GetBinContent(self.eff_data.FindBin(ptX,etaX))
        return eff

    def getEffMC(self,pt,eta):
        ptX = pt
        etaX = abs(eta)
        if ptX<self.minPt: ptX = self.minPt+0.01
        if ptX>self.maxPt: ptX = self.maxPt-0.01
        if etaX<self.minEta: etaX = self.minEta+0.01
        if etaX>self.maxEta: etaX = self.maxEta-0.01
        eff = self.eff_mc.GetBinContent(self.eff_mc.FindBin(ptX,etaX))
        return eff

    def getBinnedSF(self,pt,eta):
        ptX = pt
        etaX = abs(eta)
        if ptX<self.minPt: ptX = self.minPt+0.01
        if ptX>self.maxPt: ptX = self.maxPt-0.01
        if etaX<self.minEta: etaX = self.minEta+0.01
        if etaX>self.maxEta: etaX = self.maxEta-0.01
        effData = self.eff_data.GetBinContent(self.eff_data.FindBin(ptX,etaX))
        effMC = self.eff_mc.GetBinContent(self.eff_mc.FindBin(ptX,etaX))
        sf = 1.0
        if effData>0 and effMC>0:
            sf = effData/effMC
        return sf
    

    
