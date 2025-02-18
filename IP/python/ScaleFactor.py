import ROOT
import os

class ScaleFactor:

    def __init__(self,**kwargs):
        self.filename = kwargs.get('filename','None')
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
        print('nbinsPt  = %1i'%(self.nbinsPt))
        print('nbinsEta = %1i'%(self.nbinsEta))
        print('minPt    = %2.0f'%(self.minPt))
        print('maxPt    = %2.0f'%(self.maxPt))
        print('minEta   = %4.2f'%(self.minEta))
        print('maxEta   = %4.2f'%(self.maxEta))
        
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

    def getSF(self,pt,eta):
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
    

    
