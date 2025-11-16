import ROOT
from array import array
import os

class FakeFactor:

    def __init__(self,**kwargs):
        self.filename = kwargs.get('filename','None')
        print('FakeFactor class ->')
        if os.path.isfile(self.filename):
            print('Loading file with scale factors for IPSig cut : %s'%(self.filename))
        else:
            print('No file %s is found'%(self.filename))
            print('Quitting')
            exit()
        self.ff_file = ROOT.TFile(self.filename,'read')
        self.typ_labels = ['qcd','wj','ss_antiiso','mc_top','mc_wj']
        self.dm_labels = ['pi','rho','a1_1pr','a1_3pr']
        self.eta_labels = ['barrel','endcap','all']
        self.njets_labels = ['njets0','njets1','njets2']
        self.sys_labels = ['nom','up','down']
        self.ff_funcs = {}
        for typ in self.typ_labels:
            for dm in self.dm_labels:
                for njets in self.njets_labels:
                    for eta in self.eta_labels:
                        for sys in self.sys_labels:
                            suffix = '%s_%s_%s_%s_%s'%(typ,dm,njets,eta,sys)
                            name = 'fitFunc_'+suffix
                            self.ff_funcs[suffix] = self.ff_file.Get(name)
        self.ptMin = 20.
        self.ptMax = 120.
        
    def getFF(self,pt,**kwargs):
        eta = kwargs.get('eta','barrel')
        dm = kwargs.get('dm','pi')
        typ = kwargs.get('typ','qcd')
        njets = kwargs.get('njets','njets0')
        sys = kwargs.get('sys','nom')
        ff = 1
        if eta not in self.eta_labels:
            print('unknown eta specified: %s. returning 1.'%(eta))
            return ff
        if dm not in self.dm_labels:
            print('unknown decay mode specified: %s. returning 1.'%(dm))
            return ff
        if njets not in self.njets_labels:
            print('unknown njets specified %s. returning 1.'%(njets))
            return ff
        if typ not in self.typ_labels:
            print('unknown typ specified %s. returning 1.'%(typ))
            return ff
        if sys not in self.sys_labels:
            print('unknown sys specified %s. returning 1.'%(sys))
            return ff

        Pt = pt
        if Pt>self.ptMax: Pt = self.ptMax
        if Pt<self.ptMin: Pt = self.ptMin
        suffix = '%s_%s_%s_%s_%s'%(typ,dm,njets,eta,sys)
        ff = self.ff_funcs[suffix].Eval(Pt)
        if ff<0: ff=0
        return ff
        
