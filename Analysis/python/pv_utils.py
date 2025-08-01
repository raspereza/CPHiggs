### utilities to compute phi(CP)
import ROOT
import math
from array import array
import numpy as np
import os
import CPHiggs.Analysis.utils as utils
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1

####################################################################
# TauDirRho : routine returns direction of fully reconstructed tau
#             in tau->rho+nu decay channel 
# Inputs ->
# Pt -float : tau pT from fastMTT
# Pvis - TLorentzVector : visible 4-momentum of tau (pi(+/-) + pi0
# Ppion - TLorentzVector : visible 4-momentum of charged pion
# IP - TVector3 : impact parameter vector
# PGen - TLorentzVector : generator tau momentum
#                         for exploratory MC study only,
#                         assumes perfect knowledge of neutrino momentum
#                         (any dummy TLorentzVector can be passed for data)
# method - string : options ->
#                   'reco'      - reconstruction info is used,
#                                 method uses fastMTT and decay
#                                 length information
#                   'collinear' - reconstruction into is used,
#                                 neutrino momentum is assumed
#                                 to be collinear with tau momentum
#                   'reco_gen'  - reconstruction info is used
#                                  but two-fold ambiguity in
#                                  solutions is resolved by
#                                  mathcing to generator tau momentum 
#                   'gen'       - solely generator info is used
#                                 (for exploratory studies)
#
# Output ->
# dir_tau - TVector3 : direction of tau momentum (unit vector)
#########################################################
def TauDirRho(Pt,Pvis,Ppion,IP,PGen,method):

    Ptau = ROOT.TLorentzVector()
    Ptau.SetPtEtaPhiM(Pt,Pvis.Eta(),Pvis.Phi(),utils.tau_mass)
    if method=='collinear':
        tau_dir = Ptau.Vect().Unit()
        return tau_dir
    elif method=='gen':
        tau_dir = PGen.Vect().Unit()
        return tau_dir

    
    n = Pvis.Vect().Unit()
    p = Ppion.Vect().Unit()
    r = IP.Unit()    
    l = p.Cross(r).Unit()
    prod = l.Dot(n)
    if prod<0:
        l = -l

    dirGen = PGen.Vect().Unit()
        
#    prod = l.Dot(n)
#    print('p*l = %6.4f'%(prod))
        
    vistau_mass = Pvis.M()
    vistau_mom = Pvis.P()
    vistau_E = Pvis.E()
    tau_mom = Ptau.P()
    tau_E = Ptau.E()
    
#    if Ptau.P()==0:
#        print('Ptau = (%8.6f,%8.6f,%8.6f)'%(Ptau.X(),Ptau.Y(),Ptau.Z()))
#    if Pvis.P()==0:
#        print('Pvis = (%8.6f,%8.6f,%8.6f)'%(Pvis.X(),Pvis.Y(),Pvis.Z()))
    cosThetaGJ = max(-1.0,min(1.0,0.5*(2*tau_E*vistau_E-utils.tau_mass*utils.tau_mass-vistau_mass*vistau_mass)/(tau_mom*vistau_mom)))
    sinThetaGJ = math.sqrt(1.0-cosThetaGJ*cosThetaGJ)
#    print('cosThetaGJ = %10.8f    sinThetaGJ = %10.8f'%(cosThetaGJ,sinThetaGJ))
    
    A = l.Z()*cosThetaGJ/(n.Y()*l.Z()-n.Z()*l.Y())
    B = (n.Z()*l.X()-n.X()*l.Z())/(n.Y()*l.Z()-n.Z()*l.Y())
    C = -A*l.Y()/l.Z()
    D = -l.X()/l.Z() - B*l.Y()/l.Z()

    P = (A*B+C*D)/(B*B+D*D+1)
    Q = (A*A+C*C-1)/(B*B+D*D+1)

    Det = P*P-Q

    # initialization of tau direction
    q = ROOT.TVector3()
#    print('D = %10.8f'%(D))
    if Det<0:
        sin_GJ_max = 0.5*(utils.tau_mass*utils.tau_mass-vistau_mass*vistau_mass)/(utils.tau_mass*vistau_mom)
        cos_GJ_max = ROOT.TMath.Sqrt(1-sin_GJ_max*sin_GJ_max)
        
        s = l.Cross(n).Unit()
        m = n.Cross(s).Unit()
        q = cos_GJ_max*n-sin_GJ_max*m
        
#        print('Det = %10.8f'%(Det))
    else:
        q1_x = -P+math.sqrt(Det)
        q1_y = A + B*q1_x
        q1_z = C + D*q1_x

        q2_x = -P-math.sqrt(Det)
        q2_y = A + B*q2_x
        q2_z = C + D*q2_x
        
        q1 = ROOT.TVector3(q1_x,q1_y,q1_z)
        q2 = ROOT.TVector3(q2_x,q2_y,q2_z)
#        mag1 = math.sqrt(q1_x*q1_x+q1_y*q1_y+q1_z*q1_z)
#        mag2 = math.sqrt(q2_x*q2_x+q2_y*q2_y+q2_z*q2_z)
#        print('m1 = %10.8f    m2 = %10.8f'%(q1.Mag(),q2.Mag()))

        q1 = q1.Unit()
        q2 = q2.Unit()
        px1 = q1.Dot(dirGen)
        px2 = q2.Dot(dirGen)
        if px1>px2:
            q = q1
        else:
            q = q2
        if method=='reco':
            cosOmega1 = q1.Dot(r)        
            cosOmega2 = q2.Dot(r)
            IP_mag = IP.Mag()
            L1 = IP_mag/cosOmega1
            L2 = IP_mag/cosOmega2
            if L1>0 and L2>0:
                beta = utils.tau_mass/(utils.ctau*Ptau.P())
                P1 = 1.0-ROOT.TMath.Exp(-L2*beta)
                P2 = ROOT.TMath.Exp(-L1*beta)
                if P1>P2:
                    q = q1
                else:
                    q = q2
#                s = l.Cross(n).Unit()
#                m = n.Cross(s).Unit()
#                q = cosThetaGJ*n-sinThetaGJ*m
            else:
                if L1>L2:
                    q = q1
                else:
                    q = q2

    dir_tau = q.Unit()
#    cosTJ = dir_tau.Dot(n)
#    perp = dir_tau.Dot(l)
#    print('q*n = %10.8f   q*l = %10.8f'%(cosTJ,perp))
    
    return dir_tau
    
################################################################
# PolVectRho : routine returns full tau 4-momentum and
# polarimetric vector in tau->rho+nu decay mode
# Inputs ->
# Pt  - float : tau pT from fastMTT
# Pi  - TLorentzVector : 4-momentum of charged pion
# Pi0 - TLorentzVector : 4-momentum of neutral pion
# IP  - TLorentzVector : IP vector of charged pion
# PGen - TLorentzVector : generator 4-momentum of tau
#                         (for exploratory studies)
# method - string : options ->
#                   'reco'      - reconstruction,
#                                 method uses fastMTT and decay
#                                 length information
#                   'collinear' - reconstruction,
#                                 neutrino momentum is assumed
#                                 to be collinear with tau momentum
#                   'reco_gen'  - reconstruction,
#                                 two-fold ambiguity in
#                                 solutions is resolved by
#                                 mathcing to generator tau momentum 
#                   'gen'       - solely generator info is used
#                                 (for exploratory studies)
#
# Outputs ->
# Tau4P - TLorentzVector : full tau momentum
# pv    - TLorentzVector : polarimetric vector
#################################################################
def PolVectRho(Pt,Pi,Pi0,IP,PGen,method):

    visTau = Pi+Pi0

    tauP4 = ROOT.TLorentzVector()
    if method=='gen':
        tauP4.SetXYZT(PGen.X(),PGen.Y(),PGen.Z(),PGen.T())
    else:
        tau = ROOT.TLorentzVector()
        tau.SetPtEtaPhiM(Pt,visTau.Eta(),visTau.Phi(),utils.tau_mass)
        tau_mom = tau.P()
        dir_tau = TauDirRho(Pt,visTau,Pi,IP,PGen,method)
        tauP4.SetXYZT(tau_mom*dir_tau.X(),tau_mom*dir_tau.Y(),tau_mom*dir_tau.Z(),tau.E())
        
    tempPi = ROOT.TLorentzVector()
    tempPi.SetXYZT(Pi.X(),Pi.Y(),Pi.Z(),Pi.T())

    tempPi0 = ROOT.TLorentzVector()
    tempPi0.SetXYZT(Pi0.X(),Pi0.Y(),Pi0.Z(),Pi0.T())

#    boost = -tau4P.BoostVector()
#    tau4P.Boost(boost)
#    tempPi.Boost(boost)
#    tempPi0.Boost(boost)

    tempQ = tempPi - tempPi0
    NN = tauP4 - tempPi - tempPi0
    tempN = ROOT.TLorentzVector()
#    print('neutrino P = %5.3f    M = %10.8f'%(NN.P(),NN.M()))
#    print('')
    tempN.SetPtEtaPhiM(NN.Pt(),NN.Eta(),NN.Phi(),0.)

    X1 = tempQ.E()*tempN.E()-tempQ.Px()*tempN.Px()-tempQ.Py()*tempN.Py()-tempQ.Pz()*tempN.Pz()
    X2 = tempQ.E()*tempQ.E()-tempQ.Px()*tempQ.Px()-tempQ.Py()*tempQ.Py()-tempQ.Pz()*tempQ.Pz()
    tempPV = X1*tempQ-X2*tempN

    pv = ROOT.TLorentzVector()
    pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),tempPV.T())
#    pv.Boost(-boost)
#    tau4P.Boost(-boost)
    
    return tauP4,pv


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


#    print('a1 sorted -> ')
#    print('pi_os  (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p1.Pt(),p1.Eta(),p1.Phi(),p1.M()))
#    print('pi_ss1 (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p2.Pt(),p2.Eta(),p2.Phi(),p2.M()))
#    print('pi_ss2 (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p3.Pt(),p3.Eta(),p3.Phi(),p3.M()))
#    print('M(os_ss1) = %5.4f -- M(os_ss2) = %5.4f -- m(rho) = %5.4f'%((p1+p2).M(),(p1+p3).M(),utils.rho_mass))

    
    return p1,p2,p3


###################################################
# Rotation of initial tau momentum
# (along visible momentum) towards reconstructed
# direction (e.g. SV-PV) to match theta_GJ
###################################################
def rotateToGJ(visTau,tau):

    vistau_P = visTau.Vect()
    tau_P = tau.Vect()
    mass_vis_tau = visTau.M()
    
    tau_mom = tau.P()
    vistau_mom = visTau.P()

    tau_dir = tau_P.Unit()

    cos_theta_GJX = 0.5*(2.0*tau.E()*visTau.E()-utils.tau_mass*utils.tau_mass-mass_vis_tau*mass_vis_tau)/(tau_mom*vistau_mom)
    cos_theta_GJ = max(-1.0,min(1.0,cos_theta_GJX))
    sin_theta_GJ = ROOT.TMath.Sqrt(1.0-cos_theta_GJ*cos_theta_GJ)

    #    print('E(tau) = %5.3f  E(vis) = %5.3f  m(tau) = %5.3f  m(vis) = %5.3f'%(tau.E(),visTau.E(),utils.tau_mass,mass_vis_tau))
    #    print('cos_theta_GJ_1 = %10.8f   cos_theta_GJ_2 = %10.8f'%(cos_theta_GJX_1,cos_theta_GJX_2)) 
    
    # orthogonal basis
    
    n1 = vistau_P.Unit()
    n3 = n1.Cross(tau_dir).Unit()
    n2 = n3.Cross(n1).Unit()

    Pperp  = tau_P.Dot(n3)
    Pparal = ROOT.TMath.Sqrt(tau_mom*tau_mom-Pperp*Pperp) 
    p1 = n3*Pperp + (n1*cos_theta_GJ+n2*sin_theta_GJ)*Pparal
    p2 = n3*Pperp + (n1*cos_theta_GJ-n2*sin_theta_GJ)*Pparal
    cos1 = p1.Dot(tau_dir)
    cos2 = p2.Dot(tau_dir)
    newTauP = p1
    if cos2>cos1:
        newTauP = p2
    
    new_tau = ROOT.TLorentzVector()
    new_tau.SetXYZT(newTauP.X(),newTauP.Y(),newTauP.Z(),tau.E())
    return new_tau


###########################################################
# Rotation of visible tau momentum towards reconstructed  #
# direction (e.g. SV-PV) to match theta_GJmax             #
# vistau - 4-momentum of visible tau decay products       #
# tau - full 4-momentum of tau                            #
###########################################################
def rotateToGJMax(vistau,tau):
    vistau_P = vistau.Vect()
    tau_P = tau.Vect()
    vistau_dir = vistau_P.Unit()
    tau_dir = tau_P.Unit()
    vistau_mom = vistau.P()
    tau_mom = tau.P()
    tau_E = tau.E()
    vistau_E = tau.E()
    mass_vis_tau = vistau.M()

    cos_GJ = max(-1.0,min(1.0,vistau_dir.Dot(tau_dir)))
    theta_GJ = math.acos(cos_GJ) 
    sin_GJ_max = 0.5*(utils.tau_mass*utils.tau_mass-mass_vis_tau*mass_vis_tau)/(utils.tau_mass*vistau_mom)
    theta_GJ_max = math.asin(max(-1.0,min(1.0,sin_GJ_max)))

    new_tau = ROOT.TLorentzVector()
    new_tau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())
#    print('')
#    print('Routine rotateToGJMax ->')
    if theta_GJ>theta_GJ_max:
#       print('theta_GJ (%5.3f) > theta_GJ_max (%5.3f) -> rotation of tau 3-momentum'%(theta_GJ,theta_GJ_max))
        n1 = vistau_dir
        n3 = n1.Cross(tau_dir).Unit()
        n2 = n3.Cross(n1).Unit()
        
        P_perp  = tau_P.Dot(n3)
        P_paral = math.sqrt(tau_mom*tau_mom-P_perp*P_perp) 
        
        #            cos_theta_GJX = 0.5*(tau_E*vistau_E-utils.tau_mass*utils.tau_mass-mass_vis_tau*mass_vis_tau)/(tau_mom*vistau_mom)
        #            cos_theta_GJ = max(-1.0,min(1.0,cos_theta_GJX))
        #            sin_theta_GJ = math.sqrt(1-cos_theta_GJ*cos_theta_GJ)
        #            n_paral_1 = (n1*cos_theta_GJ-n2*sin_theta_GJ).Unit()
        #            n_paral_2 = (n1*cos_theta_GJ+n2*sin_theta_GJ).Unit()

        n_paral_1 = (n1*math.cos(theta_GJ_max)-n2*math.sin(theta_GJ_max)).Unit()
        n_paral_2 = (n1*math.cos(theta_GJ_max)+n2*math.sin(theta_GJ_max)).Unit()
        new_dir_1 = (n3*P_perp+n_paral_1*P_paral).Unit()
        new_dir_2 = (n3*P_perp+n_paral_2*P_paral).Unit()
        
        condition = new_dir_1.Dot(tau_dir) < new_dir_2.Dot(tau_dir)
        
        new_dir = new_dir_1
        if condition:
            new_dir = new_dir_2
            
        new_P = tau_mom*new_dir
        new_tau.SetXYZT(new_P.X(),new_P.Y(),new_P.Z(),tau.E())
#    print('Initial tau (x,y,z,t) = (%5.3f,%5.3f,%5.3f,%5.3f)'%(tau.X(),tau.Y(),tau.Z(),tau.T()))
#    print('Rotated tau (x,y,z,t) = (%5.3f,%5.3f,%5.3f,%5.3f)'%(new_tau.X(),new_tau.Y(),new_tau.Z(),new_tau.T()))
#    print('')
    return new_tau

###################################################################
# PolVectA1 : routine returns corrected full
# tau 4-momentum and polarimetric vector
# in tau->a1(3-prong)+nu decay channel
# Inputs ->
# PV - TLorentzVector : primary vertex
# SV - TLorentzVector : secondary vertex
# Pt - float : tau pT from fastMTT
# P1 - TLorentzVector : 4-momentum of os charged pion
# P2 - TLorentzVector : 4-momentum of ss1 charged pion
# P3 - TLorentzVector : 4-momentum of ss2 charged pion
# PGen - TLorentzVector : full generator 4-momentum of tau
#                         (only for exploratory studies with MC)
# q  - tau charge
# method - string : options ->
#                   'recoDESY' - DESY-wise computation (preferred)
#                   'recoIC'   - IC-wise computation
#                   'gen'      - generator information is used
# Outputs ->
# tau - TLorentzVector : corrected full 4-momentum of tau
# pv  - TLorentzVector : polarimetric vector
###################################################################
def PolVectA1(PV,SV,
              Pt,P1,P2,P3,PGen,
              q,method):    

#    P1,P2,P3 = sortA1(Pi1,Pi2,Pi3,q1,q2,q3)    
#    print('PolVectA1 : q1 = %3.1f  q2 = %3.1f  q3 = %3.1f   q = %3.1f'%(q1,q2,q3,q))
#    print('PolVectA1 : PV(x) = %6.4f  PV(y) = %6.4f  PV(z) = %6.4f'%(PV.X(),PV.Y(),PV.Z()))
#    print('PolVectA1 : SV(x) = %6.4f  SV(y) = %6.4f  SV(z) = %6.4f'%(SV.X(),SV.Y(),SV.Z()))

    n = SV - PV
    
    visTau = P1+P2+P3
    tau = ROOT.TLorentzVector()
    if method=='gen':
        tau = PGen
    else:
        Ptau = ROOT.TLorentzVector()
        Ptau.SetPtEtaPhiM(Pt,n.Eta(),n.Phi(),utils.tau_mass)
        tau = rotateToGJMax(visTau,Ptau,method)

    tempTau = ROOT.TLorentzVector()
    tempTau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())
    tempP1 = ROOT.TLorentzVector()
    tempP1.SetXYZT(P1.X(),P1.Y(),P1.Z(),P1.T())
    tempP2 = ROOT.TLorentzVector()
    tempP2.SetXYZT(P2.X(),P2.Y(),P2.Z(),P2.T())
    tempP3 = ROOT.TLorentzVector()
    tempP3.SetXYZT(P3.X(),P3.Y(),P3.Z(),P3.T())
    
    boost = -tau.BoostVector()
    if method=='IC':
        tempTau.Boost(boost)
        tempP1.Boost(boost)
        tempP2.Boost(boost)
        tempP3.Boost(boost)
    
    PA1 = PolarimetricA1(tempTau,
                         tempP1,
                         tempP2,
                         tempP3,
                         q)

    xTau = ROOT.TLorentzVector(tempTau.X(),tempTau.Y(),tempTau.Z(),tempTau.T())
    xP1 = ROOT.TLorentzVector(tempP1.X(),tempP1.Y(),tempP1.Z(),tempP1.T())
    xP2 = ROOT.TLorentzVector(tempP2.X(),tempP2.Y(),tempP2.Z(),tempP2.T())
    xP3 = ROOT.TLorentzVector(tempP3.X(),tempP3.Y(),tempP3.Z(),tempP3.T())
    
    PA1_tau = PolarimetricA1(xTau,
                             xP1,
                             xP2,
                             xP3,
                             q)

    tempPV = -PA1.PVC()
    xPV = - PA1_tau.PVC()
#    print('')
#    print('pol. vector from PolarimetricA1.PVC() ->')
#    print('pol. vec. lab frame : (X,Y,Z,T)=(%6.4f,%6.4f,%6.4f,%6.4f)'%(tempPV.X(),tempPV.Y(),tempPV.Z(),tempPV.T()))
#    print('pol. vec. tau frame : (X,Y,Z,T)=(%6.4f,%6.4f,%6.4f,%6.4f)'%(xPV.X(),xPV.Y(),xPV.Z(),xPV.T()))
#    print('')
    
    pv = ROOT.TLorentzVector()
    if method=='IC':
        pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),0.)
#        print('Pol. vector before boost = (%5.3f,%5.3f,%5.3f)'%(pv.X(),pv.Y(),pv.Z()))
        pv.Boost(-boost)
#        print('Pol. vector after boost = (%5.3f,%5.3f,%5.3f)'%(pv.X(),pv.Y(),pv.Z()))
    else:
        pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),tempPV.T())
        
    return tau,pv
    
######################################################################
# acoCP : routine returns acoplanarity angle phi(CP)
# Inputs -> 
# P1_input, P2_input, R1_input, R2_input - TLorentzVectors
# firstNeg - boolean : True if first lepton is negative,
#                      otherwise False
# method_1, method_2 - string : available options
#                               - Impact-Parameter
#                               - Decay-Plane
#                               - Decay-Plane-a1
#                               - PV
# Output -> acop (float): this is tau decay mode specific phi(CP) 
######################################################################
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

            
def CosAlpha(pt,eta,phi,mass,Rx,Ry,Rz):

    P = ROOT.TLorentzVector()
    P.SetPtEtaPhiM(pt,eta,phi,mass)
    ez = ROOT.TVector3(0.,0.,1.)
    p  = ROOT.TVector3(P.Px(),P.Py(),P.Pz())
    n  = ROOT.TVector3(Rx,Ry,Rz)

    ez = ez.Unit();
    p = p.Unit();
    n = n.Unit();

    denom = ez.Cross(p).Mag() * n.Cross(p).Mag()

    cosa = 1.0/math.sqrt(2.0)
    if denom > 0.: 
        cosa = abs( ez.Cross(p).Dot(n.Cross(p)) / denom )
    
    return cosa


