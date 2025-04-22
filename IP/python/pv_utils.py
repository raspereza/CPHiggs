### utils to compute phi(CP)
import ROOT 
import math
from array import array
import numpy as np
import os
import CPHiggs.IP.utils as utils
from CPHiggs.IP.ScaleFactor import ScaleFactor
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1

# Ptau - TLorentzVector : total momentum of tau (from fastMTT)
# Pvis - TLorentzVector : 4p visible tau
# Ppion - TLorentzVector : 4p of charged pion
# IP - TVector3 : impact parameter vector
# PGen - TLorentzVector : generator tau momentum
def TauDirRho(Ptau,Pvis,Ppion,IP,PGen,method):

    n = Pvis.Vect().Unit()
    p = Ppion.Vect().Unit()
    r = IP.Unit()    
    l = p.Cross(r).Unit()
    prod = l.Dot(n)
    if prod<0:
        l = -l

    if method=='GenV1':
        dirGen = PGen.Vect().Unit()
        return dirGen
        
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

        """
        q1 = q1.Unit()
        q2 = q2.Unit()
        px1 = q1.Dot(dirGen)
        px2 = q2.Dot(dirGen)
        if px1>px2:
            q = q1
        else:
            q = q2
        cosOmega1 = q1.Dot(r)        
        cosOmega2 = q2.Dot(r)
        IP_mag = IP.Mag()
        L1 = IP_mag/cosOmega1
        L2 = IP_mag/cosOmega2
#        print('cos(1) = %10.8f    cos(2) = %10.8f : L1 = %10.8f    L2 = %10.8f'%(px1,px2,L1,L2))
#        print('L1 = %10.8f     L2 = %10.8f'%(L1,L2))
        if L1>0 and L2>0:
            beta = utils.tau_mass/(utils.ctau*Ptau.P())
            P1 = 1.0-ROOT.TMath.Exp(-L2*beta)
            P2 = ROOT.TMath.Exp(-L1*beta)
            if P1>P2:
                q = q1
            else:
                q = q2
            s = l.Cross(n).Unit()
            m = n.Cross(s).Unit()
            q = cosThetaGJ*n-sinThetaGJ*m
        else:
            if L1>L2:
                q = q1
            else:
                q = q2
        """                
    dir_tau = q.Unit()
#    cosTJ = dir_tau.Dot(n)
#    perp = dir_tau.Dot(l)
#    print('q*n = %10.8f   q*l = %10.8f'%(cosTJ,perp))
    
    return dir_tau
    
# returns full tau 4-momentum and polarimetric vector
def PolVectRho(Pt,Pi,Pi0,IP,PGen):

    visTau = Pi+Pi0
    tau = ROOT.TLorentzVector()
    GenPt = PGen.Pt()
    tau.SetPtEtaPhiM(Pt,visTau.Eta(),visTau.Phi(),utils.tau_mass)
    vistau_mom = visTau.P()
    tau_mom = tau.P()
    vistau_mass = visTau.M()

#    dir_tau = TauDirRho(tau,visTau,Pi,IP,PGen)
#    mag_tau = dir_tau.Mag()
#    print('Magnitude = %6.4f'%(mag_tau))
    
    tempPi = ROOT.TLorentzVector()
    tempPi.SetXYZT(Pi.X(),Pi.Y(),Pi.Z(),Pi.T())

    tempPi0 = ROOT.TLorentzVector()
    tempPi0.SetXYZT(Pi0.X(),Pi0.Y(),Pi0.Z(),Pi0.T())

    tempTau = ROOT.TLorentzVector()
#    tempTau.SetXYZT(tau_mom*dir_tau.X(),tau_mom*dir_tau.Y(),tau_mom*dir_tau.Z(),tau.E())
    tempTau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())
#    tempTau.SetXYZT(PGen.X(),PGen.Y(),PGen.Z(),PGen.T())

#    boost = -tempTau.BoostVector()
#    tempTau.Boost(boost)
#    tempPi.Boost(boost)
#    tempPi0.Boost(boost)

    tempQ = tempPi - tempPi0
    NN = tempTau - tempPi - tempPi0
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
    
    return tempTau,pv


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

def rotateToGJMax(vistau,tau,method):
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
    if theta_GJ>theta_GJ_max:
        if method=='IC':
            r_vistau = vistau_P.X()/vistau_P.Y()            
            nx = 1.0/math.sqrt(1+r_vistau*r_vistau)
            ny = -nx*r_vistau
            n = ROOT.TVector3(nx,ny,0.)
            m = n.Cross(vistau_dir)
            phi1 = math.atan(tau_dir.Dot(m)/tau_dir.Dot(n))
            phi2 = phi1 + ROOT.TMath.Pi()
            new_dir_1 = math.cos(theta_GJ_max)*vistau_dir+math.sin(theta_GJ_max)*(math.cos(phi1)*n+math.sin(phi1)*m)
            new_dir_2 = math.cos(theta_GJ_max)*vistau_dir+math.sin(theta_GJ_max)*(math.cos(phi2)*n+math.sin(phi2)*m)

            condition = new_dir_1.Dot(tau_dir) > new_dir_2.Dot(tau_dir)
            new_dir = new_dir_1
            if condition:
                new_dir = new_dir_2
            
            new_P = tau_mom*new_dir
            new_tau.SetXYZT(new_P.X(),new_P.Y(),new_P.Z(),tau.E())
            
        else:
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
    print('Method = %s'%(method))
    print('Initial tau (x,y,z,t) = (%5.3f,%5.3f,%5.3f,%5.3f)'%(tau.X(),tau.Y(),tau.Z(),tau.T()))
    print('Rotated tau (x,y,z,t) = (%5.3f,%5.3f,%5.3f,%5.3f)'%(new_tau.X(),new_tau.Y(),new_tau.Z(),new_tau.T()))
    print('')
    return new_tau


# returns corrected full tau 4-momentum and polarimetric vector in 
def PolVectA1(PV,SV,
              Pt,P1,P2,P3,
              q,method):    

#    P1,P2,P3 = sortA1(Pi1,Pi2,Pi3,q1,q2,q3)    
#    print('PolVectA1 : q1 = %3.1f  q2 = %3.1f  q3 = %3.1f   q = %3.1f'%(q1,q2,q3,q))
#    print('PolVectA1 : PV(x) = %6.4f  PV(y) = %6.4f  PV(z) = %6.4f'%(PV.X(),PV.Y(),PV.Z()))
#    print('PolVectA1 : SV(x) = %6.4f  SV(y) = %6.4f  SV(z) = %6.4f'%(SV.X(),SV.Y(),SV.Z()))

    n = SV - PV
    
    Ptau = ROOT.TLorentzVector()
    Ptau.SetPtEtaPhiM(Pt,n.Eta(),n.Phi(),utils.tau_mass)
    visTau = P1+P2+P3
#    tau = rotateToGJ(visTau,Ptau)
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

    tempPV = -PA1.PVC()
    pv = ROOT.TLorentzVector()
    if method=='IC':
        pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),0.)
        print('Pol. vector before boost = (%5.3f,%5.3f,%5.3f)'%(pv.X(),pv.Y(),pv.Z()))
        pv.Boost(-boost)
        print('Pol. vector after boost = (%5.3f,%5.3f,%5.3f)'%(pv.X(),pv.Y(),pv.Z()))
    else:
        pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),tempPV.T())
        
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


