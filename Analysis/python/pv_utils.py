### utilities to compute phi(CP)
import ROOT
import math
from array import array
import numpy as np
import os
import CPHiggs.Analysis.utils as utils
from CPHiggs.PolarimetricVector.PolarimetricA1 import PolarimetricA1
from KinFit import kinfit_3pr, kinfit_3pr_3pr  

def FastMTT_3pr_cpp(P1,P2,metPt,metPhi,metcovXX,metcovXY,metcovYY,
                    sv,
                    svcovXX,svcovXY,svcovYY,
                    svcovXZ,svcovYZ,svcovZZ,
                    phiScan,
                    massX):

    Pt_1 = np.zeros(1,dtype=np.float64)
    Eta_1 = np.zeros(1,dtype=np.float64)
    Phi_1 = np.zeros(1,dtype=np.float64)
    Mass_1 = np.zeros(1,dtype=np.float64)
    
    Pt_2 = np.zeros(1,dtype=np.float64)
    Eta_2 = np.zeros(1,dtype=np.float64)
    Phi_2 = np.zeros(1,dtype=np.float64)
    Mass_2 = np.zeros(1,dtype=np.float64)

    svX = np.zeros(1,dtype=np.float64)
    svY = np.zeros(1,dtype=np.float64)
    svZ = np.zeros(1,dtype=np.float64)
    
    Pt_1[0] = P1.Pt()
    Eta_1[0] = P1.Eta()
    Phi_1[0] = P1.Phi()
    Mass_1[0] = P1.M()
    
    Pt_2[0] = P2.Pt()
    Eta_2[0] = P2.Eta()
    Phi_2[0] = P2.Phi()
    Mass_2[0] = P2.M()

    svX[0] = sv.X()
    svY[0] = sv.Y()
    svZ[0] = sv.Z()
    
    results = kinfit_3pr(Pt_1,Eta_1,Phi_1,Mass_1,
                         Pt_2,Eta_2,Phi_2,Mass_2,
                         metPt,metPhi,metcovXX,metcovXY,metcovYY,
                         svX,svY,svZ,
                         svcovXX,svcovXY,svcovYY,svcovXZ,svcovYZ,svcovZZ,
                         phiScan,
                         massX)
    
    Px1 = results['px_1'][0]
    Py1 = results['py_1'][0]
    Pz1 = results['pz_1'][0]
    E1 = math.sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1+utils.tau_mass*utils.tau_mass)
    
    Px2 = results['px_2'][0]
    Py2 = results['py_2'][0]
    Pz2 = results['pz_2'][0]
    E2 = math.sqrt(Px2*Px2+Py2*Py2+Pz2*Pz2+utils.tau_mass*utils.tau_mass)

    P1 = ROOT.TLorentzVector(Px1,Py1,Pz1,E1)
    P2 = ROOT.TLorentzVector(Px2,Py2,Pz2,E2)
    
    chi2 = results['chi2'][0]
    #    print('Mass = %5.1f : %5.1f chi2 = %5.1f'%(massX,(P1+P2).M(),chi2))
    return P1,P2,chi2


def FastMTT_3pr_3pr_cpp(P1,P2,metPt,metPhi,
                        metcovXX,metcovXY,metcovYY,
                        sv1,
                        sv1covXX,sv1covXY,sv1covYY,
                        sv1covXZ,sv1covYZ,sv1covZZ,
                        sv2,
                        sv2covXX,sv2covXY,sv2covYY,
                    	sv2covXZ,sv2covYZ,sv2covZZ,
                        massX):

    Pt_1 = np.zeros(1,dtype=np.float64)
    Eta_1 = np.zeros(1,dtype=np.float64)
    Phi_1 = np.zeros(1,dtype=np.float64)
    Mass_1 = np.zeros(1,dtype=np.float64)
    
    Pt_2 = np.zeros(1,dtype=np.float64)
    Eta_2 = np.zeros(1,dtype=np.float64)
    Phi_2 = np.zeros(1,dtype=np.float64)
    Mass_2 = np.zeros(1,dtype=np.float64)

    sv1X = np.zeros(1,dtype=np.float64)
    sv1Y = np.zeros(1,dtype=np.float64)
    sv1Z = np.zeros(1,dtype=np.float64)
    
    sv2X = np.zeros(1,dtype=np.float64)
    sv2Y = np.zeros(1,dtype=np.float64)
    sv2Z = np.zeros(1,dtype=np.float64)
    
    Pt_1[0] = P1.Pt()
    Eta_1[0] = P1.Eta()
    Phi_1[0] = P1.Phi()
    Mass_1[0] = P1.M()
    
    Pt_2[0] = P2.Pt()
    Eta_2[0] = P2.Eta()
    Phi_2[0] = P2.Phi()
    Mass_2[0] = P2.M()

    sv1X[0] = sv1.X()
    sv1Y[0] = sv1.Y()
    sv1Z[0] = sv1.Z()
    
    sv2X[0] = sv2.X()
    sv2Y[0] = sv2.Y()
    sv2Z[0] = sv2.Z()
    
    verbosity = False
    results = kinfit_3pr_3pr(Pt_1,Eta_1,Phi_1,Mass_1,
                             Pt_2,Eta_2,Phi_2,Mass_2,
                             metPt,metPhi,metcovXX,metcovXY,metcovYY,
                             sv1X,sv1Y,sv1Z,
                             sv1covXX,sv1covXY,sv1covYY,sv1covXZ,sv1covYZ,sv1covZZ,
                             sv2X,sv2Y,sv2Z,
                             sv2covXX,sv2covXY,sv2covYY,sv2covXZ,sv2covYZ,sv2covZZ,
                             verbosity,
                             massX)
    
    Px1 = results['px_1'][0]
    Py1 = results['py_1'][0]
    Pz1 = results['pz_1'][0]
    E1 = math.sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1+utils.tau_mass*utils.tau_mass)
    
    Px2 = results['px_2'][0]
    Py2 = results['py_2'][0]
    Pz2 = results['pz_2'][0]
    E2 = math.sqrt(Px2*Px2+Py2*Py2+Pz2*Pz2+utils.tau_mass*utils.tau_mass)

    P1 = ROOT.TLorentzVector(Px1,Py1,Pz1,E1)
    P2 = ROOT.TLorentzVector(Px2,Py2,Pz2,E2)
    
    chi2 = results['chi2'][0]

    return P1,P2,chi2


def tauMomGJ(E_vis,m_vis,cos_theta_GJ):

    P_vis = math.sqrt(E_vis*E_vis-m_vis*m_vis)
    
    A = 4 * (E_vis**2-P_vis*P_vis*cos_theta_GJ*cos_theta_GJ)
    B = 4 * P_vis*(utils.tau_mass*utils.tau_mass+m_vis*m_vis)*cos_theta_GJ
    C = 4*utils.tau_mass*utils.tau_mass*E_vis*E_vis - (utils.tau_mass*utils.tau_mass+m_vis*m_vis)**2
    
#    B = 4 * E_vis*(utils.tau_mass*utils.tau_mass+m_vis*m_vis)
#    C = (utils.tau_mass*utils.tau_mass+m_vis*m_vis)**2 + 4*utils.tau_mass*utils.tau_mass*cos_theta_GJ*cos_theta_GJ*P_vis*P_vis

    discriminant = B**2 - 4*A*C
    solutions = []
    if discriminant < 0:
        return solutions
    
    sqrt_disc = math.sqrt(discriminant)
    
    for sign in [+1, -1]:
        p_tau = (B + sign * sqrt_disc) / (2 * A)
        if p_tau > 0:
            solutions.append(p_tau)
            
    return solutions

def Basis(P,SV):
    nz = P.Vect().Unit()
    sv = SV.Unit()
    ny = nz.Cross(sv).Unit()
    nx = ny.Cross(nz)
    return nx,ny,nz

def ThetaGJ(P,Pvis,mass):
    E = math.sqrt(P*P+utils.tau_mass*utils.tau_mass)
    Evis = math.sqrt(Pvis*Pvis+mass*mass)
    num = 2*E*Evis - mass*mass - utils.tau_mass*utils.tau_mass
    den = 2*P*Pvis
    cosTheta = max(-1.0,min(1.0,num/den))
    theta = math.acos(cosTheta)
    return theta
    
# Computes angle between visible
# momentum of tau and secondary vertex
# also computes maximal Gottfried-Jackson
# angle.
def Theta(P,SV):
    vistau_dir = P.Vect().Unit()
    mass_vis_tau = P.M()
    vistau_mom = P.P()
    tau_dir = SV.Unit()
    cos_theta = max(-1.0,min(1.0,vistau_dir.Dot(tau_dir)))
    theta = math.acos(cos_theta)
    sin_theta_GJ_max = max(-1.0,min(1.0,0.5*(utils.tau_mass*utils.tau_mass-mass_vis_tau*mass_vis_tau)/(utils.tau_mass*vistau_mom)))
    theta_GJ_max = math.asin(sin_theta_GJ_max)
    return theta,theta_GJ_max


# Kinematic fit with the vertex constraint
# Decay tau->X + tau->3prong
def FastMTT_3pr(P1, P2, # TLorentzVectors of P1 and P2 (visible)
                MET, # Vector of MET 
                met_covxx,met_covxy,met_covyy,
                SV2,
                sv2_covxx,sv2_covyx,sv2_covyy,
                sv2_covzx,sv2_covzy,sv2_covzz):

    x2_min = P2.M()*P2.M()/(utils.tau_mass*utils.tau_mass)
    x1_min = P1.M()*P1.M()/(utils.tau_mass*utils.tau_mass)

    Evis = P2.E()
    Pvis = P2.P()
    mvis = P2.M()
    
    mass_vis = (P1+P2).M()
    mass_ratio = mass_vis*mass_vis/(utils.massH*utils.massH)

    theta_sv2, theta_GJ_max = Theta(P2,SV2)
    Nx,Ny,Nz = Basis(P2,SV2)

    metcov_det = met_covxx[0]*met_covyy[0] - met_covxy[0]*met_covxy[0]
    illDefinedMetCov = False
    if abs(metcov_det)<1e-10:
        metcov_det = met_covxx[0]*met_covyy[0]
        invalidMetCov = True

    chi2_min = 1.0e+12
    x1_opt = 0.5
    x2_opt = 0.5
    theta_opt = min(theta_sv2,theta_GJ_max)
    sv2_opt = math.sin(theta_opt)*Nx+math.cos(theta_opt)*Nz
    
    for i in range(1,100):
        x2 = float(i)*0.01
        if x2>1.0 or x2<x2_min: continue
        p2 = P2.P()/x2
        theta = ThetaGJ(p2,Pvis,mvis)
        if theta>theta_GJ_max: continue
        x1 = mass_ratio/x2
        if x1>1.0 or x1<x1_min: continue
        p1 = P1.P()/x1

        nuX = P1.X()*(1.0/x1-1.0) + P2.X()*(1.0/x2-1.0)
        nuY = P1.Y()*(1.0/x1-1.0) + P2.Y()*(1.0/x2-1.0)
        ResidualMet_X = nuX - MET.X()
        ResidualMet_Y = nuY - MET.Y()
        chi2_met = 0.
        if illDefinedMetCov:
            chi2_met += residualX*residualX/met_covxx[0]
            chi2_met += residualY*residualY/met_covyy[0]
        else:
            chi2_met += ResidualMet_X*(met_covxx[0]*ResidualMet_X+met_covxy[0]*ResidualMet_Y)
            chi2_met += ResidualMet_Y*(met_covxy[0]*ResidualMet_X+met_covyy[0]*ResidualMet_Y)
        chi2_met /= metcov_det

        for sign in [-1,1]:
            N_unit = sign*math.sin(theta)*Nx+math.cos(theta)*Nz
            Residual_N = SV2.Mag()*SV2.Mag()*(sv2_opt.Unit() - N_unit)
            chi2_sv2 = Residual_N.X()*Residual_N.X()/sv2_covxx[0]
            chi2_sv2 += Residual_N.Y()*Residual_N.Y()/sv2_covyy[0]
            chi2_sv2 += Residual_N.Z()*Residual_N.Z()/sv2_covzz[0]        
            chi2_tot = chi2_sv2 + chi2_met
            if chi2_tot<chi2_min:
                chi2_min = chi2_tot
                x1_opt = x1
                x2_opt = x2
                theta_opt = theta
            
    P1_mom = P1.P()/x1_opt
    P1_x   = P1.X()/x1_opt
    P1_y   = P1.Y()/x1_opt
    P1_z   = P1.Z()/x1_opt
    P1_en  = math.sqrt(P1_mom*P1_mom+utils.tau_mass*utils.tau_mass)
    P1_out = ROOT.TLorentzVector(P1_x,P1_y,P1_z,P1_en)
    
    P2_mom  = P2.P()/x2_opt
    P2_unit = Nx*math.sin(theta_opt)+Nz*math.cos(theta_opt)
    P2_x    = P2_mom*P2_unit.X()
    P2_y    = P2_mom*P2_unit.Y()
    P2_z    = P2_mom*P2_unit.Z()
    P2_en   = math.sqrt(P2_mom*P2_mom+utils.tau_mass*utils.tau_mass)
    P2_out  = ROOT.TLorentzVector(P2_x,P2_y,P2_z,P2_en)
    
    return P1_out,P2_out,chi2_min

# Kinematic fit with the vertex constraint
def FastMTT_3pr_3pr(P1, P2,
                    MET,  
                    met_covxx,met_covxy,met_covyy,
                    SV1,
                    sv1_covxx,sv1_covyx,sv1_covyy,
                    sv1_covzx,sv1_covzy,sv1_covzz,                
                    SV2,
                    sv2_covxx,sv2_covyx,sv2_covyy,
                    sv2_covzx,sv2_covzy,sv2_covzz):

    Evis1 = P1.E()
    Pvis1 = P1.P()    
    mvis1 = P1.M()
    
    Evis2 = P2.E()
    Pvis2 = P2.P()
    mvis2 = P2.M()

    x1_min = mvis1*mvis1/(utils.tau_mass*utils.tau_mass)
    x2_min = mvis2*mvis2/(utils.tau_mass*utils.tau_mass)
    
    mass_vis = (P1+P2).M()
    mass_ratio = mass_vis*mass_vis/(utils.massH*utils.massH)

    metcov_det = met_covxx[0]*met_covyy[0] - met_covxy[0]*met_covxy[0]
    illDefinedMetCov = False
    if abs(metcov_det)<1e-10:
        illDefinedMetCov = True
        
    #
    # angles
    #
    theta_sv1, theta_GJ_max1 = Theta(P1,SV1)
    theta_sv2, theta_GJ_max2 = Theta(P2,SV2)

    ###########################
    Nx1,Ny1,Nz1 = Basis(P1,SV1)
    Nx2,Ny2,Nz2 = Basis(P2,SV2)
    ###########################
    
    chi2_min = 1.0e+12
    theta1_opt = min(theta_sv1,theta_GJ_max1)
    theta2_opt = min(theta_sv2,theta_GJ_max2)
    x1_opt = 0.5
    x2_opt = 0.5

    sv1_opt = math.sin(theta1_opt)*Nx1+math.cos(theta1_opt)*Nz1
    sv2_opt = math.sin(theta2_opt)*Nx2+math.cos(theta2_opt)*Nz2
    
    # scan over x2
    for i in range(1,100):
        x2 = float(i)*0.01
        if x2>1.0 or x2<x2_min: continue
        x1 = mass_ratio/x2
        if x1>1.0 or x1<x2_min: continue
        
        p2 = P2.P()/x2
        theta2 = ThetaGJ(p2,Pvis2,mvis2)
        if theta2>theta_GJ_max2:continue

        p1 = P1.P()/x1
        theta1 = ThetaGJ(p1,Pvis1,mvis1)
        if theta1>theta_GJ_max1: continue

        nuX = P1.X()*(1.0/x1-1.0) + P2.X()*(1.0/x2-1.0)
        nuY = P1.Y()*(1.0/x1-1.0) + P2.Y()*(1.0/x2-1.0)
        ResidualMet_X = nuX - MET.X()
        ResidualMet_Y = nuY - MET.Y()

        chi2_met = 0.
        if illDefinedMetCov:
            chi2_met += residualX*residualX/met_covxx[0]
            chi2_met += residualY*residualY/met_covyy[0]
        else:
            chi2_met += ResidualMet_X*(met_covxx[0]*ResidualMet_X+met_covxy[0]*ResidualMet_Y)
            chi2_met += ResidualMet_Y*(met_covxy[0]*ResidualMet_X+met_covyy[0]*ResidualMet_Y)
        chi2_met /= metcov_det

        for sign1 in [-1,1]:
            for sign2 in [-1,1]:        
                N1_unit = sign1*math.sin(theta1)*Nx1+math.cos(theta1)*Nz1
                N2_unit = sign2*math.sin(theta2)*Nx2+math.cos(theta2)*Nz2
                Residual_N1 = SV1.Mag()*SV1.Mag()*(sv1_opt.Unit() - N1_unit)
                Residual_N2 = SV2.Mag()*SV2.Mag()*(sv2_opt.Unit() - N2_unit)

                chi2_sv1 = Residual_N1.X()*Residual_N1.X()/sv1_covxx[0]
                chi2_sv1 += Residual_N1.Y()*Residual_N1.Y()/sv1_covyy[0]
                chi2_sv1 += Residual_N1.Z()*Residual_N1.Z()/sv1_covzz[0]

                chi2_sv2 = Residual_N2.X()*Residual_N2.X()/sv2_covxx[0]
                chi2_sv2 += Residual_N2.Y()*Residual_N2.Y()/sv2_covyy[0]
                chi2_sv2 += Residual_N2.Z()*Residual_N2.Z()/sv2_covzz[0]
                
                chi2_tot = chi2_sv1 + chi2_sv2 + chi2_met
        
                if chi2_tot<chi2_min:
                    chi2_min = chi2_tot
                    x1_opt = x1
                    x2_opt = x2
                    theta1_opt = sign1*theta1
                    theta2_opt = sign2*theta2

    P1_mom = P1.P()/x1_opt
    P1_unit = Nx1*math.sin(theta1_opt)+Nz1*math.cos(theta1_opt)
    P1_x = P1_mom*P1_unit.X()
    P1_y = P1_mom*P1_unit.Y()
    P1_z = P1_mom*P1_unit.Z()
    P1_en = math.sqrt(P1_mom*P1_mom+utils.tau_mass*utils.tau_mass)    
    P1_out = ROOT.TLorentzVector(P1_x,P1_y,P1_z,P1_en)
    
    P2_mom = P2.P()/x2_opt
    P2_unit = Nx2*math.sin(theta2_opt)+Nz2*math.cos(theta2_opt)
    P2_x = P2_mom*P2_unit.X()
    P2_y = P2_mom*P2_unit.Y()
    P2_z = P2_mom*P2_unit.Z()
    P2_en = math.sqrt(P2_mom*P2_mom+utils.tau_mass*utils.tau_mass)
    P2_out = ROOT.TLorentzVector(P2_x,P2_y,P2_z,P2_en)
    
    return P1_out, P2_out, chi2_min

################################################################
# PolVectRho : routine returns full tau 4-momentum and
# polarimetric vector in tau->rho+nu decay mode
# Inputs ->
# P   - TLorentzVector : 4-momentum of tau
# Pi  - TLorentzVector : 4-momentum of charged pion
# Pi0 - TLorentzVector : 4-momentum of neutral pion
# pv  - TLorentzVector : polarimetric vector
#################################################################
def PolVectRho(P,Pi,Pi0):

    visTau = Pi+Pi0

    tauP4 = ROOT.TLorentzVector(P.X(),P.Y(),P.Z(),P.T())
        
    tempPi = ROOT.TLorentzVector()
    tempPi.SetXYZT(Pi.X(),Pi.Y(),Pi.Z(),Pi.T())

    tempPi0 = ROOT.TLorentzVector()
    tempPi0.SetXYZT(Pi0.X(),Pi0.Y(),Pi0.Z(),Pi0.T())

    tempQ = tempPi - tempPi0
    NN = tauP4 - tempPi - tempPi0
    tempN = ROOT.TLorentzVector()
    tempN.SetPtEtaPhiM(NN.Pt(),NN.Eta(),NN.Phi(),0.)

    X1 = tempQ.E()*tempN.E()-tempQ.Px()*tempN.Px()-tempQ.Py()*tempN.Py()-tempQ.Pz()*tempN.Pz()
    X2 = tempQ.E()*tempQ.E()-tempQ.Px()*tempQ.Px()-tempQ.Py()*tempQ.Py()-tempQ.Pz()*tempQ.Pz()
    tempPV = X1*tempQ-X2*tempN

    pv = ROOT.TLorentzVector()
    pv.SetXYZT(tempPV.X(),tempPV.Y(),tempPV.Z(),tempPV.T())
    
    return pv


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

#    print('')
#    print('a1 sorted -> ')
#    print('pi_os  (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p1.Pt(),p1.Eta(),p1.Phi(),p1.M()))
#    print('pi_ss1 (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p2.Pt(),p2.Eta(),p2.Phi(),p2.M()))
#    print('pi_ss2 (pt,eta,phi,mass) = (%5.1f,%5.3f,%5.3f,%5.3f)'%(p3.Pt(),p3.Eta(),p3.Phi(),p3.M()))
#    print('M(os_ss1) = %5.4f -- M(os_ss2) = %5.4f -- m(rho) = %5.4f'%((p1+p2).M(),(p1+p3).M(),utils.rho_mass))

    
    return p1,p2,p3

###########################################################
# Rotation of visible tau momentum towards reconstructed  #
# direction (e.g. SV-PV) to match theta_GJmax             #
# vistau - 4-momentum of visible tau decay products       #
# tau - full 4-momentum of tau                            #
###########################################################
def rotateToGJ(Pvis,Pt,sv):
    
    tau = ROOT.TLorentzVector()
    SV = ROOT.TLorentzVector(sv.X(),sv.Y(),sv.Z(),0.)
    tau.SetPtEtaPhiM(Pt,SV.Eta(),SV.Phi(),utils.tau_mass)

    vistau_P = Pvis.Vect()
    tau_P = tau.Vect()
    vistau_dir = vistau_P.Unit()
    tau_dir = tau_P.Unit()
    vistau_mom = Pvis.P()
    tau_mom = tau.P()
    tau_E = tau.E()
    vistau_E = tau.E()
    mass_vis_tau = Pvis.M()

    cos_GJ = max(-1.0,min(1.0,vistau_dir.Dot(tau_dir)))
    theta_GJ = math.acos(cos_GJ) 
    sin_GJ_max = 0.5*(utils.tau_mass*utils.tau_mass-mass_vis_tau*mass_vis_tau)/(utils.tau_mass*vistau_mom)
    theta_GJ_max = math.asin(max(-1.0,min(1.0,sin_GJ_max)))

    new_tau = ROOT.TLorentzVector()
    new_tau.SetXYZT(tau.X(),tau.Y(),tau.Z(),tau.T())
    #    print('')
    #    print('Routine rotateToGJMax ->')
    if theta_GJ>theta_GJ_max:
        #        print('theta_GJ (%5.3f) > theta_GJ_max (%5.3f) -> rotation of tau 3-momentum'%(theta_GJ,theta_GJ_max))
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
    return new_tau,theta_GJ,theta_GJ_max

def PolVectA1(P,P1,P2,P3,q):
    PA1_tau = PolarimetricA1(P,
                             P1,
                             P2,
                             P3,
                             q)
    PV = - PA1_tau.PVC()
    return PV


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

    #    print('After boost into ZMF of visible products')
    #    print('P1 (x,y,z)=(%5.3f,%5.3f,%5.3f)'%(vecP1.X(),vecP1.Y(),vecP1.Z()))
    #    print('P2 (x,y,z)=(%5.3f,%5.3f,%5.3f)'%(vecP2.X(),vecP2.Y(),vecP2.Z()))
    #    print('R1 (x,y,z)=(%5.3f,%5.3f,%5.3f)'%(vecR1.X(),vecR1.Y(),vecR1.Z()))
    #    print('R2 (x,y,z)=(%5.3f,%5.3f,%5.3f)'%(vecP1.X(),vecR2.Y(),vecR2.Z()))
    
    
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


    if math.isnan(acop): acop = -9999.
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


def AnalysisGenerator(PtFastMTT,Ptau,Pvis,PGen,sv):
    
    variables = {}


    dPt = Ptau.Pt()/PGen.Pt()
    dPt_FastMTT = PtFastMTT/PGen.Pt()
    variables['dPt_KinFit'] = dPt
    variables['dPt'] = dPt_FastMTT
    
    P_unit = Ptau.Vect().Unit()
    Pvis_unit = Pvis.Vect().Unit()
    PGen_unit = PGen.Vect().Unit()
    PSV,thetaSV,thetaMax = rotateToGJ(Pvis,PtFastMTT,sv)
    PSV_unit = PSV.Vect().Unit()
    dAlpha = math.acos(max(-1.,min(1.,P_unit.Dot(PGen_unit))))
    dAlphaSV = math.acos(max(-1.,min(1.,PSV_unit.Dot(PGen_unit))))
    dAlphaKinFitSV = math.acos(max(-1.,min(1.,PSV_unit.Dot(P_unit))))
    variables['dAlpha'] = dAlphaSV
    variables['dAlpha_KinFit'] = dAlpha
    variables['dAlpha_KinFitSV'] = dAlphaKinFitSV
    
    thetaGJ_gen = ThetaGJ(PGen.P(),Pvis.P(),Pvis.M())
    cosThetaSV = max(-1.,min(1.,Pvis_unit.Dot(PSV_unit)))
    thetaSV = math.acos(cosThetaSV)
    cosThetaGen = max(-1.,min(1.,Pvis_unit.Dot(PGen_unit)))
    thetaGen = math.acos(cosThetaGen)
    cosTheta = max(-1.,min(1.,Pvis_unit.Dot(P_unit)))
    theta = math.acos(cosTheta)
    variables['dTheta'] = thetaSV-thetaGJ_gen
    variables['dTheta_KinFit'] = theta-thetaGJ_gen
    variables['dTheta_KinFitSV'] = theta-thetaSV
    
    return variables

def GeneratorVector(PGen,Pvis,sv):

    #    p_unit = PGen.Vect().Unit()
    #    pvis_unit = Pvis.Vect().Unit()
    theta = ThetaGJ(PGen.P(),Pvis.P(),Pvis.M())
    #    cosTheta = max(-1.,min(1.,p_unit.Dot(pvis_unit)))
    #    theta = math.acos(cosTheta)
    Nx,Ny,Nz = Basis(Pvis,sv)
    Direction = math.sin(theta)*Nx+math.cos(theta)*Nz
    N_unit = Direction.Unit()
    P_x = PGen.P()*N_unit.X()
    P_y = PGen.P()*N_unit.Y()
    P_z = PGen.P()*N_unit.Z()
    P_en = PGen.T()
    
    Pout = ROOT.TLorentzVector(P_x,P_y,P_z,P_en)
    return Pout
