#####################################################################################
# Based on original version of FastMTT algorithm                                    #
# https://github.com/SVfit/ClassicSVfit/blob/fastMTT_19_02_2019                     #
# http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_032_v3.pdf         #
# Author : Artur Kalinowski                                                         #
#                                                                                   #
# Modified version includes soft mass constraint (mZ or mH).                        #
# Only those (x1,x2) solutions are considered which yield di-tau                    #
# mass within certain mass window: Mass - Width < m(tautau) < Mass + Width          #
# Default values: Mass = 125 GeV, Width = 4 GeV -> window = [123,127] GeV           #
# was applied in the A->Zh(125)->(ll)(tautau) search, arXiv 2501.14825 (HIG-22-004) #
# Authors : Gage DeZoort and Alexei Raspereza                                       #
#                                                                                   #
# Addition (Alexei Raspereza, 6.07.2025) : implemented hard mass constraint         #
# Scan is performed over variable x2. At each scanned values of x2, x1 is fixed     #
# by mass constraint x1 = mvis^2/(x2*mX), where mX is the mass of the resonance     # 
#####################################################################################

import os
import sys
import math
import vector as vec
import numpy as np
import numba as nb
import awkward as ak
import ROOT
from numpy.linalg import inv

# Convolution of TVector3 V and TMatrix M:
# V^T * M * V
def VectorMatrixConv(V,M):
    # vector V^T
    XT = V.X()
    YT = V.Y()
    ZT = V.Z()
    # vector M * V
    X = M(0,0)*V.X()+M(0,1)*V.Y()+M(0,2)*V.Z()
    Y = M(1,0)*V.X()+M(1,1)*V.Y()+M(1,2)*V.Z()
    Z = M(2,0)*V.X()+M(2,1)*V.Y()+M(2,2)*V.Z()

    result = XT*X + YT*Y + ZT*Z
    return result
    
# tautau -> X+a1 decay: scan over x2 and phi to compute best
# estimate of 
@nb.jit(nopython=True, parallel=False)
def fastmtt_a1(pt_1, eta_1, phi_1, mass_1, decay_type_1,
               pt_2, eta_2, phi_2, mass_2, decay_type_2,
               met_x, met_y, metcov_xx, metcov_xy, metcov_yx, metcov_yy,
               PV_x, PV_y, PV_z, SV_x, SV_y, SV_z, # primary and secondary vertexes
               svcov_xx, svcov_xy, svcov_xz, svcov_yy, svcov_yz, svcov_zz, # cov. of SV
               verbosity=-1, delta=1/1.15, reg_order=6,
               Mass=125.0,
               Width=2.0,
               npoints_x2=50.,
               npoints_phi=50.):
    
    # initialize global parameters
    step_x2 = 1.0/npoints_x2
    twopi = 2*math.pi
    step_phi = twopi/npoints_phi
    m_ele = 0.51100e-3
    m_muon = 0.10566
    m_tau = 1.77685
    m_pion = 0.13957
    mass_dict = {0: m_ele, 1: m_muon}

    mass_min = Mass - Width
    mass_max = Mass + Width
    
    # initialize higgs->ditau decays, tau decay types
    N = len(pt_1)
    px_out = np.zeros(N, dtype=np.float32)
    py_out = np.zeros(N, dtype=np.float32)
    pz_out = np.zeros(N, dtype=np.float32)
    
    # loop over all events, calculate corrected ditau mass
    for i in range(N):
        
        # grab the correct masses based on tau decay type
        # tau decay_type: 0 ==> leptonic to electron, 
        #                 1 ==> leptonic to muon, 
        #                 2 ==> leptonic to hadronic
        if (decay_type_1[i] != 2): m1 = mass_dict[decay_type_1[i]]
        else: m1 = mass1[i]
        if (decay_type_2[i] != 2): m2 = mass_dict[decay_type_2[i]]
        else: m2 = mass2[i]
            
        # store visible masses
        m_vis_1 = m1
        m_vis_2 = m2
        
        # determine minimum and maximum possible masses
        m_vis_min_1, m_vis_max_1 = 0, 0
        m_vis_min_2, m_vis_max_2 = 0, 0
        if (decay_type_1[i] == 0): m_vis_min_1, m_vis_max_1 = m_ele, m_ele
        if (decay_type_1[i] == 1): m_vis_min_1, m_vis_max_1 = m_muon, m_muon
        if (decay_type_1[i] == 2): m_vis_min_1, m_vis_max_1 = m_pion, 1.5
        if (decay_type_2[i] == 0): m_vis_min_2, m_vis_max_2 = m_ele, m_ele
        if (decay_type_2[i] == 1): m_vis_min_2, m_vis_max_2 = m_muon, m_muon
        if (decay_type_2[i] == 2): m_vis_min_2, m_vis_max_2 = m_pion, 1.5
        if (m_vis_1 < m_vis_min_1): m_vis_1 = m_vis_min_1
        if (m_vis_1 > m_vis_max_1): m_vis_1 = m_vis_max_1
        if (m_vis_2 < m_vis_min_2): m_vis_2 = m_vis_min_2
        if (m_vis_2 > m_vis_max_2): m_vis_2 = m_vis_max_2
          
        # store both tau candidate four vectors
        leg1 = vec.obj(pt=pt_1[i], eta=eta_1[i], phi=phi_1[i], mass=m_vis_1)
        leg2 = vec.obj(pt=pt_2[i], eta=eta_2[i], phi=phi_2[i], mass=m_vis_2)

        # store visible mass of ditau pair
        m_vis = math.sqrt(2*leg1.pt*leg2.pt*(math.cosh(leg1.eta - leg2.eta) - 
                                             math.cos(leg1.phi - leg2.phi)))
        m_tt_vis[i] = m_vis
            
        # correct initial visible masses
        if (decay_type_1[i] == 2 and m_vis_1 > 1.5): m_vis_1 = 0.3
        if (decay_type_2[i] == 2 and m_vis_2 > 1.5): m_vis_2 = 0.3

        # invert met covariance matrix, calculate determinant
        metcovinv_xx, metcovinv_yy = metcov_yy[i], metcov_xx[i]
        metcovinv_xy, metcovinv_yx = -metcov_xy[i], -metcov_yx[i]
        metcovinv_det = (metcovinv_xx*metcovinv_yy -
                         metcovinv_yx*metcovinv_xy)
        if (metcovinv_det<1e-10): 
                print("Warning! Ill-conditioned MET covariance at event index", i)
                continue

        # define estimate of tau direction as unit vector
        # in the direction pointing from PV to SV
        # and define inverted matrix for the direction
        pv = ROOT.TVector3(PV_x[i],PV_y[i],PV_z[i])
        sv = ROOT.TVector3(SV_x[i],SV_y[i],PV_z[i])
        ntau = sv - pv
        magnitute_ntau = ntau.Mag()
        magnitude_ntau2 = magnitute_ntau*magnitute_ntau
        ntau /= magnitude_ntau
        ntau_cov_xx = svcov_xx[i]/magnitude_ntau2
        ntau_cov_xy = svcov_xy[i]/magnitude_ntau2
        ntau_cov_xz = svcov_xz[i]/magnitude_ntau2
        ntau_cov_yy = svcov_yy[i]/magnitude_ntau2
        ntau_cov_yz = svcov_yz[i]/magnitude_ntau2
        ntau_cov_zz = svcov_zz[i]/magnitude_ntau2
        ntau_cov = ROOT.TMatrixD()
        ntau_cov(0,0) = ntau_cov_xx
        ntau_cov(0,1) = ntau_cov_xy
        ntau_cov(0,2) = ntau_cov_xz
        ntau_cov(1,0) = ntau_cov_xy
        ntau_cov(1,1) = ntau_cov_yz
        ntau_cov(1,2) = ntau_cov_yz
        ntau_cov(2,0) = ntau_cov_xz
        ntau_cov(2,1) = ntau_cov_yz
        ntau_cov(2,2) = ntau_cov_zz

        ntau_inv = ntau_cov.Invert()

        # Define local coordinate system
        # with z' ort oriented along visible
        # tau2 momentum
        Pvis = ROOT.TVector3(leg2.px,leg2.py,leg2.pz)
        zprime = PVis.Unit()
        # beam axis :
        zaxis = ROOT.TVector(0.,0.,1.)
        # define orts y' and x'
        yprime = zprime.Cross(zaxis).Unit()
        xprime = yprime.Cross(zprime).Unit()
        
        # perform likelihood scan 
        # see http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_032_v3.pdf
        met_const = 1/(2*math.pi*math.sqrt(metcovinv_det))
        # some initializations
        min_likelihood = 1e9
        PX = (leg1+leg2).px
        PY = (leg1+leg2).py
        PZ = (leg1+leg2).pz
        x1_opt, x2_opt = step_x2, step_x2
        px_opt, py_opt, pz_opt = PX, PY, PZ 
        
        # scan over x2 variable
        for x2 in np.arange(step_x2, 1.00, step_x2):
            # hard mass constraintm
            x1 = m_vis*m_vis/(Mass*Mass*x2)
            x1_min = min(1, math.pow((m_vis_1/m_tau), 2))
            x2_min = min(1, math.pow((m_vis_2/m_tau),2))
            if ((x1 < x1_min) or (x2 < x2_min)): 
                continue
        
            # test weighted four-vectors
            leg1_x1, leg2_x2 = leg1*(1/x1), leg2*(1/x2)
            ditau_test = vec.obj(px=leg1_x1.px+leg2_x2.px,
                                 py=leg1_x1.py+leg2_x2.py,
                                 pz=leg1_x1.pz+leg2_x2.pz,
                                 E=leg1_x1.E+leg2_x2.E)
            nu_test = vec.obj(px=ditau_test.px-leg1.px-leg2.px, 
                              py=ditau_test.py-leg1.py-leg2.py,
                              pz=ditau_test.pz-leg1.pz-leg2.pz,
                              E=ditau_test.E-leg1.E-leg2.E)
            test_mass = ditau_test.mass
            #  obsolete option
            #                if (((test_mass < 124.0) or
            #                     (test_mass > 126.0)) and
            #                    constrain): continue 
            
            # calculate mass likelihood integral 
            m_shift = test_mass * delta
            if (m_shift < m_vis): continue 
            x1_min = min(1.0, math.pow((m_vis_1/m_tau),2))
            x2_min = max(math.pow((m_vis_2/m_tau),2), 
                         math.pow((m_vis/m_shift),2))
            x2_max = min(1.0, math.pow((m_vis/m_shift),2)/x1_min)
            if (x2_max < x2_min): continue
            J = 2*math.pow(m_vis,2) * math.pow(m_shift, -reg_order)
            I_x2 = math.log(x2_max) - math.log(x2_min)
            I_tot = I_x2
            if (decay_type_1[i] != 2):
                I_m_nunu_1 = math.pow((m_vis/m_shift),2) * (math.pow(x2_max,-1) - math.pow(x2_min,-1))
                I_tot += I_m_nunu_1
            if (decay_type_2[i] != 2):
                I_m_nunu_2 = math.pow((m_vis/m_shift),2) * I_x2 - (x2_max - x2_min)
                I_tot += I_m_nunu_2
            mass_likelihood = math.log(J) + math.log(I_tot)
                
            # calculate MET transfer function 
            residual_x = met_x[i] - nu_test.x
            residual_y = met_y[i] - nu_test.y
            pull2 = (residual_x*(metcovinv_xx*residual_x + 
                                 metcovinv_xy*residual_y) +
                     residual_y*(metcovinv_yx*residual_x +
                                 metcovinv_yy*residual_y))
            pull2 /= metcovinv_det
            met_transfer = math.log(met_const) - 0.5*pull2
                
            likelihood_fastmtt = -met_transfer - mass_likelihood

            # perform scan over phi - azimuth around direction
            # of visible tau2 momentum (a1 resonance)
            Ptau = math.sqrt(leg2_x2.E*leg2_x2.E-m_tau*m_tau)
            Pvis = math.sqrt(leg2.E*leg2.E-leg2.mass*leg2.mass)
            cosTheta_GJ = 0.5*(2*leg2_x2.E*leg2.E - leg2.mass*leg2.mass - m_tau*m_tau)/(Ptau*Pvis)
            sinTheta_GJ = math.sqrt(1.0-cosTheta_GJ*cosTheta_GJ)
            for phi in np.arange(step_phi,twopi,step_phi):
                # rotate direction around visible momentum
                # (zprime ort defined earlier in this code) 
                # produce unit vector
                ndir_zprime = cosTheta_GJ
                ndir_xprime = sinTheta_GJ*math.cos(phi)
                ndir_yprime = sinTheta_GJ*math.sin(phi)
                ndir = ndir_xprime*xprime + ndir_yprime*yprime + ndir_zprime*zprime
                ndiff = ndir - ntau
                conv = VectorMatrixConv(ndiff,ntau_inv)
                ntau_likelihood = -0.5*conv
                likelihood = -met_transfer - mass_likelihood - ntau_likelihood
                if (likelihood < min_likelihood):
                    min_likelihood = likelihood
                    PVecTau = Ptau*ndir
                    pX = PVecTau.X()
                    pY = PVecTau.Y()
                    pZ = PVecTau.Z()
                    x1_opt, x2_opt, px_opt, py_opt, py_opt = x1, x2, pX, pY, pZ
                
        px_out[i] = px_opt
        py_out[i] = py_opt
        pz_out[i] = pz_opt

        
    return {'x1': x1_opt,
            'x2': x2_opt
            'px2': px_out,
            'py2': py_out,
            'pz2': pz_out}

# fastmtt with mass constraint
@nb.jit(nopython=True, parallel=False)
def fastmtt(pt_1, eta_1, phi_1, mass_1, decay_type_1,
            pt_2, eta_2, phi_2, mass_2, decay_type_2,
            met_x, met_y, metcov_xx, metcov_xy, metcov_yx, metcov_yy,
            verbosity=-1, delta=1/1.15, reg_order=6,
            Mass=125.0,
            Width=2.0,
            nsteps=50.):

    # step size in scan
    step = 1.0/nsteps
    
    # initialize global parameters
    m_ele = 0.51100e-3
    m_muon = 0.10566
    m_tau = 1.77685
    m_pion = 0.13957
    mass_dict = {0: m_ele, 1: m_muon}

    mass_min = Mass - Width
    mass_max = Mass + Width
    
    # initialize higgs->ditau decays, tau decay types
    N = len(pt_1)
    m_tt_vis = np.zeros(N, dtype=np.float32)
    m_tt_opt = np.zeros(N, dtype=np.float32)
    m_tt_opt_cons = np.zeros(N, dtype=np.float32)
    x1_out = np.zeros(N, dtype=np.float32)
    x2_out = np.zeros(N, dtype=np.float32)
    x1_out_cons = np.zeros(N, dtype=np.float32)
    x2_out_cons = np.zeros(N, dtype=np.float32)
    x1_out_hard = np.zeros(N, dtype=np.float32)
    x2_out_hard = np.zeros(N, dtype=np.float32)
    
    # loop over all events, calculate corrected ditau mass
    for i in range(N):
        
        # grab the correct masses based on tau decay type
        # tau decay_type: 0 ==> leptonic to electron, 
        #                 1 ==> leptonic to muon, 
        #                 2 ==> leptonic to hadronic
        if (decay_type_1[i] != 2): m1 = mass_dict[decay_type_1[i]]
        else: m1 = mass1[i]
        if (decay_type_2[i] != 2): m2 = mass_dict[decay_type_2[i]]
        else: m2 = mass2[i]
            
        # store visible masses
        m_vis_1 = m1
        m_vis_2 = m2
        
        # determine minimum and maximum possible masses
        m_vis_min_1, m_vis_max_1 = 0, 0
        m_vis_min_2, m_vis_max_2 = 0, 0
        if (decay_type_1[i] == 0): m_vis_min_1, m_vis_max_1 = m_ele, m_ele
        if (decay_type_1[i] == 1): m_vis_min_1, m_vis_max_1 = m_muon, m_muon
        if (decay_type_1[i] == 2): m_vis_min_1, m_vis_max_1 = m_pion, 1.5
        if (decay_type_2[i] == 0): m_vis_min_2, m_vis_max_2 = m_ele, m_ele
        if (decay_type_2[i] == 1): m_vis_min_2, m_vis_max_2 = m_muon, m_muon
        if (decay_type_2[i] == 2): m_vis_min_2, m_vis_max_2 = m_pion, 1.5
        if (m_vis_1 < m_vis_min_1): m_vis_1 = m_vis_min_1
        if (m_vis_1 > m_vis_max_1): m_vis_1 = m_vis_max_1
        if (m_vis_2 < m_vis_min_2): m_vis_2 = m_vis_min_2
        if (m_vis_2 > m_vis_max_2): m_vis_2 = m_vis_max_2
          
        # store both tau candidate four vectors
        leg1 = vec.obj(pt=pt_1[i], eta=eta_1[i], phi=phi_1[i], mass=m_vis_1)
        leg2 = vec.obj(pt=pt_2[i], eta=eta_2[i], phi=phi_2[i], mass=m_vis_2)

        # store visible mass of ditau pair
        m_vis = math.sqrt(2*leg1.pt*leg2.pt*(math.cosh(leg1.eta - leg2.eta) - 
                                             math.cos(leg1.phi - leg2.phi)))
        m_tt_vis[i] = m_vis
            
        # correct initial visible masses
        if (decay_type_1[i] == 2 and m_vis_1 > 1.5): m_vis_1 = 0.3
        if (decay_type_2[i] == 2 and m_vis_2 > 1.5): m_vis_2 = 0.3

        # invert met covariance matrix, calculate determinant
        metcovinv_xx, metcovinv_yy = metcov_yy[i], metcov_xx[i]
        metcovinv_xy, metcovinv_yx = -metcov_xy[i], -metcov_yx[i]
        metcovinv_det = (metcovinv_xx*metcovinv_yy -
                         metcovinv_yx*metcovinv_xy)
        if (metcovinv_det<1e-10): 
                print("Warning! Ill-conditioned MET covariance at event index", i)
                continue
               
        # perform likelihood scan 
        # see http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_032_v3.pdf
        met_const = 1/(2*math.pi*math.sqrt(metcovinv_det))
        x1_opt, x2_opt, mass_opt = step, step, -1.0
        x1_opt_BW, x2_opt_BW = step, step
        x1_opt_hard, x2_opt_hard = step, step
        min_likelihood, min_likelihood_BW, min_likelihood_hard = 1e9, 1e9, 1e9
        mass_likelihood, met_transfer = 0., 0.
        
        # scan over weights for each ditau four-vector
        intialise = False
        for x1 in np.arange(step, 1.00, step):
            # first implementing hard mass constraint
            x2 = m_vis*m_vis/(x1*Mass*Mass)
            x1_min = min(1, math.pow((m_vis_1/m_tau), 2))
            x2_min = min(1, math.pow((m_vis_2/m_tau),2))
            if ((x1 > x1_min) and (x2 > x2_min)): 
                # test weighted four-vectors
                leg1_x1, leg2_x2 = leg1*(1/x1), leg2*(1/x2)
                ditau_test = vec.obj(px=leg1_x1.px+leg2_x2.px,
                                     py=leg1_x1.py+leg2_x2.py,
                                     pz=leg1_x1.pz+leg2_x2.pz,
                                     E=leg1_x1.E+leg2_x2.E)
                nu_test = vec.obj(px=ditau_test.px-leg1.px-leg2.px, 
                                  py=ditau_test.py-leg1.py-leg2.py,
                                  pz=ditau_test.pz-leg1.pz-leg2.pz,
                                  E=ditau_test.E-leg1.E-leg2.E)
                test_mass = ditau_test.mass

                m_shift = test_mass * delta
                if (m_shift < m_vis): continue 
                x1_min = min(1.0, math.pow((m_vis_1/m_tau),2))
                x2_min = max(math.pow((m_vis_2/m_tau),2), 
                             math.pow((m_vis/m_shift),2))
                x2_max = min(1.0, math.pow((m_vis/m_shift),2)/x1_min)
                if (x2_max < x2_min): continue
                J = 2*math.pow(m_vis,2) * math.pow(m_shift, -reg_order)
                I_x2 = math.log(x2_max) - math.log(x2_min)
                I_tot = I_x2
                if (decay_type_1[i] != 2):
                    I_m_nunu_1 = math.pow((m_vis/m_shift),2) * (math.pow(x2_max,-1) - math.pow(x2_min,-1))
                    I_tot += I_m_nunu_1
                if (decay_type_2[i] != 2):
                    I_m_nunu_2 = math.pow((m_vis/m_shift),2) * I_x2 - (x2_max - x2_min)
                    I_tot += I_m_nunu_2
                mass_likelihood = math.log(J) + math.log(I_tot)
                
                # calculate MET transfer function 
                residual_x = met_x[i] - nu_test.x
                residual_y = met_y[i] - nu_test.y
                pull2 = (residual_x*(metcovinv_xx*residual_x + 
                                     metcovinv_xy*residual_y) +
                         residual_y*(metcovinv_yx*residual_x +
                                     metcovinv_yy*residual_y))
                pull2 /= metcovinv_det
                met_transfer = math.log(met_const) - 0.5*pull2
                
                # calculate final likelihood, store if minimum
                likelihood = -met_transfer - mass_likelihood
                if (likelihood < min_likelihood_hard):
                    min_likelihood_hard = likelihood
                    x1_opt_hard, x2_opt_hard = x1, x2
            
            for x2 in np.arange(step, 1.00, step):
                x1_min = min(1, math.pow((m_vis_1/m_tau), 2))
                x2_min = min(1, math.pow((m_vis_2/m_tau),2))
                if ((x1 < x1_min) or (x2 < x2_min)): 
                    continue
        
                # test weighted four-vectors
                leg1_x1, leg2_x2 = leg1*(1/x1), leg2*(1/x2)
                ditau_test = vec.obj(px=leg1_x1.px+leg2_x2.px,
                                     py=leg1_x1.py+leg2_x2.py,
                                     pz=leg1_x1.pz+leg2_x2.pz,
                                     E=leg1_x1.E+leg2_x2.E)
                nu_test = vec.obj(px=ditau_test.px-leg1.px-leg2.px, 
                                  py=ditau_test.py-leg1.py-leg2.py,
                                  pz=ditau_test.pz-leg1.pz-leg2.pz,
                                  E=ditau_test.E-leg1.E-leg2.E)
                test_mass = ditau_test.mass
                #  obsolete option
                #                if (((test_mass < 124.0) or
                #                     (test_mass > 126.0)) and
                #                    constrain): continue 
            
                # calculate mass likelihood integral 
                m_shift = test_mass * delta
                if (m_shift < m_vis): continue 
                x1_min = min(1.0, math.pow((m_vis_1/m_tau),2))
                x2_min = max(math.pow((m_vis_2/m_tau),2), 
                             math.pow((m_vis/m_shift),2))
                x2_max = min(1.0, math.pow((m_vis/m_shift),2)/x1_min)
                if (x2_max < x2_min): continue
                J = 2*math.pow(m_vis,2) * math.pow(m_shift, -reg_order)
                I_x2 = math.log(x2_max) - math.log(x2_min)
                I_tot = I_x2
                if (decay_type_1[i] != 2):
                    I_m_nunu_1 = math.pow((m_vis/m_shift),2) * (math.pow(x2_max,-1) - math.pow(x2_min,-1))
                    I_tot += I_m_nunu_1
                if (decay_type_2[i] != 2):
                    I_m_nunu_2 = math.pow((m_vis/m_shift),2) * I_x2 - (x2_max - x2_min)
                    I_tot += I_m_nunu_2
                mass_likelihood = math.log(J) + math.log(I_tot)
                
                # calculate MET transfer function 
                residual_x = met_x[i] - nu_test.x
                residual_y = met_y[i] - nu_test.y
                pull2 = (residual_x*(metcovinv_xx*residual_x + 
                                     metcovinv_xy*residual_y) +
                         residual_y*(metcovinv_yx*residual_x +
                                     metcovinv_yy*residual_y))
                pull2 /= metcovinv_det
                met_transfer = math.log(met_const) - 0.5*pull2
                
                # calculate final likelihood, store if minimum
                likelihood = -met_transfer - mass_likelihood
                if (likelihood < min_likelihood):
                    min_likelihood = likelihood
                    x1_opt, x2_opt, mass_opt = x1, x2, test_mass
                if (test_mass>mass_min and test_mass<mass_max):
                    if (likelihood < min_likelihood_BW):
                        min_likelihood_BW = likelihood
                        x1_opt_BW, x2_opt_BW = x1, x2
                
        m_tt_opt[i] = mass_opt
        x1_out[i] = x1_opt
        x2_out[i] = x2_opt
        x1_out_cons[i] = x1_opt_BW
        x2_out_cons[i] = x2_opt_BW
        x1_out_hard[i] = x1_opt_hard
        x2_out_hard[i] = x2_opt_hard
        
    return {'mtt': m_tt_opt,
            'x_1': x1_out,
            'x_2': x2_out,
            'x_1_cons': x1_out_cons,
            'x_2_cons': x2_out_cons,
            'x_1_hard': x1_out_hard,
            'x_2_hard': x2_out_hard}



