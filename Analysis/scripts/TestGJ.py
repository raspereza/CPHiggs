#! /usr/bin/env python3
# Author: Alexei Raspereza (January 2026)
import CPHiggs.Analysis.pv_utils as pv_utils
import CPHiggs.Analysis.utils as utils
import ROOT
import math

if __name__ == "__main__":


    E_vis = 30
    m_vis = 1.2
    P_vis = math.sqrt(E_vis*E_vis-m_vis*m_vis)

    theta_GJ_max = math.asin(max(-1.0,min(1.0,0.5*(utils.tau_mass*utils.tau_mass-m_vis*m_vis)/(utils.tau_mass*P_vis))))

    for i in range(0,10):
        theta = theta_GJ_max*(1.0-0.1*float(i))
        cos_theta = math.cos(theta)
        p = pv_utils.tauMomGJ(E_vis,m_vis,cos_theta)
        print(f"{theta} {p}")
