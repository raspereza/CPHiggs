import numpy as np
import awkward
import coffea
from higgs_dna.tools.htautau_cp.a1_resonance_utilities import tauPairMomentumSolutions
from higgs_dna.tools.htautau_cp.PolarimetricA1 import PolarimetricA1
import logging

logger = logging.getLogger(__name__)


def acoplanarity_calculator(pairs: awkward.Array, method_leg_1: str, method_leg_2: str, mode_leg_1:str, mode_leg_2:str, use_svfit=False, use_fastmtt=False) -> awkward.Array:
    candidate_dict = {}
    candidate_dict["pv"] = awkward.zip({"x": pairs.PVBS_x, "y": pairs.PVBS_y, "z": pairs.PVBS_z, "t": awkward.zeros_like(pairs.PVBS_x)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["ip_1"] = awkward.zip({"x": pairs.ip_x_1, "y": pairs.ip_y_1, "z": pairs.ip_z_1, "t": awkward.zeros_like(pairs.ip_x_1)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["ip_2"] = awkward.zip({"x": pairs.ip_x_2, "y": pairs.ip_y_2, "z": pairs.ip_z_2, "t": awkward.zeros_like(pairs.ip_x_2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["sv_1"] = awkward.zip({"x": pairs.sv_x_1, "y": pairs.sv_y_1, "z": pairs.sv_z_1, "t": awkward.zeros_like(pairs.sv_x_1)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["sv_2"] = awkward.zip({"x": pairs.sv_x_2, "y": pairs.sv_y_2, "z": pairs.sv_z_2, "t": awkward.zeros_like(pairs.sv_x_2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["p1"] = awkward.zip({"pt":pairs['obj_1'].pt, "eta":pairs['obj_1'].eta, "phi":pairs['obj_1'].phi, "mass":pairs['obj_1'].mass, "charge": pairs['obj_1'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["p2"] = awkward.zip({"pt":pairs['obj_2'].pt, "eta":pairs['obj_2'].eta, "phi":pairs['obj_2'].phi, "mass":pairs['obj_2'].mass, "charge": pairs['obj_2'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi_1"] = awkward.zip({"pt":pairs.pi_pt_1, "eta":pairs.pi_eta_1, "phi":pairs.pi_phi_1, "mass":pairs.pi_mass_1, "charge": pairs.pi_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi_2"] = awkward.zip({"pt":pairs.pi_pt_2, "eta":pairs.pi_eta_2, "phi":pairs.pi_phi_2, "mass":pairs.pi_mass_2, "charge": pairs.pi_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi2_1"] = awkward.zip({"pt":pairs.pi2_pt_1, "eta":pairs.pi2_eta_1, "phi":pairs.pi2_phi_1, "mass":pairs.pi2_mass_1, "charge": pairs.pi2_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi2_2"] = awkward.zip({"pt":pairs.pi2_pt_2, "eta":pairs.pi2_eta_2, "phi":pairs.pi2_phi_2, "mass":pairs.pi2_mass_2, "charge": pairs.pi2_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi3_1"] = awkward.zip({"pt":pairs.pi3_pt_1, "eta":pairs.pi3_eta_1, "phi":pairs.pi3_phi_1, "mass":pairs.pi3_mass_1, "charge": pairs.pi3_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi3_2"] = awkward.zip({"pt":pairs.pi3_pt_2, "eta":pairs.pi3_eta_2, "phi":pairs.pi3_phi_2, "mass":pairs.pi3_mass_2, "charge": pairs.pi3_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi0_1"] = awkward.zip({"pt":pairs.pi0_pt_1, "eta":pairs.pi0_eta_1, "phi":pairs.pi0_phi_1, "mass":pairs.pi0_mass_1, "charge": pairs.pi0_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi0_2"] = awkward.zip({"pt":pairs.pi0_pt_2, "eta":pairs.pi0_eta_2, "phi":pairs.pi0_phi_2, "mass":pairs.pi0_mass_2, "charge": pairs.pi0_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["sv1-pv"] = candidate_dict["sv_1"] - candidate_dict["pv"]
    candidate_dict["sv2-pv"] = candidate_dict["sv_2"] - candidate_dict["pv"]

    obj_list = ['p1','p2','pi_1','pi_2','pi2_1','pi2_2','pi3_1','pi3_2','pi0_1','pi0_2']
    if use_svfit:
        candidate_dict["svfit_tau_1"] = awkward.zip({"pt":pairs.svfit_tau1_pt, "eta":candidate_dict["sv1-pv"].eta, "phi":candidate_dict["sv1-pv"].phi, "mass": awkward.ones_like(pairs.svfit_tau1_pt) * 1.777, "charge": pairs['obj_1'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
        candidate_dict["svfit_tau_2"] = awkward.zip({"pt":pairs.svfit_tau2_pt, "eta":candidate_dict["sv2-pv"].eta, "phi":candidate_dict["sv2-pv"].phi, "mass": awkward.ones_like(pairs.svfit_tau2_pt) * 1.777, "charge": pairs['obj_2'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
        obj_list.append('svfit_tau_1')
        obj_list.append('svfit_tau_2')

    if use_fastmtt:
        candidate_dict["fastmtt_tau_1"] = awkward.zip({"pt":pairs.FastMTT_Tau1Pt, "eta":candidate_dict["sv1-pv"].eta, "phi":candidate_dict["sv1-pv"].phi, "mass": awkward.ones_like(pairs.FastMTT_Tau1Pt) * 1.777, "charge": pairs['obj_1'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
        candidate_dict["fastmtt_tau_2"] = awkward.zip({"pt":pairs.FastMTT_Tau2Pt, "eta":candidate_dict["sv2-pv"].eta, "phi":candidate_dict["sv2-pv"].phi, "mass": awkward.ones_like(pairs.FastMTT_Tau2Pt) * 1.777, "charge": pairs['obj_2'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
        obj_list.append('fastmtt_tau_1')
        obj_list.append('fastmtt_tau_2')

    for obj in obj_list:
        candidate_dict[obj] = convert_to_cartesian(candidate_dict[obj])

    R1, P1 = prepareVectors(candidate_dict, method_leg_1, mode_leg_1, "1")
    R2, P2 = prepareVectors(candidate_dict, method_leg_2, mode_leg_2, "2")

    if method_leg_1 == "Polarimetric-Vector-A1A1" and method_leg_2 == "Polarimetric-Vector-A1A1":

        R1_new = awkward.where(pairs['obj_1'].charge == -1, R1, R2)
        R2_new = awkward.where(pairs['obj_2'].charge == 1, R2, R1)
        P1_new = awkward.where(pairs['obj_1'].charge == -1, P1, P2)
        P2_new = awkward.where(pairs['obj_2'].charge == 1, P2, P1)

        del R1, R2, P1, P2

        R1, R2, P1, P2 = compute_PV_A1A1(R1_new, R2_new, P1_new, P2_new)

    elif (method_leg_1 == "Polarimetric-Vector-SVFIT-A1" and use_svfit) or (method_leg_1 == "Polarimetric-Vector-FASTMTT-A1" and use_fastmtt):

        R1, P1 = compute_PV(R1, P1, method_leg_1, mode_leg_1)

    elif (method_leg_2 == "Polarimetric-Vector-SVFIT-A1" and use_svfit) or (method_leg_2 == "Polarimetric-Vector-FASTMTT-A1" and use_fastmtt):

        R2, P2 = compute_PV(R2, P2, method_leg_2, mode_leg_2)

    return computeAcoplanarity(R1, R2, P1, P2, method_leg_1, method_leg_2)


def computeAcoplanarity(R1: awkward.Array, R2: awkward.Array, P1: awkward.Array, P2: awkward.Array, method_leg_1:str, method_leg_2:str) -> awkward.Array:
    if method_leg_1 == "Polarimetric-Vector-A1A1" and method_leg_2 == "Polarimetric-Vector-A1A1":
        return computeAcoplanarity_PV(R1, R2, P1, P2)
    else:
        boost = (P1 + P2).boostvec

        if method_leg_1 == "Impact-Parameter":
            R1 = R1.unit
        if method_leg_2 == "Impact-Parameter":
            R2 = R2.unit

        R1_boosted = R1.boost(boost.negative())
        R2_boosted = R2.boost(boost.negative())
        P1_boosted = P1.boost(boost.negative())
        P2_boosted = P2.boost(boost.negative())

        n1 = R1_boosted.pvec - R1_boosted.pvec.dot(P1_boosted.pvec.unit) * P1_boosted.pvec.unit
        n2 = R2_boosted.pvec - R2_boosted.pvec.dot(P2_boosted.pvec.unit) * P2_boosted.pvec.unit

        n1 = n1.unit
        n2 = n2.unit

        angle = np.arccos(n1.dot(n2))
        sign = P2_boosted.pvec.unit.dot(n1.cross(n2))

        angle = awkward.where(sign >= 0, angle, 2 * np.pi - angle)

        cp_sign = np.ones_like(angle)

        if method_leg_1 == "Decay-Plane":
            Y1 = (P1_boosted.t - R1_boosted.t) / (P1_boosted.t + R1_boosted.t)
            cp_sign = cp_sign * awkward.to_numpy(Y1)
        elif method_leg_1 == "Decay-Plane-a1":
            Y1 = (R1_boosted.t - P1_boosted.t) / (R1_boosted.t + P1_boosted.t)
            cp_sign = cp_sign * awkward.to_numpy(Y1)

        if method_leg_2 == "Decay-Plane":
            Y2 = (P2_boosted.t - R2_boosted.t) / (P2_boosted.t + R2_boosted.t)
            cp_sign = cp_sign * awkward.to_numpy(Y2)
        elif method_leg_2 == "Decay-Plane-a1":
            Y2 = (R2_boosted.t - P2_boosted.t) / (R2_boosted.t + P2_boosted.t)
            cp_sign = cp_sign * awkward.to_numpy(Y2)

        if method_leg_1 == "Decay-Plane" or method_leg_2 == "Decay-Plane" or method_leg_1 == "Decay-Plane-a1" or method_leg_2 == "Decay-Plane-a1":
            angle = awkward.where(cp_sign < 0, awkward.where(angle < np.pi, angle + np.pi, angle - np.pi), angle)
        return angle


def computeAcoplanarity_PV(R1: awkward.Array, R2: awkward.Array, P1: awkward.Array, P2: awkward.Array) -> awkward.Array:

    h1 = R1.unit
    h2 = R2.unit
    n1 = P1.pvec.unit
    n2 = P2.pvec.unit

    k1 = h1.cross(n1).unit
    k2 = h2.cross(n2).unit

    angle = np.arctan2(k1.cross(k2).absolute(), k1.dot(k2))
    sign = h1.cross(h2).dot(n1)

    angle = awkward.where(sign <= 0, angle, 2 * np.pi - angle)

    return angle


def compute_PV(R: awkward.Array, P: awkward.Array, method_leg: str, mode_leg: str) -> awkward.Array:
    if mode_leg == "a1":
        tau = R
        visible_tau = P["0"] + P["1"] + P["2"]
        tau = rotate_to_GJMax(visible_tau, tau)
        frame = tau.boostvec
        os_pi_HRF = P["0"].boost(frame.negative())
        ss1_pi_HRF = P["1"].boost(frame.negative())
        ss2_pi_HRF = P["2"].boost(frame.negative())
        a1_pol = PolarimetricA1(tau.boost(frame.negative()), os_pi_HRF, ss1_pi_HRF, ss2_pi_HRF, tau.charge)
        pv = -a1_pol.PVC().pvec
        pv = awkward.zip({"x": pv.x, "y": pv.y, "z": pv.z, "t": awkward.zeros_like(pv.x)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
        pv = pv.boost(frame)
        return pv, tau


def compute_PV_A1A1(R1: awkward.Array, R2: awkward.Array, P1: awkward.Array, P2: awkward.Array) -> awkward.Array:

    a1_1 = P1["0"] + P1["1"] + P1["2"]
    a1_2 = P2["0"] + P2["1"] + P2["2"]

    Tauminus, Tauplus = tauPairMomentumSolutions(R1, a1_1, R2, a1_2)

    frame = (Tauminus + Tauplus).boostvec

    os_pi_HRF_1 = P1["0"].boost(frame.negative())
    ss1_pi_HRF_1 = P1["1"].boost(frame.negative())
    ss2_pi_HRF_1 = P1["2"].boost(frame.negative())

    os_pi_HRF_2 = P2["0"].boost(frame.negative())
    ss1_pi_HRF_2 = P2["1"].boost(frame.negative())
    ss2_pi_HRF_2 = P2["2"].boost(frame.negative())

    a1_pol_1 = PolarimetricA1(Tauminus.boost(frame.negative()), os_pi_HRF_1, ss1_pi_HRF_1, ss2_pi_HRF_1, -1)
    pv_1 = -a1_pol_1.PVC().pvec

    a1_pol_2 = PolarimetricA1(Tauplus.boost(frame.negative()), os_pi_HRF_2, ss1_pi_HRF_2, ss2_pi_HRF_2, +1)
    pv_2 = -a1_pol_2.PVC().pvec

    return pv_1, pv_2, Tauminus.boost(frame.negative()), Tauplus.boost(frame.negative())


def alphaAngle_calculator(pairs: awkward.Array, method_leg_1: str, method_leg_2: str, mode_leg_1:str, mode_leg_2:str, index: str) -> awkward.Array:

    candidate_dict = {}
    candidate_dict["pv"] = awkward.zip({"x": pairs.PVBS_x, "y": pairs.PVBS_y, "z": pairs.PVBS_z, "t": awkward.zeros_like(pairs.PVBS_x)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["ip_1"] = awkward.zip({"x": pairs.ip_x_1, "y": pairs.ip_y_1, "z": pairs.ip_z_1, "t": awkward.zeros_like(pairs.ip_x_1)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["ip_2"] = awkward.zip({"x": pairs.ip_x_2, "y": pairs.ip_y_2, "z": pairs.ip_z_2, "t": awkward.zeros_like(pairs.ip_x_2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["sv_1"] = awkward.zip({"x": pairs.sv_x_1, "y": pairs.sv_y_1, "z": pairs.sv_z_1, "t": awkward.zeros_like(pairs.sv_x_1)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["sv_2"] = awkward.zip({"x": pairs.sv_x_2, "y": pairs.sv_y_2, "z": pairs.sv_z_2, "t": awkward.zeros_like(pairs.sv_x_2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["p1"] = awkward.zip({"pt":pairs['obj_1'].pt, "eta":pairs['obj_1'].eta, "phi":pairs['obj_1'].phi, "mass":pairs['obj_1'].mass, "charge": pairs['obj_1'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["p2"] = awkward.zip({"pt":pairs['obj_2'].pt, "eta":pairs['obj_2'].eta, "phi":pairs['obj_2'].phi, "mass":pairs['obj_2'].mass, "charge": pairs['obj_2'].charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi_1"] = awkward.zip({"pt":pairs.pi_pt_1, "eta":pairs.pi_eta_1, "phi":pairs.pi_phi_1, "mass":pairs.pi_mass_1, "charge": pairs.pi_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi_2"] = awkward.zip({"pt":pairs.pi_pt_2, "eta":pairs.pi_eta_2, "phi":pairs.pi_phi_2, "mass":pairs.pi_mass_2, "charge": pairs.pi_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi2_1"] = awkward.zip({"pt":pairs.pi2_pt_1, "eta":pairs.pi2_eta_1, "phi":pairs.pi2_phi_1, "mass":pairs.pi2_mass_1, "charge": pairs.pi2_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi2_2"] = awkward.zip({"pt":pairs.pi2_pt_2, "eta":pairs.pi2_eta_2, "phi":pairs.pi2_phi_2, "mass":pairs.pi2_mass_2, "charge": pairs.pi2_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi3_1"] = awkward.zip({"pt":pairs.pi3_pt_1, "eta":pairs.pi3_eta_1, "phi":pairs.pi3_phi_1, "mass":pairs.pi3_mass_1, "charge": pairs.pi3_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi3_2"] = awkward.zip({"pt":pairs.pi3_pt_2, "eta":pairs.pi3_eta_2, "phi":pairs.pi3_phi_2, "mass":pairs.pi3_mass_2, "charge": pairs.pi3_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi0_1"] = awkward.zip({"pt":pairs.pi0_pt_1, "eta":pairs.pi0_eta_1, "phi":pairs.pi0_phi_1, "mass":pairs.pi0_mass_1, "charge": pairs.pi0_charge_1}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    candidate_dict["pi0_2"] = awkward.zip({"pt":pairs.pi0_pt_2, "eta":pairs.pi0_eta_2, "phi":pairs.pi0_phi_2, "mass":pairs.pi0_mass_2, "charge": pairs.pi0_charge_2}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)

    for obj in ['p1','p2','pi_1','pi_2','pi2_1','pi2_2','pi3_1','pi3_2','pi0_1','pi0_2']:
        candidate_dict[obj] = convert_to_cartesian(candidate_dict[obj])

    if method_leg_1 == "Impact-Parameter" and (mode_leg_1 == "e" or mode_leg_1 == "mu"):
        R1 = candidate_dict["ip_1"].pvec.unit
        P1 = candidate_dict["p1"].pvec
        if method_leg_2 == "Impact-Parameter" and mode_leg_2 == "pi":
            R2 = candidate_dict["ip_2"].pvec.unit
            P2 = candidate_dict["pi_2"].pvec

        if index == 1:
            alphaAngle = computeAlphaAngle(R1, P1)
        elif index == 2:
            alphaAngle = computeAlphaAngle(R2, P2)
        else:
            raise RuntimeError("Invalid index for alpha angle calculation")
        return alphaAngle
    else:
        raise RuntimeError("Only Impact-Parameter method is implemented for alpha angle calculation")


def computeAlphaAngle(R: awkward.Array, P: awkward.Array) -> awkward.Array:
    z = awkward.zip({"x": awkward.zeros_like(R.x), "y": awkward.zeros_like(R.y), "z": awkward.ones_like(R.z)}, with_name="ThreeVector", behavior=coffea.nanoevents.methods.vector.behavior)
    R = R.unit
    P = P.unit
    cross1 = z.cross(P)
    cross2 = R.cross(P)
    angle = np.arccos(np.abs(cross1.dot(cross2) / (cross1.absolute() * cross2.absolute())))
    return angle


def prepareVectors(candidate_dict: dict, method_leg:str, mode_leg: str, index_leg: str) -> awkward.Array:
    if method_leg == "Impact-Parameter":
        R = candidate_dict[f"ip_{index_leg}"]
        if mode_leg == "e" or mode_leg == "mu":
            P = candidate_dict[f"p{index_leg}"]
        elif mode_leg == "pi":
            P = candidate_dict[f"pi_{index_leg}"]
        elif mode_leg == "rho" or mode_leg == "a1":
            raise RuntimeError("Rho and A1 decay modes not implemented")
        else:
            raise RuntimeError("Invalid leg mode")

    elif method_leg == "Decay-Plane":
        R = candidate_dict[f"pi0_{index_leg}"]
        if mode_leg == "e" or mode_leg == "mu":
            P = candidate_dict[f"p{index_leg}"]
        elif mode_leg == "rho":
            P = candidate_dict[f"pi_{index_leg}"]
        elif mode_leg == "a1" or mode_leg == "pi":
            raise RuntimeError("A1 and Pi decay modes not implemented")

    elif method_leg == "Decay-Plane-a1":
        sorted_hads = sortA1(candidate_dict[f"pi_{index_leg}"], candidate_dict[f"pi2_{index_leg}"], candidate_dict[f"pi3_{index_leg}"])
        R = sorted_hads["0"]
        P = sorted_hads["1"]

    elif method_leg == "Polarimetric-Vector-A1A1":
        sorted_hads = sortA1(candidate_dict[f"pi_{index_leg}"], candidate_dict[f"pi2_{index_leg}"], candidate_dict[f"pi3_{index_leg}"])
        R = candidate_dict[f"sv_{index_leg}"] - candidate_dict["pv"]
        R = R.pvec.unit
        P = sorted_hads

    elif method_leg == "Polarimetric-Vector-SVFIT-A1":
        sorted_hads = sortA1(candidate_dict[f"pi_{index_leg}"], candidate_dict[f"pi2_{index_leg}"], candidate_dict[f"pi3_{index_leg}"])
        R = candidate_dict[f"svfit_tau_{index_leg}"]
        P = sorted_hads

    elif method_leg == "Polarimetric-Vector-FASTMTT-A1":
        sorted_hads = sortA1(candidate_dict[f"pi_{index_leg}"], candidate_dict[f"pi2_{index_leg}"], candidate_dict[f"pi3_{index_leg}"])
        R = candidate_dict[f"fastmtt_tau_{index_leg}"]
        P = sorted_hads

    return R, P


def sortA1(pi: awkward.Array, pi2: awkward.Array, pi3: awkward.Array) -> list:

    hads = awkward.zip({"pi": pi, "pi2": pi2, "pi3": pi3})

    condition1 = (hads["pi"].charge != hads["pi2"].charge) & (hads["pi"].charge != hads["pi3"].charge)
    condition2 = (hads["pi2"].charge != hads["pi"].charge) & (hads["pi2"].charge != hads["pi3"].charge)

    sorted_hads = awkward.where(condition1, awkward.zip([hads["pi"], hads["pi2"], hads["pi3"]]), awkward.where(condition2, awkward.zip([hads["pi2"], hads["pi"], hads["pi3"]]), awkward.zip([hads["pi3"], hads["pi2"], hads["pi"]])))

    rho_mass = 0.7755
    dM1 = abs((sorted_hads["0"] + sorted_hads["1"]).mass - rho_mass)
    dM2 = abs((sorted_hads["0"] + sorted_hads["2"]).mass - rho_mass)

    condition_3 = dM2 < dM1
    sorted_hads = awkward.where(condition_3, awkward.zip([sorted_hads["0"], sorted_hads["2"], sorted_hads["1"]]), sorted_hads)
    return sorted_hads


def convert_to_cartesian(obj):

    pt = awkward.where(obj.pt == -9999.0, np.nan, obj.pt)
    eta = awkward.where(obj.eta == -9999.0, np.nan, obj.eta)
    phi = awkward.where(obj.phi == -9999.0, np.nan, obj.phi)
    mass = awkward.where(obj.mass == -9999.0, np.nan, obj.mass)

    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    E = np.sqrt(mass**2 + pt**2 * (np.cosh(eta)**2))

    cartesian_obj = awkward.zip({
        'x': awkward.where(px == np.nan, -9999.0, px),
        'y': awkward.where(py == np.nan, -9999.0, py),
        'z': awkward.where(pz == np.nan, -9999.0, pz),
        't': awkward.where(E == np.nan, -9999.0, E),
        'charge': obj.charge
    }, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)

    return cartesian_obj


def rotate_to_GJMax(visible_tau: awkward.Array, tau: awkward.Array) -> awkward.Array:

    tau_direction = tau.pvec.unit
    visible_tau_direction = visible_tau.pvec.unit

    mass_tau = tau.mass
    mass_visible_tau = visible_tau.mass

    theta_GJ = np.arccos(np.clip(tau_direction.dot(visible_tau_direction), -1, 1))
    theta_GJ_max = np.arcsin(np.clip((mass_tau**2 - mass_visible_tau**2) / (2 * mass_tau * visible_tau.p), -1, 1))

    mask = theta_GJ > theta_GJ_max

    n_1_x = 1 / np.sqrt(1 + (visible_tau.x / visible_tau.y)**2)
    n_1_y = -n_1_x * visible_tau.x / visible_tau.y
    n_1 = awkward.zip({"x": n_1_x, "y": n_1_y, "z": awkward.zeros_like(n_1_x)}, with_name="ThreeVector", behavior=coffea.nanoevents.methods.vector.behavior)
    n_2 = n_1.cross(visible_tau.pvec.unit)

    phi_opt_1 = np.arctan(tau_direction.dot(n_2) / tau_direction.dot(n_1))
    new_dir_1 = np.cos(theta_GJ_max) * visible_tau_direction + np.sin(theta_GJ_max) * (np.cos(phi_opt_1) * n_1 + np.sin(phi_opt_1) * n_2)
    phi_opt_2 = phi_opt_1 + np.pi
    new_dir_2 = np.cos(theta_GJ_max) * visible_tau_direction + np.sin(theta_GJ_max) * (np.cos(phi_opt_2) * n_1 + np.sin(phi_opt_2) * n_2)

    mask_dir = new_dir_1.dot(tau_direction) > new_dir_2.dot(tau_direction)
    new_dir = awkward.where(mask_dir, new_dir_1, new_dir_2)

    new_pt = tau.p * np.sin(np.arctan(np.exp(calculate_pseudorapidity(-new_dir))) * 2)
    new_tau = awkward.zip({"pt": new_pt, "eta": calculate_pseudorapidity(new_dir), "phi": calculate_phi(new_dir), "mass": tau.mass, "charge": tau.charge}, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    new_tau = convert_to_cartesian(new_tau)
    corrected_tau = awkward.where(mask, new_tau, tau)

    return corrected_tau


def calculate_pseudorapidity(three_vector: awkward.Array) -> awkward.Array:
    mag = np.sqrt(three_vector.x**2 + three_vector.y**2 + three_vector.z**2)
    cos_theta = three_vector.z / mag
    cos_theta = awkward.where(mag == 0, 1, cos_theta)

    eta = -0.5 * np.log((1 - cos_theta) / (1 + cos_theta))

    eta = awkward.where((cos_theta * cos_theta) < 1, eta
                        , awkward.where(three_vector.z == 0, 0,
                                        awkward.where(three_vector.z > 0, 10**10, -10**10)))
    return eta


def calculate_phi(three_vector: awkward.Array) -> awkward.Array:
    phi = np.arctan2(three_vector.y, three_vector.x)
    phi = awkward.where((three_vector.x == 0) & (three_vector.y == 0), 0, phi)
    return phi
