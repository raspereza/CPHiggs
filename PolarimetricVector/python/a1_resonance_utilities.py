import awkward
import coffea.nanoevents
import numpy as np
import coffea.nanoevents.methods.vector


def _quadratic_alternate(a: awkward.Array, b: awkward.Array, c: awkward.Array):
    D = b**2 - 4.0 * a * c

    masked_D = awkward.where(D < 0, 0, D)
    q = np.where(D < 0, -0.5 * b, -0.5 * (b + np.copysign(np.sqrt(masked_D), b)))

    x_0 = c / q
    x_1 = q / a

    return x_0, x_1


def _mag2(v: awkward.Array, type) -> awkward.Array:
    if type == "ThreeVector":
        return v.x**2 + v.y**2 + v.z**2
    elif type == "LorentzVector":
        return v.t**2 - (v.x**2 + v.y**2 + v.z**2)
    return None


def _dot(v1: awkward.Array, v2: awkward.Array) -> awkward.Array:
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z


def _compute_angle(v1: awkward.Array, v2: awkward.Array) -> awkward.Array:
    ptot2 = _mag2(v1, "ThreeVector") * _mag2(v2, "ThreeVector")
    dot_product = _dot(v1, v2)
    arg = awkward.where(ptot2 > 0, dot_product / np.sqrt(ptot2), 0)
    arg = np.clip(arg, -1.0, 1.0)
    angle = np.arccos(arg)
    return awkward.where(ptot2 <= 0, 0.0, angle)


def tauMomentumSolutions(tauDir: awkward.Array, a1LV: awkward.Array):
    tauMass = 1.77682
    thetaGJ = _compute_angle(a1LV, tauDir)
    tauDirUnit = tauDir.unit

    a = 4.0 * (_mag2(a1LV, "LorentzVector") + np.sqrt(_mag2(a1LV, "ThreeVector"))**2 * np.sin(thetaGJ)**2)
    b = -4.0 * (_mag2(a1LV, "LorentzVector") + tauMass**2) * np.sqrt(_mag2(a1LV, "ThreeVector")) * np.cos(thetaGJ)
    c = 4.0 * tauMass**2 * (_mag2(a1LV, "LorentzVector") + np.sqrt(_mag2(a1LV, "ThreeVector"))**2) - (_mag2(a1LV, "LorentzVector") + tauMass**2)**2

    tauMomentumSmall, tauMomentumLarge = _quadratic_alternate(a,b,c)
    tauMomentumMean = (tauMomentumSmall + tauMomentumLarge) / 2.0

    tauSmall = awkward.zip({'x': tauMomentumSmall * tauDirUnit.x, 'y': tauMomentumSmall * tauDirUnit.y, 'z': tauMomentumSmall * tauDirUnit.z, 't': np.sqrt(tauMomentumSmall**2 + tauMass**2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    tauLarge = awkward.zip({'x': tauMomentumLarge * tauDirUnit.x, 'y': tauMomentumLarge * tauDirUnit.y, 'z': tauMomentumLarge * tauDirUnit.z, 't': np.sqrt(tauMomentumLarge**2 + tauMass**2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
    tauMean = awkward.zip({'x': tauMomentumMean * tauDirUnit.x, 'y': tauMomentumMean * tauDirUnit.y, 'z': tauMomentumMean * tauDirUnit.z, 't': np.sqrt(tauMomentumMean**2 + tauMass**2)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)

    solutions = [tauSmall, tauLarge, tauMean]
    return solutions


def tauPairConstraint(tau1Dir: awkward.Array, a1LV1: awkward.Array, tau2Dir: awkward.Array, a1LV2: awkward.Array):
    Hmass = 125.10

    tau1Solutions = tauMomentumSolutions(tau1Dir, a1LV1)
    tau2Solutions = tauMomentumSolutions(tau2Dir, a1LV2)

    mass0 = (tau1Solutions[0] + tau2Solutions[0]).mass  # tau1:small solution -- tau2:small solution
    mass1 = (tau1Solutions[0] + tau2Solutions[1]).mass  # tau1:small solution -- tau2:large solution
    mass2 = (tau1Solutions[1] + tau2Solutions[0]).mass  # tau1:large solution -- tau2:small solution
    mass3 = (tau1Solutions[1] + tau2Solutions[1]).mass  # tau1:large solution -- tau2:large solution

    massPairs = awkward.Array([mass0, mass1, mass2, mass3])

    massConstraint = [np.abs(massPairs[0] - Hmass), np.abs(massPairs[1] - Hmass), np.abs(massPairs[2] - Hmass), np.abs(massPairs[3] - Hmass)]
    bestCouple = np.argmin(massConstraint, axis=0)

    tau1PairConstraintLV = awkward.where(bestCouple == 0, tau1Solutions[0], awkward.where(bestCouple == 1, tau1Solutions[0], awkward.where(bestCouple == 2, tau1Solutions[1], tau1Solutions[1])))

    tau2PairConstraintLV = awkward.where(bestCouple == 0, tau2Solutions[0], awkward.where(bestCouple == 1, tau2Solutions[1], awkward.where(bestCouple == 2, tau2Solutions[0], tau2Solutions[1])))

    return tau1PairConstraintLV, tau2PairConstraintLV


def tauPairMomentumSolutions(tau1Dir: awkward.Array, a1LV1: awkward.Array, tau2Dir: awkward.Array, a1LV2: awkward.Array):
    tau1PairConstraintLV, tau2PairConstraintLV = tauPairConstraint(tau1Dir, a1LV1, tau2Dir, a1LV2)
    return tau1PairConstraintLV, tau2PairConstraintLV
