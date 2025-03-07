import awkward
from higgs_dna.tools.htautau_cp.acoplanarity_calculator import acoplanarity_calculator, alphaAngle_calculator
import logging
import warnings

logger = logging.getLogger(__name__)


def acoplanarity_producer(pairs: awkward.Array, channel: str, use_svfit: bool, use_fastmtt: bool) -> awkward.Array:
    logger.debug(f"Computing acoplanarity angles for channel {channel}")

    # Suppress warnings associated with non-sense LorentzVectors (As Advised by coffea devs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        # Acoplanarity angle: channel --> etau
        if channel == "et":
            pairs['aco_e_pi'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "e", "pi")
            pairs['aco_e_rho'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Decay-Plane", "e", "rho")
            pairs['aco_e_a1'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Decay-Plane-a1", "e", "a1")

            pairs['alphaAngle_e_pi_1'] = alphaAngle_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "e", "pi", 1)
            pairs['alphaAngle_e_pi_2'] = alphaAngle_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "e", "pi", 2)

        # Acoplanarity angle: channel --> mutau
        if channel == "mt":
            pairs['aco_mu_pi'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "mu", "pi")
            pairs['aco_mu_rho'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Decay-Plane", "mu", "rho")
            pairs['aco_mu_a1'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Decay-Plane-a1", "mu", "a1")

            pairs['alphaAngle_mu_pi_1'] = alphaAngle_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "mu", "pi", 1)
            pairs['alphaAngle_mu_pi_2'] = alphaAngle_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "mu", "pi", 2)

        # Acoplanarity angle: channel --> tautau
        if channel == "tt":
            pairs['aco_pi_pi'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Impact-Parameter", "pi", "pi")
            pairs['aco_pi_rho'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Decay-Plane", "pi", "rho")
            pairs['aco_pi_a1'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Decay-Plane-a1", "pi", "a1")
            pairs['aco_rho_pi'] = acoplanarity_calculator(pairs, "Decay-Plane", "Impact-Parameter", "rho", "pi")
            pairs['aco_rho_rho'] = acoplanarity_calculator(pairs, "Decay-Plane", "Decay-Plane", "rho", "rho")
            pairs['aco_rho_a1'] = acoplanarity_calculator(pairs, "Decay-Plane", "Decay-Plane-a1", "rho", "a1")
            pairs['aco_a1_pi'] = acoplanarity_calculator(pairs, "Decay-Plane-a1", "Impact-Parameter", "a1", "pi")
            pairs['aco_a1_rho'] = acoplanarity_calculator(pairs, "Decay-Plane-a1", "Decay-Plane", "a1", "rho")
            pairs['aco_a1_a1'] = acoplanarity_calculator(pairs, "Polarimetric-Vector-A1A1", "Polarimetric-Vector-A1A1", "a1", "a1")

            if use_svfit:
                pairs['aco_pi_a1_SVFIT'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Polarimetric-Vector-SVFIT-A1", "pi", "a1", use_svfit=use_svfit)
                pairs['aco_rho_a1_SVFIT'] = acoplanarity_calculator(pairs, "Decay-Plane", "Polarimetric-Vector-SVFIT-A1", "rho", "a1", use_svfit=use_svfit)
                pairs['aco_a1_pi_SVFIT'] = acoplanarity_calculator(pairs, "Polarimetric-Vector-SVFIT-A1", "Impact-Parameter", "a1", "pi", use_svfit=use_svfit)
                pairs['aco_a1_rho_SVFIT'] = acoplanarity_calculator(pairs, "Polarimetric-Vector-SVFIT-A1", "Decay-Plane", "a1", "rho", use_svfit=use_svfit)

            if use_fastmtt:
                pairs['aco_pi_a1_FASTMTT'] = acoplanarity_calculator(pairs, "Impact-Parameter", "Polarimetric-Vector-FASTMTT-A1", "pi", "a1", use_fastmtt=use_fastmtt)
                pairs['aco_rho_a1_FASTMTT'] = acoplanarity_calculator(pairs, "Decay-Plane", "Polarimetric-Vector-FASTMTT-A1", "rho", "a1", use_fastmtt=use_fastmtt)
                pairs['aco_a1_pi_FASTMTT'] = acoplanarity_calculator(pairs, "Polarimetric-Vector-FASTMTT-A1", "Impact-Parameter", "a1", "pi", use_fastmtt=use_fastmtt)
                pairs['aco_a1_rho_FASTMTT'] = acoplanarity_calculator(pairs, "Polarimetric-Vector-FASTMTT-A1", "Decay-Plane", "a1", "rho", use_fastmtt=use_fastmtt)

    return pairs
