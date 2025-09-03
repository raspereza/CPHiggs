import correctionlib
import correctionlib.schemav2 as schema
import json
import uproot
import numpy as np
import logging
import re
from itertools import product

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

path_to_tau_json = "/vols/cms/ia2318/HiggsDNA/higgs_dna/systematics/ditau/JSONs/Tau_Trigger/2022_Summer22/tau_trigger.json.gz"
tau_json_name = "tauTriggerSF"
tau_evaluator = correctionlib.CorrectionSet.from_file(path_to_tau_json)[tau_json_name]


def get_doubletau_value(pt, VSjet, VSe, VSmu, corrtype, syst, values=None, trigtype=None, pt_bins=None, pt2_bins=None, pt2=None):
    dm = -1
    if corrtype == "sf":
        return tau_evaluator.evaluate(pt, dm, "ditau", VSjet, corrtype, syst)
    elif values is not None and trigtype is not None and pt_bins is not None and pt2_bins is not None:
        # You're asking for a value from a 2D flattened array
        vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
        pt_centers = get_bin_centers(pt_bins)
        pt2_centers = get_bin_centers(pt2_bins)

        # Get nearest bin indices (clamped)
        i = max(min(range(len(pt_centers)), key=lambda idx: abs(pt - pt_centers[idx])), 0)
        j = max(min(range(len(pt2_centers)), key=lambda idx: abs((pt2 or pt) - pt2_centers[idx])), 0)  # fall back to pt if pt2 missing

        idx = i * len(pt2_centers) + j
        return vals[idx]
    else:
        raise ValueError("Invalid inputs to get_doubletau_value")

def get_tau_pog_uncertainty(pt, VSjet, systype):
    sf_nom = tau_evaluator.evaluate(pt, -1, "ditau", VSjet, "sf", "nom")
    sf_up = tau_evaluator.evaluate(pt, -1, "ditau", VSjet, "sf", "up")
    sf_down = tau_evaluator.evaluate(pt, -1, "ditau", VSjet, "sf", "down")
    return sf_nom, sf_up, sf_down


# Define a mapping of WP combinations to SF values
sf_mapping = {
    "2022": {
        ("Medium", "VVLoose", "VLoose"): 0.93,
        ("Tight", "VVLoose", "VLoose"): 0.94,
        ("Medium", "VVLoose", "Tight"): 0.93,
        ("Tight", "VVLoose", "Tight"): 0.94,
        ("Medium", "Tight", "VLoose"): 0.91,
        ("Tight", "Tight", "VLoose"): 0.95,
        ("Medium", "Tight", "Tight"): 0.97,
        ("Tight", "Tight", "Tight"): 0.98,
    },
    "2022EE": {
        ("Medium", "VVLoose", "VLoose"): 0.93,
        ("Tight", "VVLoose", "VLoose"): 0.94,
        ("Medium", "VVLoose", "Tight"): 0.93,
        ("Tight", "VVLoose", "Tight"): 0.94,
        ("Medium", "Tight", "VLoose"): 0.91,
        ("Tight", "Tight", "VLoose"): 0.95,
        ("Medium", "Tight", "Tight"): 0.97,
        ("Tight", "Tight", "Tight"): 0.98,
    },
    "2023": {
        ("Medium", "VVLoose", "VLoose"): 0.92,
        ("Tight", "VVLoose", "VLoose"): 0.96,
        ("Medium", "VVLoose", "Tight"): 0.92,
        ("Tight", "VVLoose", "Tight"): 0.98,
        ("Medium", "Tight", "VLoose"): 0.92,
        ("Tight", "Tight", "VLoose"): 0.99,
        ("Medium", "Tight", "Tight"): 0.93,
        ("Tight", "Tight", "Tight"): 0.98,
    },
    "2023BPix": {
        ("Medium", "VVLoose", "VLoose"): 0.92,
        ("Tight", "VVLoose", "VLoose"): 0.96,
        ("Medium", "VVLoose", "Tight"): 0.92,
        ("Tight", "VVLoose", "Tight"): 0.98,
        ("Medium", "Tight", "VLoose"): 0.92,
        ("Tight", "Tight", "VLoose"): 0.99,
        ("Medium", "Tight", "Tight"): 0.93,
        ("Tight", "Tight", "Tight"): 0.98,
    }
}

uncerts_mapping = {
    "2022": {
        ("Medium", "VVLoose", "VLoose"): (0.09, 0.04, 0.15),
        ("Tight", "VVLoose", "VLoose"): (0.06, 0.04, 0.13),
        ("Medium", "VVLoose", "Tight"): (0.09, 0.04, 0.15),
        ("Tight", "VVLoose", "Tight"): (0.06, 0.04, 0.14),
        ("Medium", "Tight", "VLoose"): (0.08, 0.04, 0.18),
        ("Tight", "Tight", "VLoose"): (0.07, 0.04, 0.17),
        ("Medium", "Tight", "Tight"): (0.08, 0.05, 0.19),
        ("Tight", "Tight", "Tight"): (0.04, 0.06, 0.14),
    },
    "2022EE": {
        ("Medium", "VVLoose", "VLoose"): (0.09, 0.04, 0.15),
        ("Tight", "VVLoose", "VLoose"): (0.06, 0.04, 0.13),
        ("Medium", "VVLoose", "Tight"): (0.09, 0.04, 0.15),
        ("Tight", "VVLoose", "Tight"): (0.06, 0.04, 0.14),
        ("Medium", "Tight", "VLoose"): (0.08, 0.04, 0.18),
        ("Tight", "Tight", "VLoose"): (0.07, 0.04, 0.17),
        ("Medium", "Tight", "Tight"): (0.08, 0.05, 0.19),
        ("Tight", "Tight", "Tight"): (0.04, 0.06, 0.14),
    },
    "2023": {
        ("Medium", "VVLoose", "VLoose"): (0.09, 0.05, 0.15),
        ("Tight", "VVLoose", "VLoose"): (0.08, 0.05, 0.15),
        ("Medium", "VVLoose", "Tight"): (0.08, 0.05, 0.16),
        ("Tight", "VVLoose", "Tight"): (0.08, 0.05, 0.14),
        ("Medium", "Tight", "VLoose"): (0.07, 0.05, 0.19),
        ("Tight", "Tight", "VLoose"): (0.06, 0.05, 0.17),
        ("Medium", "Tight", "Tight"): (0.07, 0.05, 0.18),
        ("Tight", "Tight", "Tight"): (0.05, 0.05, 0.15),
    },
    "2023BPix": {
        ("Medium", "VVLoose", "VLoose"): (0.09, 0.05, 0.15),
        ("Tight", "VVLoose", "VLoose"): (0.08, 0.05, 0.15),
        ("Medium", "VVLoose", "Tight"): (0.08, 0.05, 0.16),
        ("Tight", "VVLoose", "Tight"): (0.08, 0.05, 0.14),
        ("Medium", "Tight", "VLoose"): (0.07, 0.05, 0.19),
        ("Tight", "Tight", "VLoose"): (0.06, 0.05, 0.17),
        ("Medium", "Tight", "Tight"): (0.07, 0.05, 0.18),
        ("Tight", "Tight", "Tight"): (0.05, 0.05, 0.15),
    },
}


def get_bin_centers(edges):
    return [(low + high) / 2 for low, high in zip(edges[:-1], edges[1:])]


# TODO: need to clean this function up, and make it more readable
def make_syst_category(vals_nom, corrtype, pt_bins, pt2_bins, trigtype, VSjet, VSe, VSmu, era):
    pt_centers = get_bin_centers(pt_bins)
    pt2_centers = get_bin_centers(pt2_bins)
    pt_grid = [(x, y) for x in pt_centers for y in pt2_centers]

    if corrtype == "eff_mc":  # mc_eff stat uncerts are subdominant
        vals_up = vals_nom
        vals_down = vals_nom

    elif corrtype == "sf":
        vals_up = []
        vals_down = []
        if trigtype in {"dt"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                sf_nom1, sf_up1, sf_down1 = get_tau_pog_uncertainty(pt, VSjet, "ditau")
                sf_nom2, sf_up2, sf_down2 = get_tau_pog_uncertainty(pt2, VSjet, "ditau")

                sf_nom = sf_nom1 * sf_nom2

                err_up = np.sqrt((sf_up1 - sf_nom1)**2 * sf_nom2**2 +
                                (sf_up2 - sf_nom2)**2 * sf_nom1**2)

                err_down = np.sqrt((sf_nom1 - sf_down1)**2 * sf_nom2**2 +
                                (sf_nom2 - sf_down2)**2 * sf_nom1**2)

                vals_up.append(sf_nom + err_up)
                vals_down.append(sf_nom - err_down)
        elif trigtype in {"dt_s1"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                sf_dt2, sf_up2, sf_down2 = get_tau_pog_uncertainty(pt2, VSjet, "ditau")
                sf_s1 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
                stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
                sigma_s1 = (stat**2 + (syst + syst_pt * (pt / 1000.0))**2)**0.5

                # Nominal SF
                sf_nom = sf_s1 * sf_dt2

                # Propagate uncertainty
                sigma_dt2_up = sf_up2 - sf_dt2
                sigma_dt2_down = sf_dt2 - sf_down2

                err_up = np.sqrt((sf_s1**2) * (sigma_dt2_up**2) + (sf_dt2**2) * (sigma_s1**2))
                err_down = np.sqrt((sf_s1**2) * (sigma_dt2_down**2) + (sf_dt2**2) * (sigma_s1**2))

                vals_up.append(sf_nom + err_up)
                vals_down.append(sf_nom - err_down)
        elif trigtype in {"dt_s2"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                sf_dt, sf_up, sf_down = get_tau_pog_uncertainty(pt, VSjet, "ditau")
                sf_s2 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
                stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
                sigma_s2 = (stat**2 + (syst + syst_pt * (pt / 1000.0))**2)**0.5

                # Nominal SF
                sf_nom = sf_s2 * sf_dt

                # Propagate uncertainty
                sigma_dt_up = sf_up - sf_dt
                sigma_dt_down = sf_dt - sf_down

                err_up = np.sqrt((sf_s2**2) * (sigma_dt_up**2) + (sf_dt**2) * (sigma_s2**2))
                err_down = np.sqrt((sf_s2**2) * (sigma_dt_down**2) + (sf_dt**2) * (sigma_s2**2))

                vals_up.append(sf_nom + err_up)
                vals_down.append(sf_nom - err_down)
        elif trigtype in {"s1_s2", "dt_s1s2"}:
            sf_s1 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
            sf_s2 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)

            # Nominal SF
            sf_nom = sf_s1 * sf_s2

            stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
            for pt, pt2 in pt_grid:
                # Get the uncertainties for each pt
                sigma_s1 = (stat**2 + (syst + syst_pt * (pt / 1000.0))**2)**0.5
                sigma_s2 = (stat**2 + (syst + syst_pt * (pt2 / 1000.0))**2)**0.5

                # Propagate uncertainty
                err = np.sqrt((sf_s2**2) * (sigma_s1**2) + (sf_s1**2) * (sigma_s2**2))

                vals_up.append(sf_nom + err)
                vals_down.append(sf_nom - err)

        elif trigtype in {"s1", "s2"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                this_pt = pt if trigtype == "s1" else pt2
                stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
                sigma = (stat**2 + (syst + syst_pt * (this_pt / 1000.0))**2)**0.5
                vals_up.append(v + sigma)
                vals_down.append(v - sigma)

    elif corrtype == "eff_data":
        vals_up = []
        vals_down = []
        if trigtype in {"dt"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                sf_nom1, sf_up1, sf_down1 = get_tau_pog_uncertainty(pt, VSjet, "ditau")
                sf_nom2, sf_up2, sf_down2 = get_tau_pog_uncertainty(pt2, VSjet, "ditau")

                sf_nom = sf_nom1 * sf_nom2

                err_up = np.sqrt((sf_up1 - sf_nom1)**2 * sf_nom2**2 +
                                (sf_up2 - sf_nom2)**2 * sf_nom1**2)

                err_down = np.sqrt((sf_nom1 - sf_down1)**2 * sf_nom2**2 +
                                (sf_nom2 - sf_down2)**2 * sf_nom1**2)

                eff_data = v

                vals_up.append(eff_data + eff_data * (err_up / sf_nom))
                vals_down.append(eff_data - eff_data * (err_down / sf_nom))
        elif trigtype in {"dt_s1"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                sf_dt2, sf_up2, sf_down2 = get_tau_pog_uncertainty(pt2, VSjet, "ditau")
                sf_s1 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
                stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
                sigma_s1 = (stat**2 + (syst + syst_pt * (pt / 1000.0))**2)**0.5

                # Nominal SF
                sf_nom = sf_s1 * sf_dt2

                # Propagate uncertainty
                sigma_dt2_up = sf_up2 - sf_dt2
                sigma_dt2_down = sf_dt2 - sf_down2

                err_up = np.sqrt((sf_s1**2) * (sigma_dt2_up**2) + (sf_dt2**2) * (sigma_s1**2))
                err_down = np.sqrt((sf_s1**2) * (sigma_dt2_down**2) + (sf_dt2**2) * (sigma_s1**2))

                eff_data = v

                vals_up.append(eff_data + eff_data * (err_up / sf_nom))
                vals_down.append(eff_data - eff_data * (err_down / sf_nom))
        elif trigtype in {"dt_s2"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                sf_dt, sf_up, sf_down = get_tau_pog_uncertainty(pt, VSjet, "ditau")
                sf_s2 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
                stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
                sigma_s2 = (stat**2 + (syst + syst_pt * (pt2 / 1000.0))**2)**0.5

                # Nominal SF
                sf_nom = sf_s2 * sf_dt

                # Propagate uncertainty
                sigma_dt_up = sf_up - sf_dt
                sigma_dt_down = sf_dt - sf_down

                err_up = np.sqrt((sf_s2**2) * (sigma_dt_up**2) + (sf_dt**2) * (sigma_s2**2))
                err_down = np.sqrt((sf_s2**2) * (sigma_dt_down**2) + (sf_dt**2) * (sigma_s2**2))

                eff_data = v

                vals_up.append(eff_data + eff_data * (err_up / sf_nom))
                vals_down.append(eff_data - eff_data * (err_down / sf_nom))
        elif trigtype in {"s1_s2", "dt_s1s2"}:
            sf_s1 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
            sf_s2 = sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)
            # Nominal SF
            sf_nom = sf_s1 * sf_s2
            stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                # Get the uncertainties for each pt
                sigma_s1 = (stat**2 + (syst + syst_pt * (pt / 1000.0))**2)**0.5
                sigma_s2 = (stat**2 + (syst + syst_pt * (pt2 / 1000.0))**2)**0.5

                # Propagate uncertainty
                err = np.sqrt((sf_s2**2) * (sigma_s1**2) + (sf_s1**2) * (sigma_s2**2))

                eff_data = v

                vals_up.append(eff_data + eff_data * (err / sf_nom))
                vals_down.append(eff_data - eff_data * (err / sf_nom))
        elif trigtype in {"s1", "s2"}:
            for (pt, pt2), v in zip(pt_grid, vals_nom):
                this_pt = pt if trigtype == "s1" else pt2
                stat, syst, syst_pt = uncerts_mapping[era][(VSjet, VSe, VSmu)]
                sigma = (stat**2 + (syst + syst_pt * (this_pt / 1000.0))**2)**0.5

                eff_data = v

                vals_up.append(eff_data + eff_data * (sigma / sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)))
                vals_down.append(eff_data - eff_data * (sigma / sf_mapping[era].get((VSjet, VSe, VSmu), 1.0)))

    return schema.Category(
        nodetype="category",
        input="syst",
        content=[
            schema.CategoryItem(key="nom", value=schema.MultiBinning(
                nodetype="multibinning",
                inputs=["pt", "pt_2"],
                edges=[pt_bins, pt2_bins],
                content=vals_nom,
                flow="clamp"
            )),
            schema.CategoryItem(key="up", value=schema.MultiBinning(
                nodetype="multibinning",
                inputs=["pt", "pt_2"],
                edges=[pt_bins, pt2_bins],
                content=vals_up,
                flow="clamp"
            )),
            schema.CategoryItem(key="down", value=schema.MultiBinning(
                nodetype="multibinning",
                inputs=["pt", "pt_2"],
                edges=[pt_bins, pt2_bins],
                content=vals_down,
                flow="clamp"
            )),
        ]
    )


# TODO: need to clean this function up, and make it more readable
def make_2d_correction(name, description, pt_bins, pt2_bins, values, trigtypes, vsJets, vsEles, vsMus, corrtypes, era, channel):
    trigtype_cats = []
    for trigtype in trigtypes:
        vsJet_cats = []
        for vsJet in vsJets:
            vsEle_cats = []
            for vsEle in vsEles:
                vsMu_cats = []
                for vsMu in vsMus:
                    corrtype_cats = []
                    for corrtype in corrtypes:
                        print(f"Processing {trigtype} for {vsJet}, {vsEle}, {vsMu} with corrtype {corrtype}")
                        if trigtype in {"dt"}:
                            vals = []
                            pt_centers = get_bin_centers(pt_bins)
                            pt2_centers = get_bin_centers(pt2_bins)

                            if corrtype == "sf":
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        sf_dt1 = get_doubletau_value(pt, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt_bins, pt2_bins, pt2)
                                        sf_dt2 = get_doubletau_value(pt2, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt2_bins, pt_bins, pt)
                                        vals.append(sf_dt1*sf_dt2)
                            elif corrtype == "eff_mc":
                                vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
                            elif corrtype == "eff_data":
                                mc_eff = values[trigtype][VSjet][VSe][VSmu][corrtype]
                                sfs = []
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        val_pt = get_doubletau_value(pt, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt_bins, pt2_bins, pt2)
                                        val_pt2 = get_doubletau_value(pt2, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt2_bins, pt_bins, pt)
                                        sfs.append(val_pt*val_pt2)
                                vals = [mc * sf for mc, sf in zip(mc_eff, sfs)]

                        elif trigtype in {"dt_s1"}:
                            vals = []
                            pt_centers = get_bin_centers(pt_bins)
                            pt2_centers = get_bin_centers(pt2_bins)

                            if corrtype == "sf":
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        sf_dt2 = get_doubletau_value(pt2, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt2_bins, pt_bins, pt)
                                        sf_s1 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        vals.append(sf_dt2*sf_s1)
                            elif corrtype == "eff_mc":
                                vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
                            elif corrtype == "eff_data":
                                mc_eff = values[trigtype][VSjet][VSe][VSmu][corrtype]
                                sfs = []
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        sf_dt2 = get_doubletau_value(pt2, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt2_bins, pt_bins, pt)
                                        sf_s1 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        sfs.append(sf_dt2*sf_s1)
                                vals = [mc * sf for mc, sf in zip(mc_eff, sfs)]

                        elif trigtype in {"dt_s2"}:
                            vals = []
                            pt_centers = get_bin_centers(pt_bins)
                            pt2_centers = get_bin_centers(pt2_bins)

                            if corrtype == "sf":
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        sf_dt1 = get_doubletau_value(pt, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt_bins, pt2_bins, pt2)
                                        sf_s2 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        vals.append(sf_dt1*sf_s2)
                            elif corrtype == "eff_mc":
                                vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
                            elif corrtype == "eff_data":
                                mc_eff = values[trigtype][VSjet][VSe][VSmu][corrtype]
                                sfs = []
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        sf_dt1 = get_doubletau_value(pt, VSjet, VSe, VSmu, "sf", "nom", values, trigtype, pt_bins, pt2_bins, pt2)
                                        sf_s2 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        sfs.append(sf_dt1*sf_s2)
                                vals = [mc * sf for mc, sf in zip(mc_eff, sfs)]

                        elif trigtype in {"s1_s2"}:
                            vals = []
                            pt_centers = get_bin_centers(pt_bins)
                            pt2_centers = get_bin_centers(pt2_bins)
                            if corrtype == "sf":
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        s1 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        s2 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        vals.append(s1*s2)
                            elif corrtype == "eff_mc":
                                vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
                            elif corrtype == "eff_data":
                                mc_eff = values[trigtype][VSjet][VSe][VSmu][corrtype]
                                sfs = []
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        s1 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        s2 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        sfs.append(s1*s2)
                                vals = [mc * sf for mc, sf in zip(mc_eff, sfs)]

                        elif trigtype in {"dt_s1s2"}:
                            vals = []
                            pt_centers = get_bin_centers(pt_bins)
                            pt2_centers = get_bin_centers(pt2_bins)

                            if corrtype == "sf":
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        s1 =  sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        s2 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        vals.append(s1*s2)
                            elif corrtype == "eff_mc":
                                vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
                            elif corrtype == "eff_data":
                                mc_eff = values[trigtype][VSjet][VSe][VSmu][corrtype]
                                sfs = []
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        s1 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        s2 = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        sfs.append(s1*s2)
                                vals = [mc * sf for mc, sf in zip(mc_eff, sfs)]
                        elif trigtype in {"s1", "s2"}:
                            vals = []
                            pt_centers = get_bin_centers(pt_bins)
                            pt2_centers = get_bin_centers(pt2_bins)

                            if corrtype == "sf":
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        s = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        vals.append(s)
                            elif corrtype == "eff_mc":
                                vals = values[trigtype][VSjet][VSe][VSmu][corrtype]
                            elif corrtype == "eff_data":
                                mc_eff = values[trigtype][VSjet][VSe][VSmu][corrtype]
                                sfs = []
                                for pt in pt_centers:
                                    for pt2 in pt2_centers:
                                        s = sf_mapping[era].get((vsJet, vsEle, vsMu), 1.0)
                                        sfs.append(s)
                                vals = [mc * sf for mc, sf in zip(mc_eff, sfs)]
                        syst_cat = make_syst_category(vals, corrtype, pt_bins, pt2_bins, trigtype, vsJet, vsEle, vsMu, era)
                        corrtype_cats.append(schema.CategoryItem(key=corrtype, value=syst_cat))
                    vsMu_cats.append(schema.CategoryItem(key=vsMu, value=schema.Category(
                        nodetype="category",
                        input="corrtype",
                        content=corrtype_cats
                    )))
                vsEle_cats.append(schema.CategoryItem(key=vsEle, value=schema.Category(
                    nodetype="category",
                    input="vsMu",
                    content=vsMu_cats
                )))
            vsJet_cats.append(schema.CategoryItem(key=vsJet, value=schema.Category(
                nodetype="category",
                input="vsEle",
                content=vsEle_cats
            )))
        trigtype_cats.append(schema.CategoryItem(key=trigtype, value=schema.Category(
            nodetype="category",
            input="vsJet",
            content=vsJet_cats
        )))

    return schema.Correction(
        name=name,
        description=description,
        version=1,
        inputs=[
            schema.Variable(name="pt", type="real", description="leading tau pt"),
            schema.Variable(name="pt_2", type="real", description="subleading tau pt"),
            schema.Variable(name="channel", type="string", description="Decay channel (tt, et, mt)"),
            schema.Variable(name="trigtype", type="string", description="Trigger type"),
            schema.Variable(name="vsJet", type="string", description="tau VSjet working point"),
            schema.Variable(name="vsEle", type="string", description="tau VSe working point"),
            schema.Variable(name="vsMu", type="string", description="tau VSmu working point"),
            schema.Variable(name="corrtype", type="string", description="eff_data, eff_mc or sf"),
            schema.Variable(name="syst", type="string", description="nom, up, down")
        ],
        output=schema.Variable(name="weight", type="real"),
        data=schema.Category(
            nodetype="category",
            input="channel",
            content=[
                schema.CategoryItem(key=channel, value=schema.Category(
                    nodetype="category",
                    input="trigtype",
                    content=trigtype_cats
                ))
            ]
        )
    )


VSjet_wp = ["Medium", "Tight"]
VSe_wp = ["VVLoose", "Tight"]
VSmu_wp = ["VLoose", "Tight"]

for era in ["2022", "2022EE", "2023", "2023BPix"]:
    corrections = []
    sf_mapping_era = sf_mapping[era]
    uncerts_mapping_era = uncerts_mapping[era]
    for VSjet, VSe, VSmu in product(VSjet_wp, VSe_wp, VSmu_wp):
        file = uproot.open(f"/vols/cms/ia2318/TIDAL/MSSMEfficiences_new_wp/{era}/efficiency_2D_mssm_signals_tt_{VSjet}_{VSe}_{VSmu}.root")
        for hist_name in [
            "eff_mssm_signals_tt_trg_doubletau",
            "eff_mssm_signals_tt_trg_singletauandpt_1>190",
            "eff_mssm_signals_tt_trg_singletau_2andpt_2>190",
            "eff_mssm_signals_tt_trg_doubletauandtrg_singletauandpt_1>190",
            "eff_mssm_signals_tt_trg_doubletauandtrg_singletau_2andpt_2>190",
            "eff_mssm_signals_tt_trg_singletauandtrg_singletau_2andpt_1>190andpt_2>190",
            "eff_mssm_signals_tt_trg_doubletauandtrg_singletauandtrg_singletau_2andpt_1>190andpt_2>190"

        ]:
            hist = file[hist_name]
            pt_bins = hist.axis("x").edges()
            pt2_bins = hist.axis("y").edges()
            eff_mc = [1.0 if value == 0.0 else value for value in hist.values(flow=False).flatten().tolist()]
            sf_vals = [1.0 if value == 0.0 else sf_mapping_era.get((VSjet, VSe, VSmu), 1.0) for value in hist.values(flow=False).flatten().tolist()]
            eff_data = [sf * mc for sf, mc in zip(sf_vals, eff_mc)]

            if hist_name == "eff_mssm_signals_tt_trg_singletauandpt_1>190":
                trigtype = "s1"
            elif hist_name == "eff_mssm_signals_tt_trg_singletau_2andpt_2>190":
                trigtype = "s2"
            elif hist_name == "eff_mssm_signals_tt_trg_doubletauandtrg_singletauandpt_1>190":
                trigtype = "dt_s1"
            elif hist_name == "eff_mssm_signals_tt_trg_doubletauandtrg_singletau_2andpt_2>190":
                trigtype = "dt_s2"
            elif hist_name == "eff_mssm_signals_tt_trg_singletauandtrg_singletau_2andpt_1>190andpt_2>190":
                trigtype = "s1_s2"
            elif hist_name == "eff_mssm_signals_tt_trg_doubletauandtrg_singletauandtrg_singletau_2andpt_1>190andpt_2>190":
                trigtype = "dt_s1s2"
            elif hist_name == "eff_mssm_signals_tt_trg_doubletau":
                trigtype = "dt"

            values = {
                trigtype: {
                    VSjet: {
                        VSe: {
                            VSmu: {
                                "eff_data": eff_data,
                                "eff_mc": eff_mc,
                                "sf": sf_vals
                            }
                        }
                    }
                }
            }

            correction = make_2d_correction(
                name=f"tauTriggerSF_tt_{trigtype}_{VSjet}_{VSe}_{VSmu}",
                description=f"Tau trigger SF for tt - {trigtype} with {VSjet}, {VSe}, {VSmu} working points for {era}",
                pt_bins=pt_bins,
                pt2_bins=pt2_bins,
                values=values,
                trigtypes=[trigtype],
                vsJets=[VSjet],
                vsEles=[VSe],
                vsMus=[VSmu],
                corrtypes=["eff_data", "eff_mc", "sf"],
                era=era,
                channel="tt"
            )
            corrections.append(correction)

    for VSjet, VSe, VSmu in product(VSjet_wp, VSe_wp, VSmu_wp):
        file = uproot.open(f"/vols/cms/ia2318/TIDAL/MSSMEfficiences_new_wp/{era}/efficiency_2D_mssm_signals_mt_{VSjet}_{VSe}_{VSmu}.root")
        hist = file["eff_mssm_signals_mt_trg_singletauandpt_2>190"]
        pt_bins = hist.axis("x").edges()
        pt2_bins = hist.axis("y").edges()
        eff_mc = [1.0 if v == 0.0 else v for v in hist.values(flow=False).flatten().tolist()]
        sf_vals = [1.0 if v == 0.0 else sf_mapping[era].get((VSjet, VSe, VSmu), 1.0) for v in eff_mc]
        eff_data = [sf * mc for sf, mc in zip(sf_vals, eff_mc)]

        trigtype = "s2"  # I think this is the same treatment as in tt
        values = {
            trigtype: {
                VSjet: {
                    VSe: {
                        VSmu: {
                            "eff_data": eff_data,
                            "eff_mc": eff_mc,
                            "sf": sf_vals
                        }
                    }
                }
            }
        }

        correction = make_2d_correction(
            name=f"tauTriggerSF_mt_{trigtype}_{VSjet}_{VSe}_{VSmu}",
            description=f"Tau trigger SF for mt - {trigtype} with {VSjet}, {VSe}, {VSmu} working points for {era}",
            pt_bins=pt_bins,
            pt2_bins=pt2_bins,
            values=values,
            trigtypes=[trigtype],
            vsJets=[VSjet],
            vsEles=[VSe],
            vsMus=[VSmu],
            corrtypes=["eff_data", "eff_mc", "sf"],
            era=era,
            channel="mt"
        )
        corrections.append(correction)

    for VSjet, VSe, VSmu in product(VSjet_wp, VSe_wp, VSmu_wp):
        file = uproot.open(f"/vols/cms/ia2318/TIDAL/MSSMEfficiences_new_wp/{era}/efficiency_2D_mssm_signals_et_{VSjet}_{VSe}_{VSmu}.root")
        hist = file["eff_mssm_signals_et_trg_singletauandpt_2>190"]
        pt_bins = hist.axis("x").edges()
        pt2_bins = hist.axis("y").edges()
        eff_mc = [1.0 if v == 0.0 else v for v in hist.values(flow=False).flatten().tolist()]
        sf_vals = [1.0 if v == 0.0 else sf_mapping[era].get((VSjet, VSe, VSmu), 1.0) for v in eff_mc]
        eff_data = [sf * mc for sf, mc in zip(sf_vals, eff_mc)]

        trigtype = "s2"  # I think this is the same treatment as in tt
        values = {
            trigtype: {
                VSjet: {
                    VSe: {
                        VSmu: {
                            "eff_data": eff_data,
                            "eff_mc": eff_mc,
                            "sf": sf_vals
                        }
                    }
                }
            }
        }

        correction = make_2d_correction(
            name=f"tauTriggerSF_et_{trigtype}_{VSjet}_{VSe}_{VSmu}",
            description=f"Tau trigger SF for et - {trigtype} with {VSjet}, {VSe}, {VSmu} working points for {era}",
            pt_bins=pt_bins,
            pt2_bins=pt2_bins,
            values=values,
            trigtypes=[trigtype],
            vsJets=[VSjet],
            vsEles=[VSe],
            vsMus=[VSmu],
            corrtypes=["eff_data", "eff_mc", "sf"],
            era=era,
            channel="et"
        )
        corrections.append(correction)

    output = schema.CorrectionSet(schema_version=2, corrections=corrections)

    # def format_array_multiline(data, max_items_per_line=10, indent_level=32):
    #     indent = " " * indent_level
    #     lines = []
    #     for i in range(0, len(data), max_items_per_line):
    #         chunk = data[i:i + max_items_per_line]
    #         line = indent + ", ".join(f"{x:.8g}" for x in chunk) + ","
    #         lines.append(line)
    #     if lines:
    #         lines[-1] = lines[-1].rstrip(',')
    #     return "[\n" + "\n".join(lines) + "\n" + (" " * (indent_level)) + "]"

    # with open("tau_trigger_sf.json", "w") as f:
    #     raw = json.dumps(output.model_dump(), indent=2)

    #     def replace_array(match):
    #         arr = json.loads(match.group(0))
    #         return format_array_multiline(arr)

    #     pattern = re.compile(r"\[\s*(?:-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\s*,\s*)+-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\s*\]")
    #     pretty = pattern.sub(replace_array, raw)
    #     f.write(pretty)

    with open(f"tau_trigger_sf_{era}.json", "w") as f:
        json.dump(output.model_dump(), f, indent=2)

    logging.info(f"Correction saved to tau_trigger_sf_{era}.json")

    # gzip the file
    import gzip
    import shutil
    with open(f"tau_trigger_sf_{era}.json", 'rb') as f_in:
        with gzip.open(f'tau_trigger_sf_{era}.json.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    logging.info(f"Correction gzipped to tau_trigger_sf_{era}.json.gz")
