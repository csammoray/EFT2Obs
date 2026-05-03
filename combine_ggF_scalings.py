import json
from collections import OrderedDict
import argparse
import numpy as np
import os

one_to_one_conversion = {
    "cll1221": "cll1",
    "cpdc": "chdd",
    "cdp": "chbox",
    "c3pl1": "chl3",
    "c3pl2": "chl3",
    "ctp": "cthre",
    "ctg": "ctgre",
    "cpg": "chg",
}

def checkCompatibility(json1, json2):
    if json1.get("nbins") != json2.get("nbins"):
        raise ValueError(f"nbins mismatch: {json1.get('nbins')} vs {json2.get('nbins')}")
    if json1.get("bin_edges") != json2.get("bin_edges"):
        raise ValueError(f"bin_edges mismatch between {json1} and {json2}")
    # bin_labels may be empty; if both non-empty require equality
    if json1.get("bin_labels") and json2.get("bin_labels") and json1.get("bin_labels") != json2.get("bin_labels"):
        raise ValueError(f"bin_labels mismatch between {json1} and {json2}")

def formatCoeffs(coeffs_list):
    """
    coeffs_list is a list like ["ctg"] or ["ctg","ctp"] or ["ctg","ctg"].
    """
    if len(coeffs_list) == 1:
        return f"a_{coeffs_list[0]}"
    if len(coeffs_list) == 2:
        wc1, wc2 = coeffs_list
        if wc1 == wc2:
            return f"b_{wc1}_2"
        # Avoid duplicates (c1_c2 vs c2_c1)
        a, b = sorted([wc1, wc2])
        return f"b_{a}_{b}"

def buildCoeffDict(input_dict):
    """
    Converts input_dict["terms"] list into dict:
      key -> (coeffs[np.array], uncs[np.array])
    """
    nbins = input_dict["nbins"]
    out = {}
    for entry in input_dict["terms"]:
        coeffs_list, coeffs, uncs = entry
        key = formatCoeffs(coeffs_list)
        coeffs = np.asarray(coeffs, dtype=float)
        uncs = np.asarray(uncs, dtype=float)
        if len(coeffs) != nbins or len(uncs) != nbins:
            raise ValueError(f"Term {coeffs_list} has wrong nbins: {len(coeffs)}/{len(uncs)} vs {nbins}")
        out[key] = (coeffs, uncs)
    return out

def combineContributions(loop, tree, tree_loop_2, tree_loop_4):
    # Consistency check
    checkCompatibility(loop, tree)
    checkCompatibility(loop, tree_loop_2)
    checkCompatibility(loop, tree_loop_4)

    nbins = loop["nbins"]
    dict_per_contribution = [
        buildCoeffDict(loop),
        buildCoeffDict(tree),
        buildCoeffDict(tree_loop_2),
        buildCoeffDict(tree_loop_4),
    ]

    all_keys = set()
    for term in dict_per_contribution:
        all_keys |= set(term.keys())

    combined = {}
    for key in all_keys:
        coeff_sum = np.zeros(nbins, dtype=float)
        var_sum = np.zeros(nbins, dtype=float)

        for coeff_dict in dict_per_contribution:
            if key in coeff_dict:
                coeff, unc = coeff_dict[key]
                coeff_sum += coeff
                var_sum += unc * unc

        combined[key] = (coeff_sum, np.sqrt(var_sum))

    meta = dict(loop)  # keep nbins/bin_edges/sm_vals/etc
    return meta, combined


def convert_SMEFTatNLO_to_SMEFTsim(coeff_dict):
    """
    Convert keys from SMEFTatNLO naming to SMEFTsim naming and rescale ctg terms.
    Input/Output: coeff_dict[key] = (coeffs[np.array], uncs[np.array])
    """
    a_s = 0.1181
    g_s = 2 * np.sqrt(a_s * np.pi)

    converted = {}
    for key, (coeffs, uncs) in coeff_dict.items():
        new_key = key
        for old, new in one_to_one_conversion.items():
            new_key = new_key.replace(old, new)

        # Special case because degeneracy of chl3 and c3pl1/2 in SMEFTsim naming
        if new_key == "b_chl3_chl3":
            new_key = "b_chl3_2"

        # rescale ctg -> ctgre conventions
        if "ctgre" in new_key:
            if new_key.startswith("a_"):
                coeffs = coeffs / (-g_s)
                uncs = uncs / (g_s)
            else:
                if new_key.endswith("_2"):
                    coeffs = coeffs / ((-g_s) ** 2)
                    uncs = uncs / (g_s ** 2)
                else:
                    coeffs = coeffs / (-g_s)
                    uncs = uncs / (g_s)

        if new_key in converted:
            c0, u0 = converted[new_key]
            converted[new_key] = (c0 + coeffs, np.sqrt(u0 * u0 + uncs * uncs))
        else:
            converted[new_key] = (coeffs, uncs)

    return converted

def format_to_common_json(meta, coeff_dict):
    central = {}
    u_mc = {}

    for key, (coeffs, uncs) in coeff_dict.items():
        central[key] = list(map(float, coeffs))
        u_mc[key] = list(map(float, uncs))

    # infer list of coefficients appearing
    coeffs_set = set()
    for key in central.keys():
        if key.startswith("a_"):
            coeffs_set.add(key[2:])
        elif key.startswith("b_"):
            body = key[2:]
            if body.endswith("_2"):
                coeffs_set.add(body[:-2])
            else:
                wc1, wc2 = body.split("_", 1)
                coeffs_set.add(wc1)
                coeffs_set.add(wc2)

    return {
        "data": {
            "central": central,
            "u_MC": u_mc,
            "sm_xs": meta.get("sm_vals", []),
        },
        "metadata": {
            "coefficients": sorted(coeffs_set),
            "observable_names": meta.get("bin_edges", []),
        },
    }

def main(args):
    input_dir = args.input_dir
    output_dir = args.output_dir
    suffix = f"_{args.obs}" if args.obs else ""

    def load(name):
        with open(os.path.join(input_dir, name), "r") as f:
            return json.load(f)

    loop = load(f"ggH_SMEFTatNLO_loop{suffix}.json")
    tree = load(f"ggH_SMEFTatNLO_tree{suffix}.json")
    tree_loop_2 = load(f"ggH_SMEFTatNLO_tree_loop_2{suffix}.json")
    tree_loop_4 = load(f"ggH_SMEFTatNLO_tree_loop_4{suffix}.json")

    new_terms = []
    for wcs, coeffs, uncs in tree_loop_4["terms"]:
        if wcs == ["cpg", "cpg"]:
            coeffs = [0.0] * tree_loop_4["nbins"]
            uncs = [0.0] * tree_loop_4["nbins"]
        new_terms.append([wcs, coeffs, uncs])
    tree_loop_4["terms"] = new_terms

    meta, combined_map = combineContributions(loop, tree, tree_loop_2, tree_loop_4)
    converted_map = convert_SMEFTatNLO_to_SMEFTsim(combined_map)

    out_json = format_to_common_json(meta, converted_map)

    output = os.path.join(output_dir, f"ggH_SMEFTatNLO_combined{suffix}.json")
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as f:
        json.dump(out_json, f, indent=4)
    
    print(f"Combined results written to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", "-i", type=str, default="/eos/user/c/csammora/ggH_results/500_events/HIG-25-015/equations/deltaphijj/")
    parser.add_argument("--output-dir", "-o", type=str, default="/eos/user/c/csammora/ggH_results/500_events/HIG-25-015/equations/deltaphijj/")
    parser.add_argument("--obs", type=str, required=True, default="", help="Observable used for differential scaling (e.g. m_z2, pt_h, ...).",
    )
    args = parser.parse_args()
    main(args)