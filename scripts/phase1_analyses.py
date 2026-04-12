#!/usr/bin/env python3
"""
Phase 1: All remaining reviewer-requested analyses.

1A. Structure-matched permutation null (compound-identity shuffle)
1B. Species-weighting sensitivity (2× mammalian)
1C. Temporal holdout (pre-2015 train, post-2015 test)
1D. Leave-mouse-out validation
1E. Negative-result simulation (30% synthetic negatives)
1F. Clustered bootstrap by PMID
"""
from __future__ import annotations

import csv
import json
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, fisher_exact

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "outputs" / "canonical"
FIGURES = CANONICAL / "figures"

CANONICAL_WEIGHTS = np.array([0.35, 0.20, 0.15, 0.10, 0.15, 0.05])

MAMMALIAN_TAXA = {"mammal"}
MOUSE_SPECIES = {"Mus musculus"}


def load_raw() -> list[dict]:
    rows = []
    with open(ROOT / "data" / "drugage_build5_2024-11-29.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                effect = float(row["avg_lifespan_change_percent"])
                rows.append({
                    "compound": row["compound_name"].strip(),
                    "species": row["species"].strip(),
                    "effect": effect,
                    "positive": 1 if effect > 0 else 0,
                    "pmid": row.get("pubmed_id", "").strip(),
                    "itp": row.get("ITP", "").strip(),
                    "gender": row.get("gender", "").strip(),
                })
            except (ValueError, KeyError):
                pass
    return rows


def load_species_map() -> dict[str, str]:
    """Load species → taxon mapping from YAML."""
    import yaml
    yaml_path = ROOT / "config" / "species_normalization.yaml"
    if yaml_path.exists():
        with open(yaml_path) as f:
            data = yaml.safe_load(f)
        return {entry.get("species", entry.get("drugage_name", "")): entry.get("taxon", "unknown")
                for entry in (data.get("species", []) if isinstance(data, dict) else data)
                if isinstance(entry, dict)}
    return {}


def compute_compound_scores(experiments: list[dict], weights: np.ndarray = CANONICAL_WEIGHTS,
                             species_weight_fn=None) -> dict[str, dict]:
    """Compute robustness scores from raw experiments."""
    by_compound = defaultdict(list)
    for e in experiments:
        by_compound[e["compound"]].append(e)

    results = {}
    for compound, exps in by_compound.items():
        n = len(exps)
        species = set(e["species"] for e in exps)
        taxa = set()  # would need taxon mapping
        pmids = set(e["pmid"] for e in exps if e["pmid"])
        n_pos = sum(e["positive"] for e in exps)
        effects = [e["effect"] for e in exps]

        median_eff = float(np.median(effects))
        trim_frac = max(1, int(0.1 * n))
        sorted_eff = sorted(effects)
        trimmed = sorted_eff[trim_frac:n - trim_frac] if n > 2 else sorted_eff
        trimmed_mean = float(np.mean(trimmed)) if trimmed else float(np.mean(effects))

        sign_cons = n_pos / n if n > 0 else 0
        agg_stab = 1.0 if median_eff > 0 and trimmed_mean > 0 else 0.0
        mag = min(max(trimmed_mean, 0), 50) / 50

        # LOSO
        if len(species) >= 2:
            loso_pos = 0
            for sp in species:
                subset = [e["effect"] for e in exps if e["species"] != sp]
                if subset:
                    sub_sorted = sorted(subset)
                    sub_trim = max(1, int(0.1 * len(subset)))
                    sub_trimmed = sub_sorted[sub_trim:len(subset)-sub_trim] if len(subset) > 2 else sub_sorted
                    if np.mean(sub_trimmed) > 0:
                        loso_pos += 1
            loso = loso_pos / len(species)
        else:
            loso = 0.0

        breadth_n = min(n, 6) / 6
        breadth_sp = min(len(species), 4) / 4
        breadth = (breadth_n + breadth_sp) / 2  # simplified without taxa

        score_components = np.array([breadth, sign_cons, loso, 0.0, agg_stab, mag])
        R = float(score_components @ weights)

        results[compound] = {
            "score": R, "n": n, "n_species": len(species), "n_pmids": len(pmids),
            "sign_consistency": sign_cons, "median_effect": median_eff,
            "trimmed_mean": trimmed_mean, "loso": loso, "breadth": breadth,
        }
    return results


def rank_from_scores(scores: dict[str, dict]) -> list[str]:
    return sorted(scores.keys(), key=lambda c: -scores[c]["score"])


# ============================================================
# 1A. Structure-matched permutation null
# ============================================================
def structure_matched_null(experiments: list[dict], n_perms: int = 1000) -> dict:
    """Permute compound identity within species, preserving per-compound sample sizes."""
    print("  Running structure-matched null...")
    canonical_scores = compute_compound_scores(experiments)
    canonical_ranking = rank_from_scores(canonical_scores)
    top_q = len(canonical_ranking) // 4

    # Count robust-like compounds in canonical
    canonical_robust = sum(1 for c, s in canonical_scores.items()
                          if s["n"] >= 3 and s["n_species"] >= 2
                          and s["sign_consistency"] >= 0.80 and s["loso"] == 1.0
                          and s["median_effect"] > 0 and s["trimmed_mean"] > 0)

    rng = np.random.default_rng(42)
    null_robust_counts = []

    by_species = defaultdict(list)
    for e in experiments:
        by_species[e["species"]].append(e)

    for perm_i in range(n_perms):
        shuffled = []
        for species, sp_exps in by_species.items():
            compound_labels = [e["compound"] for e in sp_exps]
            rng.shuffle(compound_labels)
            for e, new_compound in zip(sp_exps, compound_labels):
                shuffled.append({**e, "compound": new_compound})

        perm_scores = compute_compound_scores(shuffled)
        perm_robust = sum(1 for c, s in perm_scores.items()
                         if s["n"] >= 3 and s["n_species"] >= 2
                         and s["sign_consistency"] >= 0.80 and s["loso"] == 1.0
                         and s["median_effect"] > 0 and s["trimmed_mean"] > 0)
        null_robust_counts.append(perm_robust)

    null_arr = np.array(null_robust_counts)
    p_value = (np.sum(null_arr >= canonical_robust) + 1) / (n_perms + 1)

    return {
        "canonical_robust_count": canonical_robust,
        "null_mean": round(float(np.mean(null_arr)), 2),
        "null_std": round(float(np.std(null_arr)), 2),
        "null_z": round((canonical_robust - np.mean(null_arr)) / max(np.std(null_arr), 0.001), 2),
        "empirical_p": round(float(p_value), 6),
        "n_permutations": n_perms,
        "null_type": "compound_identity_shuffle_within_species",
    }


# ============================================================
# 1B. Species-weighting sensitivity
# ============================================================
def species_weight_sensitivity(experiments: list[dict]) -> dict:
    """Compare rankings under equal vs 2× mammalian weighting."""
    print("  Running species-weighting sensitivity...")
    species_taxon = load_species_map()

    # Equal weight (canonical)
    canonical_scores = compute_compound_scores(experiments)
    canonical_ranking = rank_from_scores(canonical_scores)

    # 2× mammalian: duplicate mammalian experiments
    mammalian_species = set()
    for sp, taxon in species_taxon.items():
        if taxon.lower() in MAMMALIAN_TAXA or "mus" in sp.lower() or "rattus" in sp.lower():
            mammalian_species.add(sp)

    weighted_exps = []
    for e in experiments:
        weighted_exps.append(e)
        if e["species"] in mammalian_species:
            weighted_exps.append(e)  # duplicate mammalian

    weighted_scores = compute_compound_scores(weighted_exps)
    weighted_ranking = rank_from_scores(weighted_scores)

    # Compare top-10 and top-20
    top10_canonical = set(canonical_ranking[:10])
    top10_weighted = set(weighted_ranking[:10])
    top20_canonical = set(canonical_ranking[:20])
    top20_weighted = set(weighted_ranking[:20])

    tau, tau_p = kendalltau(
        [canonical_ranking.index(c) if c in canonical_ranking else 9999 for c in weighted_ranking],
        list(range(len(weighted_ranking)))
    )

    return {
        "mammalian_species_identified": len(mammalian_species),
        "original_experiments": len(experiments),
        "weighted_experiments": len(weighted_exps),
        "top10_overlap": len(top10_canonical & top10_weighted),
        "top20_overlap": len(top20_canonical & top20_weighted),
        "kendall_tau": round(float(tau), 4) if np.isfinite(tau) else None,
        "top10_canonical_only": sorted(top10_canonical - top10_weighted),
        "top10_weighted_only": sorted(top10_weighted - top10_canonical),
    }


# ============================================================
# 1C. Temporal holdout
# ============================================================
def temporal_holdout(experiments: list[dict]) -> dict:
    """Train on pre-2016 PMIDs, test on post-2016 evidence accumulation."""
    print("  Running temporal holdout...")
    # PMIDs are roughly year-sortable by magnitude (newer = larger)
    # Use 25000000 as approximate cutoff for ~2015
    PMID_CUTOFF = 25000000

    pre = [e for e in experiments if e["pmid"] and int(e["pmid"]) < PMID_CUTOFF]
    post = [e for e in experiments if e["pmid"] and int(e["pmid"]) >= PMID_CUTOFF]
    no_pmid = [e for e in experiments if not e["pmid"]]

    print(f"    Pre-cutoff: {len(pre)}, Post-cutoff: {len(post)}, No PMID: {len(no_pmid)}")

    if len(pre) < 100 or len(post) < 100:
        return {"status": "insufficient_split", "pre": len(pre), "post": len(post)}

    pre_scores = compute_compound_scores(pre)
    pre_ranking = rank_from_scores(pre_scores)
    pre_top_q = set(pre_ranking[:len(pre_ranking) // 4])

    # For post-period, compute per-compound positive fraction
    post_by_compound = defaultdict(list)
    for e in post:
        post_by_compound[e["compound"]].append(e["positive"])

    # Compounds in pre top quartile: what's their post positive rate?
    pre_top_post_rates = []
    pre_bottom_post_rates = []
    for compound, post_outcomes in post_by_compound.items():
        rate = np.mean(post_outcomes)
        if compound in pre_top_q:
            pre_top_post_rates.append(rate)
        else:
            pre_bottom_post_rates.append(rate)

    from scipy.stats import mannwhitneyu
    if pre_top_post_rates and pre_bottom_post_rates:
        u, p = mannwhitneyu(pre_top_post_rates, pre_bottom_post_rates, alternative="greater")
    else:
        u, p = 0, 1.0

    return {
        "pmid_cutoff": PMID_CUTOFF,
        "pre_experiments": len(pre),
        "post_experiments": len(post),
        "pre_compounds": len(pre_scores),
        "post_compounds_with_data": len(post_by_compound),
        "pre_top_quartile_post_mean_positive_rate": round(float(np.mean(pre_top_post_rates)), 4) if pre_top_post_rates else None,
        "pre_bottom_post_mean_positive_rate": round(float(np.mean(pre_bottom_post_rates)), 4) if pre_bottom_post_rates else None,
        "mann_whitney_u": round(float(u), 1),
        "mann_whitney_p": round(float(p), 6),
        "pre_top_compounds_with_post_data": len(pre_top_post_rates),
        "pre_bottom_compounds_with_post_data": len(pre_bottom_post_rates),
    }


# ============================================================
# 1D. Leave-mouse-out validation
# ============================================================
def leave_mouse_out(experiments: list[dict]) -> dict:
    """Score without mouse data, validate against mouse-dominant outcomes."""
    print("  Running leave-mouse-out...")
    non_mouse = [e for e in experiments if e["species"] not in MOUSE_SPECIES]
    mouse_only = [e for e in experiments if e["species"] in MOUSE_SPECIES]

    print(f"    Non-mouse: {len(non_mouse)}, Mouse: {len(mouse_only)}")

    non_mouse_scores = compute_compound_scores(non_mouse)
    non_mouse_ranking = rank_from_scores(non_mouse_scores)
    top_q = set(non_mouse_ranking[:len(non_mouse_ranking) // 4])

    # For mouse-only compounds, compute positive fraction
    mouse_by_compound = defaultdict(list)
    for e in mouse_only:
        mouse_by_compound[e["compound"]].append(e["positive"])

    top_mouse_rates = []
    bottom_mouse_rates = []
    for compound, outcomes in mouse_by_compound.items():
        rate = np.mean(outcomes)
        if compound in top_q:
            top_mouse_rates.append(rate)
        else:
            bottom_mouse_rates.append(rate)

    from scipy.stats import mannwhitneyu
    if top_mouse_rates and bottom_mouse_rates:
        u, p = mannwhitneyu(top_mouse_rates, bottom_mouse_rates, alternative="greater")
    else:
        u, p = 0, 1.0

    # ITP validation (ITP is mouse-only)
    itp_data = json.loads((CANONICAL / "itp_overlap_analysis.json").read_text())
    itp_pos = set(d["compound"] for d in itp_data["compound_details"] if d["itp_direction"] == "positive")
    itp_neg = set(d["compound"] for d in itp_data["compound_details"] if d["itp_direction"] == "negative")

    itp_pos_in_top = len(itp_pos & top_q & set(non_mouse_scores.keys()))
    itp_pos_total = len(itp_pos & set(non_mouse_scores.keys()))
    itp_neg_in_top = len(itp_neg & top_q & set(non_mouse_scores.keys()))
    itp_neg_total = len(itp_neg & set(non_mouse_scores.keys()))

    return {
        "non_mouse_experiments": len(non_mouse),
        "mouse_experiments": len(mouse_only),
        "non_mouse_compounds": len(non_mouse_scores),
        "top_quartile_mouse_positive_rate": round(float(np.mean(top_mouse_rates)), 4) if top_mouse_rates else None,
        "bottom_mouse_positive_rate": round(float(np.mean(bottom_mouse_rates)), 4) if bottom_mouse_rates else None,
        "mann_whitney_p": round(float(p), 6),
        "itp_positive_in_top_quartile": f"{itp_pos_in_top}/{itp_pos_total}",
        "itp_negative_in_top_quartile": f"{itp_neg_in_top}/{itp_neg_total}",
    }


# ============================================================
# 1E. Negative-result simulation
# ============================================================
def negative_simulation(experiments: list[dict], neg_fraction: float = 0.30) -> dict:
    """Inject synthetic negative experiments and report tier survival."""
    print(f"  Running {neg_fraction:.0%} negative-result simulation...")
    rng = np.random.default_rng(123)

    canonical_scores = compute_compound_scores(experiments)
    canonical_robust = set(c for c, s in canonical_scores.items()
                          if s["n"] >= 3 and s["n_species"] >= 2
                          and s["sign_consistency"] >= 0.80 and s["loso"] == 1.0
                          and s["median_effect"] > 0 and s["trimmed_mean"] > 0)

    n_inject = int(len(experiments) * neg_fraction)
    compounds = list(set(e["compound"] for e in experiments))
    species = list(set(e["species"] for e in experiments))

    synthetic_negatives = []
    for _ in range(n_inject):
        synthetic_negatives.append({
            "compound": rng.choice(compounds),
            "species": rng.choice(species),
            "effect": float(rng.uniform(-30, -1)),
            "positive": 0,
            "pmid": "",
            "itp": "No",
            "gender": "Unknown",
        })

    augmented = experiments + synthetic_negatives
    aug_scores = compute_compound_scores(augmented)
    aug_robust = set(c for c, s in aug_scores.items()
                     if s["n"] >= 3 and s["n_species"] >= 2
                     and s["sign_consistency"] >= 0.80 and s["loso"] == 1.0
                     and s["median_effect"] > 0 and s["trimmed_mean"] > 0)

    survivors = canonical_robust & aug_robust
    lost = canonical_robust - aug_robust

    return {
        "negative_fraction": neg_fraction,
        "synthetic_negatives_injected": n_inject,
        "original_robust_count": len(canonical_robust),
        "augmented_robust_count": len(aug_robust),
        "survivors": len(survivors),
        "survival_rate": round(len(survivors) / max(len(canonical_robust), 1), 4),
        "lost_compounds": sorted(lost)[:10],
    }


# ============================================================
# 1F. Clustered bootstrap by PMID
# ============================================================
def clustered_pmid_bootstrap(experiments: list[dict], n_boot: int = 1000) -> dict:
    """Bootstrap by PMID cluster instead of by row."""
    print("  Running PMID-clustered bootstrap...")
    rng = np.random.default_rng(77)

    by_pmid = defaultdict(list)
    no_pmid = []
    for e in experiments:
        if e["pmid"]:
            by_pmid[e["pmid"]].append(e)
        else:
            no_pmid.append(e)

    pmids = list(by_pmid.keys())
    n_pmids = len(pmids)
    print(f"    {n_pmids} unique PMIDs, {len(no_pmid)} experiments without PMID")

    canonical_scores = compute_compound_scores(experiments)
    canonical_ranking = rank_from_scores(canonical_scores)

    boot_taus = []
    for _ in range(n_boot):
        sampled_pmids = rng.choice(pmids, size=n_pmids, replace=True)
        boot_exps = []
        for pmid in sampled_pmids:
            boot_exps.extend(by_pmid[pmid])
        boot_exps.extend(no_pmid)

        boot_scores = compute_compound_scores(boot_exps)
        boot_ranking = rank_from_scores(boot_scores)

        common = set(canonical_ranking) & set(boot_ranking)
        if len(common) > 10:
            canonical_order = [canonical_ranking.index(c) for c in common if c in canonical_ranking]
            boot_order = [boot_ranking.index(c) for c in common if c in boot_ranking]
            tau, _ = kendalltau(canonical_order, boot_order)
            if np.isfinite(tau):
                boot_taus.append(float(tau))

    taus = np.array(boot_taus)
    return {
        "n_bootstrap": n_boot,
        "n_pmids": n_pmids,
        "median_tau": round(float(np.median(taus)), 4) if len(taus) > 0 else None,
        "ci_low": round(float(np.percentile(taus, 2.5)), 4) if len(taus) > 0 else None,
        "ci_high": round(float(np.percentile(taus, 97.5)), 4) if len(taus) > 0 else None,
        "mean_tau": round(float(np.mean(taus)), 4) if len(taus) > 0 else None,
        "n_valid_bootstraps": len(taus),
    }


def main() -> int:
    experiments = load_raw()
    print(f"Loaded {len(experiments)} experiments\n")

    results = {}

    print("1A. Structure-matched permutation null...")
    results["structure_matched_null"] = structure_matched_null(experiments, n_perms=500)
    print(f"    Canonical robust: {results['structure_matched_null']['canonical_robust_count']}, "
          f"Null mean: {results['structure_matched_null']['null_mean']}, "
          f"z={results['structure_matched_null']['null_z']}, p={results['structure_matched_null']['empirical_p']}")

    print("\n1B. Species-weighting sensitivity...")
    results["species_weighting"] = species_weight_sensitivity(experiments)
    print(f"    Top-10 overlap: {results['species_weighting']['top10_overlap']}/10")
    print(f"    Top-20 overlap: {results['species_weighting']['top20_overlap']}/20")

    print("\n1C. Temporal holdout...")
    results["temporal_holdout"] = temporal_holdout(experiments)
    th = results["temporal_holdout"]
    if th.get("status") != "insufficient_split":
        print(f"    Pre top-Q post positive rate: {th['pre_top_quartile_post_mean_positive_rate']}")
        print(f"    Pre bottom post positive rate: {th['pre_bottom_post_mean_positive_rate']}")
        print(f"    Mann-Whitney p: {th['mann_whitney_p']}")

    print("\n1D. Leave-mouse-out...")
    results["leave_mouse_out"] = leave_mouse_out(experiments)
    lmo = results["leave_mouse_out"]
    print(f"    Top-Q mouse positive rate: {lmo['top_quartile_mouse_positive_rate']}")
    print(f"    Bottom mouse positive rate: {lmo['bottom_mouse_positive_rate']}")
    print(f"    ITP-pos in top-Q: {lmo['itp_positive_in_top_quartile']}")

    print("\n1E. Negative-result simulation (30%)...")
    results["negative_simulation"] = negative_simulation(experiments, 0.30)
    ns = results["negative_simulation"]
    print(f"    Original robust: {ns['original_robust_count']}, Augmented: {ns['augmented_robust_count']}, "
          f"Survival: {ns['survival_rate']:.0%}")

    print("\n1F. Clustered PMID bootstrap...")
    results["clustered_bootstrap"] = clustered_pmid_bootstrap(experiments, n_boot=500)
    cb = results["clustered_bootstrap"]
    print(f"    Median τ: {cb['median_tau']}, 95% CI: [{cb['ci_low']}, {cb['ci_high']}]")

    outpath = CANONICAL / "phase1_analyses.json"
    outpath.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nAll Phase 1 analyses saved to {outpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
