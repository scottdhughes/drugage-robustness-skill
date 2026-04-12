#!/usr/bin/env python3
"""
Ceiling-level analyses for the DrugAge robustness paper.

Generates all Tier 1 + Tier 2 artifacts:
- Scatter plot: n_experiments vs R (colored by tier)
- Residualized robustness score (R_adjusted = R - predicted from log(n_experiments))
- Weight sensitivity sweep (100 Dirichlet samples, Kendall tau)
- Fisher's exact test on ITP 2x2 table
- Experiment counts for ITP-positive vs ITP-negative
- Naive baseline comparison (median, trimmed mean, experiment count) with ITP concordance
- Component score correlation matrix (6x6)
- Population statistics (experiments-per-compound histogram)
- Species-weighting sensitivity (equal vs inverse-frequency)
"""
from __future__ import annotations

import csv
import json
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact, kendalltau, spearmanr, linregress

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "outputs" / "canonical"
FIGURES = CANONICAL / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

WEIGHT_NAMES = ["breadth", "sign_consistency", "loso", "loto", "aggregation_stability", "magnitude"]
CANONICAL_WEIGHTS = np.array([0.35, 0.20, 0.15, 0.10, 0.15, 0.05])


def load_rankings() -> list[dict]:
    rows = []
    with open(CANONICAL / "robustness_rankings.csv") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def load_profiles() -> list[dict]:
    rows = []
    with open(CANONICAL / "compound_evidence_profiles.csv") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def load_itp() -> dict:
    return json.loads((CANONICAL / "itp_overlap_analysis.json").read_text())


def get_component_scores(profiles: list[dict]) -> dict[str, np.ndarray]:
    compounds = [p.get("compound_name", p.get("compound", "")).strip() for p in profiles]
    breadth = np.array([float(p.get("breadth_score", 0)) for p in profiles])
    sign = np.array([float(p.get("sign_consistency", 0)) for p in profiles])
    loso = np.array([float(p.get("leave_one_species_out_stability", 0)) for p in profiles])
    loto = np.array([float(p.get("leave_one_taxon_out_stability", 0)) for p in profiles])
    agg = np.array([float(p.get("aggregation_stability", 0)) for p in profiles])
    trimmed = np.array([float(p.get("trimmed_mean_effect", 0)) for p in profiles])
    mag = np.clip(np.clip(trimmed, 0, None), 0, 50) / 50.0
    n_exp = np.array([int(p.get("num_experiments", 1)) for p in profiles])
    n_species = np.array([int(p.get("num_species", 1)) for p in profiles])
    return {
        "compounds": compounds,
        "breadth": breadth, "sign_consistency": sign,
        "loso": loso, "loto": loto,
        "aggregation_stability": agg, "magnitude": mag,
        "n_experiments": n_exp, "n_species": n_species,
    }


def compute_robustness(components: dict[str, np.ndarray], weights: np.ndarray) -> np.ndarray:
    score_matrix = np.column_stack([
        components["breadth"], components["sign_consistency"],
        components["loso"], components["loto"],
        components["aggregation_stability"], components["magnitude"],
    ])
    return score_matrix @ weights


# ============================================================
# 1. Scatter plot: n_experiments vs R
# ============================================================
def scatter_n_exp_vs_r(rankings: list[dict], profiles: list[dict]) -> dict:
    comp = get_component_scores(profiles)
    R = compute_robustness(comp, CANONICAL_WEIGHTS)

    tier_map = {}
    for r in rankings:
        name = r.get("compound", r.get("compound_name", "")).strip()
        tier_map[name] = r.get("tier", r.get("evidence_tier", "")).strip()

    itp_data = load_itp()
    itp_set = set(d["compound"] for d in itp_data["compound_details"])

    tier_colors = {"robust": "#2ecc71", "promising": "#3498db", "thin_evidence": "#f39c12", "conflicted": "#e74c3c"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, name in enumerate(comp["compounds"]):
        tier = tier_map.get(name, "conflicted")
        color = tier_colors.get(tier, "#999999")
        marker = "D" if name in itp_set else "o"
        size = 40 if name in itp_set else 12
        alpha = 0.9 if name in itp_set else 0.4
        ax.scatter(comp["n_experiments"][i], R[i], c=color, marker=marker, s=size, alpha=alpha, edgecolors="none")

    rho, pval = spearmanr(comp["n_experiments"], R)
    ax.set_xlabel("Number of experiments in DrugAge")
    ax.set_ylabel("Robustness score R")
    ax.set_title(f"Experiment count vs robustness score (Spearman ρ={rho:.3f}, p={pval:.1e})")
    ax.set_xscale("log")

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8, label=t)
        for t, c in tier_colors.items()
    ] + [Line2D([0], [0], marker="D", color="w", markerfacecolor="gray", markersize=8, label="ITP-tested")]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURES / "scatter_n_experiments_vs_robustness.png", dpi=150)
    plt.close(fig)

    result = {"spearman_rho": round(float(rho), 4), "spearman_p": float(pval), "n_compounds": len(R)}
    return result


# ============================================================
# 2. Residualized robustness score
# ============================================================
def residualized_score(rankings: list[dict], profiles: list[dict]) -> dict:
    comp = get_component_scores(profiles)
    R = compute_robustness(comp, CANONICAL_WEIGHTS)
    log_n = np.log1p(comp["n_experiments"])

    slope, intercept, r_value, p_value, std_err = linregress(log_n, R)
    R_predicted = slope * log_n + intercept
    R_adjusted = R - R_predicted

    rank_original = np.argsort(-R)
    rank_adjusted = np.argsort(-R_adjusted)

    tau, tau_p = kendalltau(rank_original, rank_adjusted)

    tier_map = {}
    for r in rankings:
        name = r.get("compound", r.get("compound_name", "")).strip()
        tier_map[name] = r.get("tier", r.get("evidence_tier", "")).strip()

    itp_data = load_itp()
    itp_details = {d["compound"]: d for d in itp_data["compound_details"]}

    top20_original = [comp["compounds"][i] for i in rank_original[:20]]
    top20_adjusted = [comp["compounds"][i] for i in rank_adjusted[:20]]

    itp_pos_in_robust_adjusted = 0
    itp_neg_in_conflicted_adjusted = 0
    for d in itp_data["compound_details"]:
        idx = comp["compounds"].index(d["compound"]) if d["compound"] in comp["compounds"] else -1
        if idx < 0:
            continue
        adj_rank = int(np.where(rank_adjusted == idx)[0][0]) + 1
        if d["itp_direction"] == "positive" and adj_rank <= len(rankings) // 4:
            itp_pos_in_robust_adjusted += 1
        if d["itp_direction"] == "negative":
            orig_tier = tier_map.get(d["compound"], "")
            if orig_tier == "conflicted":
                itp_neg_in_conflicted_adjusted += 1

    result = {
        "regression_r_squared": round(r_value**2, 4),
        "regression_slope": round(slope, 4),
        "kendall_tau_original_vs_adjusted": round(float(tau), 4),
        "kendall_tau_p": float(tau_p),
        "top20_original": top20_original,
        "top20_adjusted": top20_adjusted,
        "top20_overlap_count": len(set(top20_original) & set(top20_adjusted)),
    }
    return result


# ============================================================
# 3. Weight sensitivity sweep
# ============================================================
def weight_sensitivity(profiles: list[dict], n_samples: int = 200) -> dict:
    comp = get_component_scores(profiles)
    R_canonical = compute_robustness(comp, CANONICAL_WEIGHTS)
    rank_canonical = np.argsort(-R_canonical)

    rng = np.random.default_rng(42)
    alpha_prior = CANONICAL_WEIGHTS * 10
    taus = []
    for _ in range(n_samples):
        w = rng.dirichlet(alpha_prior)
        R_perm = compute_robustness(comp, w)
        rank_perm = np.argsort(-R_perm)
        tau, _ = kendalltau(rank_canonical, rank_perm)
        taus.append(float(tau))

    taus = np.array(taus)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(taus, bins=30, edgecolor="black", alpha=0.7)
    ax.axvline(np.median(taus), color="red", linestyle="--", label=f"median τ={np.median(taus):.3f}")
    ax.set_xlabel("Kendall τ vs canonical ranking")
    ax.set_ylabel("Count")
    ax.set_title(f"Weight sensitivity ({n_samples} Dirichlet samples)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "weight_sensitivity_tau.png", dpi=150)
    plt.close(fig)

    result = {
        "n_samples": n_samples,
        "median_tau": round(float(np.median(taus)), 4),
        "mean_tau": round(float(np.mean(taus)), 4),
        "min_tau": round(float(np.min(taus)), 4),
        "p5_tau": round(float(np.percentile(taus, 5)), 4),
        "p95_tau": round(float(np.percentile(taus, 95)), 4),
    }
    return result


# ============================================================
# 4. Fisher's exact test on ITP 2x2
# ============================================================
def fisher_itp(rankings: list[dict]) -> dict:
    itp_data = load_itp()
    tier_map = {}
    for r in rankings:
        name = r.get("compound", r.get("compound_name", "")).strip()
        tier_map[name] = r.get("tier", r.get("evidence_tier", "")).strip()

    pos_robust = 0
    pos_not_robust = 0
    neg_robust = 0
    neg_not_robust = 0
    for d in itp_data["compound_details"]:
        tier = d.get("tier", "")
        direction = d.get("itp_direction", "")
        is_robust = tier == "robust"
        if direction == "positive":
            if is_robust:
                pos_robust += 1
            else:
                pos_not_robust += 1
        elif direction == "negative":
            if is_robust:
                neg_robust += 1
            else:
                neg_not_robust += 1

    table = [[pos_robust, pos_not_robust], [neg_robust, neg_not_robust]]
    odds_ratio, p_value = fisher_exact(table, alternative="greater")

    return {
        "table": {"pos_robust": pos_robust, "pos_not_robust": pos_not_robust,
                   "neg_robust": neg_robust, "neg_not_robust": neg_not_robust},
        "odds_ratio": round(float(odds_ratio), 2) if np.isfinite(odds_ratio) else "inf",
        "fisher_p_one_sided": round(float(p_value), 6),
    }


# ============================================================
# 5. Experiment counts for ITP-pos vs ITP-neg
# ============================================================
def itp_experiment_counts(profiles: list[dict]) -> dict:
    itp_data = load_itp()
    profile_map = {}
    for p in profiles:
        name = p.get("compound", p.get("compound_name", "")).strip()
        profile_map[name] = int(p.get("num_experiments", 1))

    pos_counts = [profile_map.get(d["compound"], 0) for d in itp_data["compound_details"] if d["itp_direction"] == "positive"]
    neg_counts = [profile_map.get(d["compound"], 0) for d in itp_data["compound_details"] if d["itp_direction"] == "negative"]
    mixed_counts = [profile_map.get(d["compound"], 0) for d in itp_data["compound_details"] if d["itp_direction"] == "mixed"]

    return {
        "itp_positive_mean_experiments": round(np.mean(pos_counts), 1) if pos_counts else None,
        "itp_positive_median_experiments": round(np.median(pos_counts), 1) if pos_counts else None,
        "itp_negative_mean_experiments": round(np.mean(neg_counts), 1) if neg_counts else None,
        "itp_negative_median_experiments": round(np.median(neg_counts), 1) if neg_counts else None,
        "itp_mixed_mean_experiments": round(np.mean(mixed_counts), 1) if mixed_counts else None,
    }


# ============================================================
# 6. Naive baseline comparison
# ============================================================
def naive_baselines(rankings: list[dict], profiles: list[dict]) -> dict:
    comp = get_component_scores(profiles)
    R = compute_robustness(comp, CANONICAL_WEIGHTS)
    rank_R = np.argsort(-R)

    median_effects = np.array([float(p.get("median_effect", 0)) for p in profiles])
    trimmed_effects = np.array([float(p.get("trimmed_mean_effect", 0)) for p in profiles])
    n_exp = comp["n_experiments"].astype(float)

    baselines = {
        "median_effect": np.argsort(-median_effects),
        "trimmed_mean_effect": np.argsort(-trimmed_effects),
        "experiment_count": np.argsort(-n_exp),
    }

    itp_data = load_itp()
    itp_pos = set(d["compound"] for d in itp_data["compound_details"] if d["itp_direction"] == "positive")
    itp_neg = set(d["compound"] for d in itp_data["compound_details"] if d["itp_direction"] == "negative")

    top_q = len(rankings) // 4
    results = {}
    for name, rank_arr in baselines.items():
        tau, tau_p = kendalltau(rank_R, rank_arr)
        top_q_compounds = set(comp["compounds"][i] for i in rank_arr[:top_q])
        itp_pos_in_top_q = len(itp_pos & top_q_compounds)
        itp_neg_in_top_q = len(itp_neg & top_q_compounds)
        results[name] = {
            "kendall_tau_vs_robustness": round(float(tau), 4),
            "itp_positive_in_top_quartile": itp_pos_in_top_q,
            "itp_negative_in_top_quartile": itp_neg_in_top_q,
        }

    robustness_top_q = set(comp["compounds"][i] for i in rank_R[:top_q])
    results["robustness_score"] = {
        "kendall_tau_vs_robustness": 1.0,
        "itp_positive_in_top_quartile": len(itp_pos & robustness_top_q),
        "itp_negative_in_top_quartile": len(itp_neg & robustness_top_q),
    }

    return results


# ============================================================
# 9. Component score correlation matrix
# ============================================================
def component_correlations(profiles: list[dict]) -> dict:
    comp = get_component_scores(profiles)
    matrix = np.column_stack([
        comp["breadth"], comp["sign_consistency"],
        comp["loso"], comp["loto"],
        comp["aggregation_stability"], comp["magnitude"],
    ])
    n = matrix.shape[1]
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r, _ = spearmanr(matrix[:, i], matrix[:, j])
                corr[i, j] = r if np.isfinite(r) else 0.0

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n))
    ax.set_xticklabels(WEIGHT_NAMES, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(WEIGHT_NAMES, fontsize=9)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(corr[i, j]) > 0.5 else "black")
    fig.colorbar(im, label="Spearman ρ")
    ax.set_title("Component score correlation matrix")
    fig.tight_layout()
    fig.savefig(FIGURES / "component_correlation_matrix.png", dpi=150)
    plt.close(fig)

    corr_dict = {WEIGHT_NAMES[i]: {WEIGHT_NAMES[j]: round(corr[i, j], 4) for j in range(n)} for i in range(n)}
    return corr_dict


# ============================================================
# 14. Population statistics
# ============================================================
def population_statistics(profiles: list[dict]) -> dict:
    n_exp = [int(p.get("num_experiments", 1)) for p in profiles]
    n_sp = [int(p.get("num_species", 1)) for p in profiles]
    n_taxa = [int(p.get("num_taxa", 1)) for p in profiles]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, data, label in zip(axes, [n_exp, n_sp, n_taxa],
                                ["Experiments per compound", "Species per compound", "Taxa per compound"]):
        ax.hist(data, bins=min(30, max(data)), edgecolor="black", alpha=0.7)
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
        ax.axvline(np.median(data), color="red", linestyle="--", label=f"median={np.median(data):.0f}")
        ax.legend(fontsize=8)
    fig.suptitle("DrugAge compound population statistics")
    fig.tight_layout()
    fig.savefig(FIGURES / "population_statistics.png", dpi=150)
    plt.close(fig)

    multi_species = sum(1 for n in n_sp if n >= 2)
    multi_taxa = sum(1 for n in n_taxa if n >= 2)

    return {
        "total_compounds": len(profiles),
        "multi_species_count": multi_species,
        "multi_species_fraction": round(multi_species / len(profiles), 4),
        "multi_taxa_count": multi_taxa,
        "multi_taxa_fraction": round(multi_taxa / len(profiles), 4),
        "median_experiments": float(np.median(n_exp)),
        "mean_experiments": round(float(np.mean(n_exp)), 1),
        "max_experiments": int(np.max(n_exp)),
    }


# ============================================================
# Pipeline schematic (simple text-based)
# ============================================================
def pipeline_schematic() -> None:
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis("off")
    steps = [
        "DrugAge CSV\n(3,423 rows)",
        "Normalize\nspecies/compounds",
        "Evidence\nprofiling",
        "Robustness\nscoring (R)",
        "Tier\nassignment",
        "Perturbation\npanel",
        "Permutation\nnull (1000x)",
        "Ranked\noutput",
    ]
    n = len(steps)
    for i, step in enumerate(steps):
        x = i / (n - 1)
        color = "#3498db" if i < n - 1 else "#2ecc71"
        ax.text(x, 0.5, step, ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", facecolor=color, alpha=0.3))
        if i < n - 1:
            ax.annotate("", xy=((i + 0.7) / (n - 1), 0.5), xytext=((i + 0.3) / (n - 1), 0.5),
                         arrowprops=dict(arrowstyle="->", color="gray"))
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(FIGURES / "pipeline_schematic.png", dpi=150)
    plt.close(fig)


def main() -> int:
    rankings = load_rankings()
    profiles = load_profiles()

    results = {}

    print("1. Scatter plot: n_experiments vs R...")
    results["scatter"] = scatter_n_exp_vs_r(rankings, profiles)
    print(f"   Spearman ρ = {results['scatter']['spearman_rho']}")

    print("2. Residualized robustness score...")
    results["residualized"] = residualized_score(rankings, profiles)
    print(f"   R² of log(n_exp) → R: {results['residualized']['regression_r_squared']}")
    print(f"   Kendall τ (original vs adjusted): {results['residualized']['kendall_tau_original_vs_adjusted']}")

    print("3. Weight sensitivity sweep...")
    results["weight_sensitivity"] = weight_sensitivity(profiles, n_samples=200)
    print(f"   Median τ = {results['weight_sensitivity']['median_tau']}")

    print("4. Fisher's exact test on ITP 2x2...")
    results["fisher_itp"] = fisher_itp(rankings)
    print(f"   Odds ratio = {results['fisher_itp']['odds_ratio']}, p = {results['fisher_itp']['fisher_p_one_sided']}")

    print("5. ITP experiment counts...")
    results["itp_counts"] = itp_experiment_counts(profiles)
    print(f"   ITP-pos mean: {results['itp_counts']['itp_positive_mean_experiments']}, ITP-neg mean: {results['itp_counts']['itp_negative_mean_experiments']}")

    print("6. Naive baseline comparison...")
    results["baselines"] = naive_baselines(rankings, profiles)
    for name, vals in results["baselines"].items():
        print(f"   {name}: τ={vals['kendall_tau_vs_robustness']}, ITP-pos top-Q={vals['itp_positive_in_top_quartile']}, ITP-neg top-Q={vals['itp_negative_in_top_quartile']}")

    print("9. Component correlation matrix...")
    results["correlations"] = component_correlations(profiles)

    print("14. Population statistics...")
    results["population"] = population_statistics(profiles)
    print(f"   Multi-species: {results['population']['multi_species_count']}/{results['population']['total_compounds']} ({results['population']['multi_species_fraction']:.1%})")

    print("Pipeline schematic...")
    pipeline_schematic()

    output_path = CANONICAL / "ceiling_analyses.json"
    output_path.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nWrote {output_path}")
    print(f"Figures in {FIGURES}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
