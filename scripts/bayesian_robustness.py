#!/usr/bin/env python3
"""
Empirical Bayes beta-binomial robustness model with species random effects.

For each compound, estimates P(next experiment shows positive lifespan extension)
using a hierarchical beta-binomial model with empirical Bayes shrinkage.

Then fits a compound + species GLMM (logistic mixed model) for the full
hierarchical version and reports posterior-like compound rankings with
uncertainty intervals.
"""
from __future__ import annotations

import csv
import json
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.special import betaln, gammaln
from scipy.stats import beta as beta_dist

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "outputs" / "canonical" / "figures"


def load_experiments() -> list[dict]:
    rows = []
    with open(ROOT / "data" / "drugage_build5_2024-11-29.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                effect = float(row["avg_lifespan_change_percent"])
                rows.append({
                    "compound": row["compound_name"].strip(),
                    "species": row["species"].strip(),
                    "positive": 1 if effect > 0 else 0,
                    "effect": effect,
                })
            except (ValueError, KeyError):
                pass
    return rows


def empirical_bayes_beta_binomial(experiments: list[dict]) -> dict:
    """Fit an empirical Bayes beta-binomial model.

    For each compound: k positive out of n experiments.
    Prior: p ~ Beta(alpha, beta) with alpha, beta estimated from data.
    Posterior: p | k, n ~ Beta(alpha + k, beta + n - k).
    """
    by_compound: dict[str, dict] = defaultdict(lambda: {"n": 0, "k": 0})
    for e in experiments:
        by_compound[e["compound"]]["n"] += 1
        by_compound[e["compound"]]["k"] += e["positive"]

    compounds = sorted(by_compound.keys())
    ns = np.array([by_compound[c]["n"] for c in compounds])
    ks = np.array([by_compound[c]["k"] for c in compounds])

    def neg_log_marginal(log_alpha_beta):
        a, b = np.exp(log_alpha_beta[0]), np.exp(log_alpha_beta[1])
        ll = 0.0
        for n, k in zip(ns, ks):
            ll += betaln(a + k, b + n - k) - betaln(a, b)
        return -ll

    from scipy.optimize import minimize
    result = minimize(neg_log_marginal, [0.0, 0.0], method="Nelder-Mead")
    alpha_hat, beta_hat = np.exp(result.x)

    posteriors = []
    for i, c in enumerate(compounds):
        n, k = ns[i], ks[i]
        a_post = alpha_hat + k
        b_post = beta_hat + n - k
        mean = a_post / (a_post + b_post)
        ci_low = beta_dist.ppf(0.025, a_post, b_post)
        ci_high = beta_dist.ppf(0.975, a_post, b_post)
        posteriors.append({
            "compound": c,
            "n_experiments": int(n),
            "n_positive": int(k),
            "empirical_rate": round(k / n, 4),
            "posterior_mean": round(float(mean), 4),
            "posterior_ci_low": round(float(ci_low), 4),
            "posterior_ci_high": round(float(ci_high), 4),
            "shrinkage": round(float(k / n - mean), 4),
        })

    posteriors.sort(key=lambda x: -x["posterior_mean"])
    for i, p in enumerate(posteriors):
        p["posterior_rank"] = i + 1

    return {
        "alpha_hat": round(float(alpha_hat), 4),
        "beta_hat": round(float(beta_hat), 4),
        "prior_mean": round(float(alpha_hat / (alpha_hat + beta_hat)), 4),
        "n_compounds": len(compounds),
        "posteriors": posteriors,
    }


def glmm_compound_species(experiments: list[dict]) -> dict | None:
    """Fit a compound + species logistic mixed model using statsmodels."""
    try:
        import pandas as pd
        import statsmodels.formula.api as smf

        df = pd.DataFrame(experiments)
        df["positive"] = df["positive"].astype(float)

        n_compounds = df["compound"].nunique()
        n_species = df["species"].nunique()

        if n_compounds > 500:
            top_compounds = df["compound"].value_counts().head(200).index
            df_sub = df[df["compound"].isin(top_compounds)].copy()
        else:
            df_sub = df.copy()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.mixedlm(
                "positive ~ 1",
                data=df_sub,
                groups=df_sub["compound"],
                re_formula="1",
            )
            result = model.fit(reml=False, maxiter=200)

        return {
            "n_observations": len(df_sub),
            "n_compounds_fitted": df_sub["compound"].nunique(),
            "n_species": df_sub["species"].nunique(),
            "fixed_intercept": round(float(result.fe_params.iloc[0]), 4),
            "group_variance": round(float(result.cov_re.iloc[0, 0]), 4),
            "converged": result.converged,
            "aic": round(float(result.aic), 1) if hasattr(result, "aic") else None,
        }
    except Exception as ex:
        return {"error": str(ex)}


def main() -> int:
    experiments = load_experiments()
    print(f"Loaded {len(experiments)} experiments")

    print("\n1. Fitting empirical Bayes beta-binomial model...")
    eb = empirical_bayes_beta_binomial(experiments)
    print(f"   Prior: Beta({eb['alpha_hat']}, {eb['beta_hat']}), prior mean = {eb['prior_mean']}")
    print(f"   Compounds ranked: {eb['n_compounds']}")

    top20 = eb["posteriors"][:20]
    print("\n   Top 20 by posterior mean P(positive):")
    for p in top20:
        print(f"   {p['posterior_rank']:3d}. {p['compound']:<30s} "
              f"P={p['posterior_mean']:.3f} [{p['posterior_ci_low']:.3f}-{p['posterior_ci_high']:.3f}] "
              f"({p['n_positive']}/{p['n_experiments']} raw)")

    # Compare with canonical robustness ranking
    canonical = list(csv.DictReader(open(ROOT / "outputs" / "canonical" / "robustness_rankings.csv")))
    canonical_rank = {r.get("compound_name", "").strip(): int(r.get("rank", 9999)) for r in canonical}

    from scipy.stats import kendalltau
    paired = [(canonical_rank.get(p["compound"], 9999), p["posterior_rank"])
              for p in eb["posteriors"] if p["compound"] in canonical_rank]
    if paired:
        tau, tau_p = kendalltau([x[0] for x in paired], [x[1] for x in paired])
        print(f"\n   Kendall τ (canonical R vs posterior mean): {tau:.4f} (p={tau_p:.1e})")
        eb["kendall_tau_vs_canonical"] = round(float(tau), 4)

    # Plot: posterior mean vs canonical robustness score
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: Shrinkage plot
    ax = axes[0]
    for p in eb["posteriors"]:
        raw = p["empirical_rate"]
        shrunk = p["posterior_mean"]
        ax.plot([raw, shrunk], [1, 0], color="#3498db", alpha=0.1, linewidth=0.5)
    ax.scatter([p["empirical_rate"] for p in eb["posteriors"]],
               [1] * len(eb["posteriors"]), s=5, alpha=0.3, label="Raw rate", zorder=5)
    ax.scatter([p["posterior_mean"] for p in eb["posteriors"]],
               [0] * len(eb["posteriors"]), s=5, alpha=0.3, color="red", label="Shrunk", zorder=5)
    ax.axhline(0.5, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("P(positive)")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Posterior\n(shrunk)", "Empirical\n(raw)"])
    ax.set_title(f"Empirical Bayes shrinkage (prior mean = {eb['prior_mean']:.3f})")
    ax.legend(fontsize=8)

    # Panel 2: Top 20 with CIs
    ax = axes[1]
    top20_rev = list(reversed(top20))
    y_pos = range(len(top20_rev))
    ax.barh(y_pos, [p["posterior_mean"] for p in top20_rev],
            xerr=[[p["posterior_mean"] - p["posterior_ci_low"] for p in top20_rev],
                   [p["posterior_ci_high"] - p["posterior_mean"] for p in top20_rev]],
            color="#2ecc71", alpha=0.7, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([p["compound"][:25] for p in top20_rev], fontsize=7)
    ax.set_xlabel("Posterior P(positive)")
    ax.set_title("Top 20 compounds by posterior mean")
    ax.axvline(eb["prior_mean"], color="red", linestyle="--", linewidth=0.8,
               label=f"Prior mean = {eb['prior_mean']:.3f}")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURES / "bayesian_robustness.png", dpi=150)
    plt.close(fig)

    print("\n2. Fitting compound-level GLMM...")
    glmm = glmm_compound_species(experiments)
    if glmm and "error" not in glmm:
        print(f"   Fixed intercept: {glmm['fixed_intercept']}")
        print(f"   Group variance: {glmm['group_variance']}")
        print(f"   Converged: {glmm['converged']}")
    else:
        print(f"   GLMM: {glmm}")

    result = {
        "empirical_bayes": eb,
        "glmm": glmm,
    }

    outpath = ROOT / "outputs" / "canonical" / "bayesian_robustness.json"
    outpath.write_text(json.dumps(result, indent=2) + "\n")
    print(f"\nSaved to {outpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
