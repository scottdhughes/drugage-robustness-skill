"""Canonical DrugAge pipeline."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import SkillConfig
from .constants import CLAIM_STABILITY_PERTURBATIONS, RANKING_LEVEL_NULL_METRICS
from .plots import plot_binary_heatmap, plot_null_separation
from .utils import (
    collapse_whitespace,
    ensure_dir,
    normalize_key,
    now_timestamp,
    read_yaml,
    runtime_environment,
    sha256_file,
    trimmed_mean,
    write_csv,
    write_json,
)


def _significance_summary(values: pd.Series) -> str:
    counts = Counter(value.strip() or "blank" for value in values.fillna(""))
    return ";".join(f"{key}:{counts[key]}" for key in sorted(counts))


def _compute_loxo_stability(frame: pd.DataFrame, field: str, proportion: float) -> float:
    unique_values = sorted(frame[field].dropna().unique())
    if len(unique_values) < 2:
        return 0.0
    positive = 0
    for value in unique_values:
        remaining = frame.loc[frame[field] != value, "avg_effect"].to_numpy(dtype=float)
        if remaining.size == 0:
            continue
        if trimmed_mean(remaining, proportion) > 0:
            positive += 1
    return positive / len(unique_values)


def _breadth_score(num_experiments: int, num_species: int, num_taxa: int) -> float:
    return (min(num_experiments, 6) / 6 + min(num_species, 4) / 4 + min(num_taxa, 3) / 3) / 3


def _magnitude_score(trimmed_mean_effect: float) -> float:
    return min(max(trimmed_mean_effect, 0.0), 50.0) / 50.0


def _thin_evidence_reason(row: dict[str, Any]) -> str:
    reasons: list[str] = []
    if int(row["num_species"]) < 2:
        reasons.append("species_lt_2")
    if int(row["num_experiments"]) < 3:
        reasons.append("experiments_lt_3")
    if int(row["num_pmids"]) < 2:
        reasons.append("pmids_lt_2")
    if int(row["num_taxa"]) < 2:
        reasons.append("taxa_lt_2")
    return ";".join(reasons)


def _assign_evidence_tier(row: dict[str, Any]) -> tuple[str, bool, str]:
    robust = (
        int(row["num_experiments"]) >= 3
        and int(row["num_species"]) >= 2
        and int(row["num_pmids"]) >= 2
        and float(row["median_effect"]) > 0
        and float(row["trimmed_mean_effect"]) > 0
        and float(row["sign_consistency"]) >= 0.8
        and float(row["leave_one_species_out_stability"]) == 1.0
        and float(row["aggregation_stability"]) == 1.0
    )
    if robust:
        return "robust", True, ""
    promising = (
        float(row["median_effect"]) > 0
        and float(row["trimmed_mean_effect"]) > 0
        and float(row["sign_consistency"]) >= 0.67
        and (int(row["num_species"]) >= 2 or int(row["num_experiments"]) >= 3)
    )
    if promising:
        return "promising", False, ""
    thin = (
        float(row["median_effect"]) > 0
        and float(row["trimmed_mean_effect"]) > 0
        and float(row["sign_consistency"]) >= 0.67
    )
    if thin:
        return "thin_evidence", False, _thin_evidence_reason(row)
    return "conflicted", False, ""


def _aggregate_compounds(frame: pd.DataFrame, config: SkillConfig) -> pd.DataFrame:
    proportion = float(config.pipeline["trimmed_mean_fraction"])
    rows: list[dict[str, Any]] = []
    for compound_name, group in frame.groupby("compound_name", sort=False):
        effects = group["avg_effect"].to_numpy(dtype=float)
        median_effect = float(np.median(effects))
        trimmed_mean_effect = float(trimmed_mean(effects, proportion))
        num_experiments = int(group.shape[0])
        num_species = int(group["canonical_species"].nunique())
        num_taxa = int(group["taxon_label"].nunique())
        pmids = {value.strip() for value in group["pubmed_id"].fillna("") if value and value.strip()}
        num_pmids = len(pmids)
        sign_consistency = float((effects > 0).mean()) if effects.size else 0.0
        leave_one_species_out_stability = _compute_loxo_stability(group, "canonical_species", proportion)
        leave_one_taxon_out_stability = _compute_loxo_stability(group, "taxon_label", proportion)
        aggregation_stability = 1.0 if median_effect > 0 and trimmed_mean_effect > 0 else 0.0
        breadth_score = _breadth_score(num_experiments, num_species, num_taxa)
        magnitude_score = _magnitude_score(trimmed_mean_effect)
        robustness_score = (
            0.35 * breadth_score
            + 0.20 * sign_consistency
            + 0.15 * leave_one_species_out_stability
            + 0.10 * leave_one_taxon_out_stability
            + 0.15 * aggregation_stability
            + 0.05 * magnitude_score
        )
        row = {
            "compound_name": compound_name,
            "num_experiments": num_experiments,
            "num_species": num_species,
            "num_taxa": num_taxa,
            "num_pmids": num_pmids,
            "median_effect": median_effect,
            "trimmed_mean_effect": trimmed_mean_effect,
            "sign_consistency": sign_consistency,
            "leave_one_species_out_stability": leave_one_species_out_stability,
            "leave_one_taxon_out_stability": leave_one_taxon_out_stability,
            "aggregation_stability": aggregation_stability,
            "breadth_score": breadth_score,
            "robustness_score": robustness_score,
            "avg_significance_summary": _significance_summary(group["avg_lifespan_significance"]),
            "max_significance_summary": _significance_summary(group["max_lifespan_significance"]),
            "taxa_present": ";".join(sorted(group["taxon_label"].unique())),
            "species_present": ";".join(sorted(group["canonical_species"].unique())),
            "pmids": ";".join(sorted(pmids)),
        }
        tier, robust_flag, thin_reason = _assign_evidence_tier(row)
        row["evidence_tier"] = tier
        row["robust_flag"] = robust_flag
        row["thin_evidence_reason"] = thin_reason
        rows.append(row)
    ranked = pd.DataFrame(rows)
    tier_priority = dict(config.pipeline["tier_priority"])
    ranked["_tier_priority"] = ranked["evidence_tier"].map(tier_priority)
    ranked = ranked.sort_values(
        by=[
            "_tier_priority",
            "robustness_score",
            "num_species",
            "num_taxa",
            "num_pmids",
            "trimmed_mean_effect",
            "compound_name",
        ],
        ascending=[True, False, False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    ranked["rank"] = np.arange(1, ranked.shape[0] + 1)
    return ranked.drop(columns="_tier_priority")


def _claim_stability_certificate(ranked: pd.DataFrame, config: SkillConfig) -> tuple[dict[str, Any], pd.DataFrame]:
    top_n = int(config.pipeline["top_n_compounds"])
    evaluated = ranked.head(top_n).copy()
    rows: list[dict[str, Any]] = []
    for _, row in evaluated.iterrows():
        result = {
            "compound_name": row["compound_name"],
            "evidence_tier": row["evidence_tier"],
            "leave_one_species_out_positive": int(float(row["leave_one_species_out_stability"]) == 1.0),
            "leave_one_taxon_out_positive": int(float(row["leave_one_taxon_out_stability"]) == 1.0),
            "median_vs_trimmed_mean_positive": int(
                float(row["median_effect"]) > 0 and float(row["trimmed_mean_effect"]) > 0
            ),
            "exclude_single_pmid_compounds": int(
                int(row["num_pmids"]) >= 2 and float(row["median_effect"]) > 0 and float(row["trimmed_mean_effect"]) > 0
            ),
            "exclude_mixed_sign_compounds": int(
                float(row["sign_consistency"]) == 1.0
                and float(row["median_effect"]) > 0
                and float(row["trimmed_mean_effect"]) > 0
            ),
        }
        result["all_checks_passed"] = int(all(result[name] == 1 for name in CLAIM_STABILITY_PERTURBATIONS))
        rows.append(result)
    certificate_frame = pd.DataFrame(rows)
    summary = {
        "evaluated_compounds": int(certificate_frame.shape[0]),
        "all_perturbations": len(CLAIM_STABILITY_PERTURBATIONS),
        "passed_counts": {
            name: int(certificate_frame[name].sum()) for name in CLAIM_STABILITY_PERTURBATIONS
        },
        "all_checks_passed_count": int(certificate_frame["all_checks_passed"].sum()),
        "compounds": certificate_frame.to_dict(orient="records"),
    }
    return summary, certificate_frame


def _ranking_level_null_stats(observed: float, null_values: list[float]) -> dict[str, float | None]:
    null_array = np.asarray(null_values, dtype=float)
    exceedances = int(np.sum(null_array >= observed))
    null_std = float(null_array.std(ddof=0))
    return {
        "observed": float(observed),
        "null_mean": float(null_array.mean()),
        "null_std": null_std,
        "null_exceedance_rate": exceedances / len(null_values),
        "empirical_p_value": (1 + exceedances) / (1 + len(null_values)),
        "null_z_score": None if null_std == 0 else float((observed - float(null_array.mean())) / null_std),
    }


def _permute_effects_within_species(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    species_groups = {species: idx.to_numpy() for species, idx in frame.groupby("canonical_species").groups.items()}
    permuted = frame["avg_effect"].to_numpy(dtype=float).copy()
    for idx in species_groups.values():
        permuted[idx] = rng.permutation(permuted[idx])
    return permuted


def _compute_empirical_null(
    frame: pd.DataFrame,
    ranked: pd.DataFrame,
    config: SkillConfig,
) -> tuple[dict[str, Any], pd.DataFrame, list[float]]:
    rng = np.random.default_rng(int(config.pipeline["null_seed"]))
    num_reruns = int(config.pipeline["null_reruns"])
    null_scores: dict[str, list[float]] = {name: [] for name in ranked["compound_name"]}
    ranking_level_values = {name: [] for name in RANKING_LEVEL_NULL_METRICS}
    for _ in range(num_reruns):
        null_frame = frame.copy()
        null_frame["avg_effect"] = _permute_effects_within_species(frame, rng)
        null_ranked = _aggregate_compounds(null_frame, config)
        score_map = dict(zip(null_ranked["compound_name"], null_ranked["robustness_score"], strict=True))
        for compound_name in null_scores:
            null_scores[compound_name].append(float(score_map[compound_name]))
        ranking_level_values["top10_mean_robustness_score"].append(
            float(null_ranked.head(10)["robustness_score"].mean())
        )
        ranking_level_values["robust_compound_count"].append(
            float((null_ranked["evidence_tier"] == "robust").sum())
        )
        ranking_level_values["top10_mean_breadth_score"].append(float(null_ranked.head(10)["breadth_score"].mean()))
    observed_scores = ranked.set_index("compound_name")
    null_records: list[dict[str, Any]] = []
    for compound_name, values in null_scores.items():
        observed = float(observed_scores.loc[compound_name, "robustness_score"])
        null_array = np.asarray(values, dtype=float)
        exceedances = int(np.sum(null_array >= observed))
        null_std = float(null_array.std(ddof=0))
        null_records.append(
            {
                "compound_name": compound_name,
                "observed_robustness_score": observed,
                "null_mean_score": float(null_array.mean()),
                "null_std_score": null_std,
                "null_exceedance_rate": exceedances / num_reruns,
                "empirical_p_value": (1 + exceedances) / (1 + num_reruns),
                "null_z_score": None if null_std == 0 else float((observed - float(null_array.mean())) / null_std),
            }
        )
    null_frame = pd.DataFrame(null_records).sort_values("compound_name").reset_index(drop=True)
    ranking_summary = {
        metric: _ranking_level_null_stats(
            float(
                ranked.head(10)["robustness_score"].mean()
                if metric == "top10_mean_robustness_score"
                else (ranked["evidence_tier"] == "robust").sum()
                if metric == "robust_compound_count"
                else ranked.head(10)["breadth_score"].mean()
            ),
            values,
        )
        for metric, values in ranking_level_values.items()
    }
    summary = {
        "null_model_type": str(config.pipeline["null_model_type"]),
        "num_null_reruns": num_reruns,
        "null_seed": int(config.pipeline["null_seed"]),
        "ranking_level_metrics": ranking_summary,
        "top_compounds": ranked.head(int(config.pipeline["top_n_compounds"]))[
            ["compound_name", "evidence_tier", "robustness_score"]
        ].to_dict(orient="records"),
    }
    return summary, null_frame, ranking_level_values["top10_mean_robustness_score"]


def _representative_compound_names(frame: pd.DataFrame) -> dict[str, str]:
    names: dict[str, str] = {}
    for raw in frame["compound_name_raw"]:
        key = normalize_key(raw)
        if key not in names:
            names[key] = collapse_whitespace(raw)
    return names


def _load_and_prepare_data(config: SkillConfig) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    dataset_path = config.root_dir / config.dataset["path"]
    if sha256_file(dataset_path) != config.dataset["sha256"]:
        raise ValueError("Canonical DrugAge snapshot hash mismatch")
    frame = pd.read_csv(dataset_path, dtype=str).fillna("")
    missing_columns = sorted(set(config.dataset["required_columns"]) - set(frame.columns))
    if missing_columns:
        raise ValueError(f"Missing required DrugAge columns: {missing_columns}")

    species_config = read_yaml(config.root_dir / "config" / "species_normalization.yaml")["species"]
    species_map = {normalize_key(key): value for key, value in species_config.items()}
    compound_overrides = read_yaml(config.root_dir / "config" / "compound_synonyms.yaml").get("synonyms", {})
    compound_map = {normalize_key(key): collapse_whitespace(value) for key, value in compound_overrides.items()}
    representative_names = _representative_compound_names(frame.rename(columns={"compound_name": "compound_name_raw"}))

    included_rows: list[dict[str, Any]] = []
    drop_reasons = Counter()
    species_overrides_used: dict[str, dict[str, str]] = {}
    compound_overrides_used: dict[str, str] = {}

    for _, record in frame.iterrows():
        raw_compound = collapse_whitespace(str(record["compound_name"]))
        raw_species = collapse_whitespace(str(record["species"]))
        effect_text = collapse_whitespace(str(record["avg_lifespan_change_percent"]))

        if not raw_compound:
            drop_reasons["missing_compound_name"] += 1
            continue
        if not raw_species:
            drop_reasons["missing_species"] += 1
            continue
        try:
            avg_effect = float(effect_text)
        except ValueError:
            drop_reasons["missing_or_non_numeric_avg_lifespan_change_percent"] += 1
            continue

        species_key = normalize_key(raw_species)
        if species_key not in species_map:
            raise ValueError(f"Species normalization missing for {raw_species}")
        species_entry = species_map[species_key]
        canonical_species = collapse_whitespace(str(species_entry["canonical_species"]))
        taxon_label = collapse_whitespace(str(species_entry["taxon_label"]))
        if raw_species != canonical_species:
            species_overrides_used[raw_species] = {
                "canonical_species": canonical_species,
                "taxon_label": taxon_label,
            }

        compound_key = normalize_key(raw_compound)
        canonical_compound = compound_map.get(compound_key, representative_names[compound_key])
        if raw_compound != canonical_compound:
            compound_overrides_used[raw_compound] = canonical_compound

        included_rows.append(
            {
                "compound_name": canonical_compound,
                "compound_name_raw": raw_compound,
                "species_raw": raw_species,
                "canonical_species": canonical_species,
                "taxon_label": taxon_label,
                "avg_effect": avg_effect,
                "pubmed_id": collapse_whitespace(str(record["pubmed_id"])),
                "avg_lifespan_significance": collapse_whitespace(str(record["avg_lifespan_significance"])),
                "max_lifespan_significance": collapse_whitespace(str(record["max_lifespan_significance"])),
                "max_lifespan_change_percent": collapse_whitespace(str(record["max_lifespan_change_percent"])),
                "dosage": collapse_whitespace(str(record["dosage"])),
                "gender": collapse_whitespace(str(record["gender"])),
                "ITP": collapse_whitespace(str(record["ITP"])),
            }
        )

    included = pd.DataFrame(included_rows)
    inclusion_summary = {
        "raw_rows": int(frame.shape[0]),
        "included_rows": int(included.shape[0]),
        "excluded_rows": int(frame.shape[0] - included.shape[0]),
        "drop_reasons": dict(sorted(drop_reasons.items())),
    }
    normalization_audit = {
        "species_overrides_used": species_overrides_used,
        "compound_overrides_used": compound_overrides_used,
        "species_taxon_labels": {
            raw: {
                "canonical_species": value["canonical_species"],
                "taxon_label": value["taxon_label"],
            }
            for raw, value in species_config.items()
        },
    }
    return included, inclusion_summary, normalization_audit


def run_pipeline(config: SkillConfig, out_dir: str | Path) -> dict[str, Any]:
    out_path = ensure_dir(Path(out_dir))
    included, inclusion_summary, normalization_audit = _load_and_prepare_data(config)

    # Optional gender filter
    gender_filter = config.pipeline.get("gender_filter")
    if gender_filter:
        before = len(included)
        included = included[included["gender"] == gender_filter].reset_index(drop=True)
        inclusion_summary["gender_filter"] = gender_filter
        inclusion_summary["rows_before_gender_filter"] = before
        inclusion_summary["rows_after_gender_filter"] = len(included)

    ranked = _aggregate_compounds(included, config)
    claim_summary, claim_frame = _claim_stability_certificate(ranked, config)
    empirical_summary, null_frame, null_top10_distribution = _compute_empirical_null(included, ranked, config)

    ranked = ranked.merge(
        null_frame[
            [
                "compound_name",
                "empirical_p_value",
                "null_exceedance_rate",
                "null_z_score",
            ]
        ],
        on="compound_name",
        how="left",
    )

    manifest = {
        "metadata": {
            "title": config.title,
            "schema_version": config.contract_version,
            "run_timestamp": now_timestamp(),
        },
        "input_provenance": {
            "path": config.dataset["path"],
            "sha256": config.dataset["sha256"],
            "source_url": config.dataset["source_url"],
            "build": config.dataset["build"],
            "release_date": str(config.dataset["release_date"]),
            "row_count": int(inclusion_summary["raw_rows"]),
            "required_columns": list(config.dataset["required_columns"]),
        },
        "environment": runtime_environment(
            package_names=list(config.runtime["package_versions"]),
            pythonhashseed=int(config.runtime["pythonhashseed"]),
        ),
        "inclusion_summary": inclusion_summary,
        "normalization": {
            "species_overrides_used": normalization_audit["species_overrides_used"],
            "compound_overrides_used": normalization_audit["compound_overrides_used"],
        },
        "ranking_contract": {
            "contract_version": config.contract_version,
            "trimmed_mean_fraction": float(config.pipeline["trimmed_mean_fraction"]),
            "tier_priority": dict(config.pipeline["tier_priority"]),
            "null_model_type": config.pipeline["null_model_type"],
        },
        "evidence_tiers": {
            tier: int((ranked["evidence_tier"] == tier).sum())
            for tier in config.pipeline["tier_priority"].keys()
        },
        "top_compounds": ranked.head(int(config.pipeline["top_n_compounds"]))[
            [
                "rank",
                "compound_name",
                "evidence_tier",
                "robust_flag",
                "robustness_score",
                "median_effect",
                "trimmed_mean_effect",
                "num_species",
                "num_taxa",
                "num_pmids",
                "empirical_p_value",
            ]
        ].to_dict(orient="records"),
        "claim_stability_summary": claim_summary,
        "empirical_null_summary": empirical_summary,
        "null_model_type": empirical_summary["null_model_type"],
        "num_null_reruns": empirical_summary["num_null_reruns"],
        "null_seed": empirical_summary["null_seed"],
        "ranking_null_summary": empirical_summary["ranking_level_metrics"],
    }

    write_csv(out_path / config.outputs["robustness_rankings"], ranked)
    write_csv(out_path / config.outputs["compound_evidence_profiles"], ranked)
    write_json(out_path / config.outputs["normalization_audit"], normalization_audit)
    write_json(out_path / config.outputs["claim_stability_certificate"], claim_summary)
    write_json(out_path / config.outputs["empirical_null_certificate"], empirical_summary)
    write_csv(out_path / config.outputs["compound_null_significance"], null_frame)
    write_json(out_path / config.outputs["manifest"], manifest)

    claim_heatmap = claim_frame.set_index("compound_name")[
        [
            "leave_one_species_out_positive",
            "leave_one_taxon_out_positive",
            "median_vs_trimmed_mean_positive",
            "exclude_single_pmid_compounds",
            "exclude_mixed_sign_compounds",
        ]
    ]
    plot_binary_heatmap(
        claim_heatmap,
        out_path / config.outputs["claim_stability_heatmap"],
        "Claim Stability Certificate",
    )
    plot_null_separation(
        null_top10_distribution,
        float(ranked.head(10)["robustness_score"].mean()),
        out_path / config.outputs["null_separation_plot"],
    )
    return manifest
