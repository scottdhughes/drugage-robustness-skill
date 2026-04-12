"""Verification logic."""

from __future__ import annotations

import math
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import pandas as pd
from jsonschema import Draft202012Validator

from .config import SkillConfig
from .constants import DEFAULT_MANIFEST_SCHEMA, REQUIRED_CANONICAL_OUTPUTS
from .pipeline import run_pipeline
from .utils import read_json, write_json


def _check(name: str, passed: bool, details: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "details": details}


def _ranking_order(frame: pd.DataFrame, top_n: int) -> list[str]:
    return frame.head(top_n)["compound_name"].tolist()


def _close_enough(left: float | None, right: float | None, tolerance: float = 1e-6) -> bool:
    if left is None or right is None:
        return left is None and right is None
    return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=tolerance)


def verify_run(config: SkillConfig, run_dir: str | Path) -> dict[str, object]:
    run_path = Path(run_dir)
    manifest = read_json(run_path / config.outputs["manifest"])
    schema = read_json(config.root_dir / DEFAULT_MANIFEST_SCHEMA)
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(manifest), key=lambda item: item.path)

    checks: list[dict[str, object]] = []
    checks.append(_check("manifest_schema", not errors, "; ".join(error.message for error in errors) or "ok"))

    required_paths = [run_path / config.outputs[name] for name in REQUIRED_CANONICAL_OUTPUTS]
    artifact_ok = all(path.exists() and path.stat().st_size > 0 for path in required_paths)
    checks.append(_check("artifacts_exist", artifact_ok, "required outputs present"))

    base_rankings = pd.read_csv(run_path / config.outputs["robustness_rankings"])
    base_profiles = pd.read_csv(run_path / config.outputs["compound_evidence_profiles"])
    base_verification_top = _ranking_order(base_rankings, int(config.pipeline["verification_top_n"]))
    checks.append(_check("top20_unique", len(set(base_verification_top)) == len(base_verification_top), ",".join(base_verification_top)))

    with TemporaryDirectory(prefix="drugage-skill-verify-") as temp_dir:
        rerun_manifest = run_pipeline(config, temp_dir)
        rerun_path = Path(temp_dir)
        rerun_rankings = pd.read_csv(rerun_path / config.outputs["robustness_rankings"])
        rerun_profiles = pd.read_csv(rerun_path / config.outputs["compound_evidence_profiles"])
        rerun_claim = read_json(rerun_path / config.outputs["claim_stability_certificate"])
        rerun_null = read_json(rerun_path / config.outputs["empirical_null_certificate"])

    same_top20 = base_verification_top == _ranking_order(rerun_rankings, int(config.pipeline["verification_top_n"]))
    checks.append(_check("top20_ranking_order", same_top20, ",".join(base_verification_top)))

    tier_match = (
        base_rankings[["compound_name", "evidence_tier"]]
        .merge(rerun_rankings[["compound_name", "evidence_tier"]], on="compound_name", suffixes=("_base", "_rerun"))
        .eval("evidence_tier_base == evidence_tier_rerun")
        .all()
    )
    checks.append(_check("tier_assignments", bool(tier_match), "evidence tiers stable"))

    base_claim = read_json(run_path / config.outputs["claim_stability_certificate"])
    claim_match = base_claim["passed_counts"] == rerun_claim["passed_counts"]
    checks.append(_check("claim_stability_summary", claim_match, str(base_claim["passed_counts"])))

    base_null = read_json(run_path / config.outputs["empirical_null_certificate"])
    null_metrics_ok = True
    for metric_name, metric_payload in base_null["ranking_level_metrics"].items():
        rerun_payload = rerun_null["ranking_level_metrics"][metric_name]
        for field in ("observed", "null_mean", "null_std", "null_exceedance_rate", "empirical_p_value", "null_z_score"):
            if not _close_enough(metric_payload[field], rerun_payload[field]):
                null_metrics_ok = False
                break
        if not null_metrics_ok:
            break
    checks.append(_check("empirical_null_summary", null_metrics_ok, "ranking-level null metrics stable"))

    numeric_columns = [
        "median_effect",
        "trimmed_mean_effect",
        "sign_consistency",
        "leave_one_species_out_stability",
        "leave_one_taxon_out_stability",
        "aggregation_stability",
        "breadth_score",
        "robustness_score",
        "empirical_p_value",
        "null_exceedance_rate",
        "null_z_score",
    ]
    merged = base_profiles.merge(rerun_profiles, on="compound_name", suffixes=("_base", "_rerun"))
    profile_match = True
    for column in numeric_columns:
        for left, right in zip(merged[f"{column}_base"], merged[f"{column}_rerun"], strict=True):
            left_value = None if pd.isna(left) else float(left)
            right_value = None if pd.isna(right) else float(right)
            if not _close_enough(left_value, right_value):
                profile_match = False
                break
        if not profile_match:
            break
    checks.append(_check("profile_numeric_metrics", profile_match, "compound-level numeric metrics stable"))

    verification = {"status": "passed" if all(check["passed"] for check in checks) else "failed", "checks": checks}
    manifest["verification"] = verification
    write_json(run_path / config.outputs["manifest"], manifest)
    write_json(run_path / config.outputs["verification"], verification)
    return verification
