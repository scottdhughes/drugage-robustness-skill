"""Project constants."""

from __future__ import annotations

from pathlib import Path

DEFAULT_CONFIG = "config/canonical_drugage.yaml"
DEFAULT_MANIFEST_SCHEMA = Path("schema/manifest.schema.json")
TIER_ORDER = ("robust", "promising", "thin_evidence", "conflicted")
CLAIM_STABILITY_PERTURBATIONS = (
    "leave_one_species_out_positive",
    "leave_one_taxon_out_positive",
    "median_vs_trimmed_mean_positive",
    "exclude_single_pmid_compounds",
    "exclude_mixed_sign_compounds",
)
RANKING_LEVEL_NULL_METRICS = (
    "top10_mean_robustness_score",
    "robust_compound_count",
    "top10_mean_breadth_score",
)
REQUIRED_CANONICAL_OUTPUTS = (
    "manifest",
    "normalization_audit",
    "robustness_rankings",
    "compound_evidence_profiles",
    "claim_stability_certificate",
    "claim_stability_heatmap",
    "empirical_null_certificate",
    "compound_null_significance",
    "null_separation_plot",
)
OPTIONAL_ANAGE_OUTPUTS = (
    "anage_species_join",
    "compound_taxonomic_context",
    "taxonomic_breadth_heatmap",
)
