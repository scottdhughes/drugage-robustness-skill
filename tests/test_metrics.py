from copy import deepcopy

import pandas as pd

from drugage_skill.config import SkillConfig, load_config
from drugage_skill.pipeline import _aggregate_compounds, _compute_empirical_null, _permute_effects_within_species


def _synthetic_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "compound_name": "A",
                "compound_name_raw": "A",
                "species_raw": "Species 1",
                "canonical_species": "Species 1",
                "taxon_label": "taxon_1",
                "avg_effect": 12.0,
                "pubmed_id": "1",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "A",
                "compound_name_raw": "A",
                "species_raw": "Species 1",
                "canonical_species": "Species 1",
                "taxon_label": "taxon_1",
                "avg_effect": 10.0,
                "pubmed_id": "1",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "A",
                "compound_name_raw": "A",
                "species_raw": "Species 2",
                "canonical_species": "Species 2",
                "taxon_label": "taxon_2",
                "avg_effect": 8.0,
                "pubmed_id": "2",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "A",
                "compound_name_raw": "A",
                "species_raw": "Species 2",
                "canonical_species": "Species 2",
                "taxon_label": "taxon_2",
                "avg_effect": 9.0,
                "pubmed_id": "2",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "B",
                "compound_name_raw": "B",
                "species_raw": "Species 3",
                "canonical_species": "Species 3",
                "taxon_label": "taxon_3",
                "avg_effect": 5.0,
                "pubmed_id": "3",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "B",
                "compound_name_raw": "B",
                "species_raw": "Species 3",
                "canonical_species": "Species 3",
                "taxon_label": "taxon_3",
                "avg_effect": 4.5,
                "pubmed_id": "4",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "B",
                "compound_name_raw": "B",
                "species_raw": "Species 3",
                "canonical_species": "Species 3",
                "taxon_label": "taxon_3",
                "avg_effect": 4.0,
                "pubmed_id": "5",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "C",
                "compound_name_raw": "C",
                "species_raw": "Species 4",
                "canonical_species": "Species 4",
                "taxon_label": "taxon_4",
                "avg_effect": 3.0,
                "pubmed_id": "6",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "C",
                "compound_name_raw": "C",
                "species_raw": "Species 4",
                "canonical_species": "Species 4",
                "taxon_label": "taxon_4",
                "avg_effect": 2.0,
                "pubmed_id": "6",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "D",
                "compound_name_raw": "D",
                "species_raw": "Species 5",
                "canonical_species": "Species 5",
                "taxon_label": "taxon_5",
                "avg_effect": 3.0,
                "pubmed_id": "7",
                "avg_lifespan_significance": "S",
                "max_lifespan_significance": "S",
            },
            {
                "compound_name": "D",
                "compound_name_raw": "D",
                "species_raw": "Species 6",
                "canonical_species": "Species 6",
                "taxon_label": "taxon_6",
                "avg_effect": -5.0,
                "pubmed_id": "8",
                "avg_lifespan_significance": "NS",
                "max_lifespan_significance": "NS",
            },
        ]
    )


def _small_null_config() -> SkillConfig:
    config = load_config("config/canonical_drugage.yaml")
    raw = deepcopy(config.raw)
    raw["pipeline"]["null_reruns"] = 8
    raw["pipeline"]["top_n_compounds"] = 3
    raw["pipeline"]["null_seed"] = 17
    return SkillConfig(root_dir=config.root_dir, raw=raw)


def test_aggregate_compounds_assigns_expected_tiers() -> None:
    ranked = _aggregate_compounds(_synthetic_frame(), _small_null_config())
    tiers = dict(zip(ranked["compound_name"], ranked["evidence_tier"], strict=True))
    assert tiers["A"] == "robust"
    assert tiers["B"] == "promising"
    assert tiers["C"] == "thin_evidence"
    assert tiers["D"] == "conflicted"
    assert ranked["compound_name"].tolist() == ["A", "B", "C", "D"]


def test_species_permutation_preserves_within_species_distribution() -> None:
    frame = _synthetic_frame()
    permuted = _permute_effects_within_species(frame, __import__("numpy").random.default_rng(5))
    permuted_frame = frame.copy()
    permuted_frame["avg_effect"] = permuted
    for species, group in frame.groupby("canonical_species", sort=False):
        original = sorted(group["avg_effect"].tolist())
        shuffled = sorted(permuted_frame.loc[permuted_frame["canonical_species"] == species, "avg_effect"].tolist())
        assert original == shuffled


def test_empirical_null_returns_probability_columns() -> None:
    config = _small_null_config()
    ranked = _aggregate_compounds(_synthetic_frame(), config)
    summary, null_frame, top10_distribution = _compute_empirical_null(_synthetic_frame(), ranked, config)
    assert summary["null_model_type"] == "species_stratified_effect_permutation"
    assert len(top10_distribution) == 8
    assert null_frame.shape[0] == ranked.shape[0]
    assert ((null_frame["empirical_p_value"] >= 0) & (null_frame["empirical_p_value"] <= 1)).all()

