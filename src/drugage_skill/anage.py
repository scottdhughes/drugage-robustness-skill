"""Optional descriptive AnAge context report."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .config import SkillConfig
from .plots import plot_binary_heatmap
from .utils import collapse_whitespace, ensure_dir, normalize_key, read_json, write_csv


def run_anage_context_report(config: SkillConfig, run_dir: str | Path) -> dict[str, Any]:
    run_path = ensure_dir(Path(run_dir))
    rankings = pd.read_csv(run_path / config.outputs["robustness_rankings"])
    anage = pd.read_csv(config.root_dir / config.anage["path"], sep="\t", dtype=str).fillna("")
    anage["anage_species"] = (
        anage["Genus"].map(collapse_whitespace) + " " + anage["Species"].map(collapse_whitespace)
    ).map(collapse_whitespace)
    anage["anage_key"] = anage["anage_species"].map(normalize_key)
    anage_lookup = anage.drop_duplicates("anage_key").set_index("anage_key", drop=False)

    normalization = read_json(run_path / config.outputs["normalization_audit"])
    species_rows: list[dict[str, Any]] = []
    for raw_species, entry in normalization["species_taxon_labels"].items():
        canonical = entry["canonical_species"]
        key = normalize_key(canonical)
        matched = key in anage_lookup.index
        row = {
            "raw_species": raw_species,
            "canonical_species": canonical,
            "taxon_label": entry["taxon_label"],
            "anage_matched": matched,
        }
        if matched:
            record = anage_lookup.loc[key]
            row.update(
                {
                    "anage_common_name": record["Common name"],
                    "anage_kingdom": record["Kingdom"],
                    "anage_phylum": record["Phylum"],
                    "anage_class": record["Class"],
                    "anage_order": record["Order"],
                    "anage_family": record["Family"],
                }
            )
        species_rows.append(row)
    species_frame = pd.DataFrame(species_rows).sort_values(["anage_matched", "canonical_species"], ascending=[False, True])
    write_csv(run_path / config.outputs["anage_species_join"], species_frame)

    compound_rows: list[dict[str, Any]] = []
    species_lookup = species_frame.groupby("canonical_species", sort=False).first().reset_index()
    top = rankings.head(20)
    for _, row in top.iterrows():
        canonical_species = [item for item in str(row["species_present"]).split(";") if item]
        matched_species = species_lookup[
            species_lookup["canonical_species"].isin(canonical_species) & species_lookup["anage_matched"]
        ]
        compound_rows.append(
            {
                "compound_name": row["compound_name"],
                "rank": int(row["rank"]),
                "evidence_tier": row["evidence_tier"],
                "num_species": int(row["num_species"]),
                "matched_anage_species": int(matched_species["canonical_species"].nunique()),
                "matched_anage_classes": ";".join(sorted(set(matched_species["anage_class"]) - {""})),
                "matched_anage_orders": ";".join(sorted(set(matched_species["anage_order"]) - {""})),
                "unmatched_species": ";".join(
                    sorted(set(canonical_species) - set(matched_species["canonical_species"]))
                ),
            }
        )
    compound_frame = pd.DataFrame(compound_rows)
    write_csv(run_path / config.outputs["compound_taxonomic_context"], compound_frame)

    heatmap_rows = []
    for _, row in compound_frame.iterrows():
        row_classes = {item for item in str(row["matched_anage_classes"]).split(";") if item}
        heatmap_rows.append(
            {
                "compound_name": row["compound_name"],
                **{class_name: 1 for class_name in row_classes},
            }
        )
    heatmap = pd.DataFrame(heatmap_rows).fillna(0)
    if not heatmap.empty:
        heatmap = heatmap.set_index("compound_name")
    plot_binary_heatmap(
        heatmap,
        run_path / config.outputs["taxonomic_breadth_heatmap"],
        "AnAge Matched Class Coverage for Top Compounds",
    )
    return {
        "matched_species": int(species_frame["anage_matched"].sum()),
        "total_species": int(species_frame.shape[0]),
        "top_compounds_reported": int(compound_frame.shape[0]),
    }
