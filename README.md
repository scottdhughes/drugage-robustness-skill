# An Evidence-Robustness Index for Longevity Interventions in DrugAge

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An automated pipeline that ranks DrugAge longevity interventions by cross-species evidence consistency rather than raw lifespan effect. Validated against the NIA Interventions Testing Program (ITP) and a temporal holdout between DrugAge Build 4 (2021) and Build 5 (2024).

**Paper:** Nguyen K, Hughes S, Claw. "An Evidence-Robustness Index for Longevity Interventions in DrugAge." clawRxiv 2026.

## Key Results

- **48 robust-tier compounds** out of 1,038 scored, separated from a 1,000-permutation null at z = 4.42
- **ITP holdout AUROC = 0.847** (95% CI [0.614, 1.000]) on a clean holdout with all ITP rows removed before scoring
- **12/12 ITP-negative compounds** correctly classified as conflicted (Fisher's exact p = 0.0036)
- **Build 4 → Build 5 forward validation**: top-quartile compounds gain positive evidence at 0.886 vs 0.587 (p < 0.0001), controlled for attention confound via regression (β = 0.319, p < 10⁻⁶) and propensity matching (p = 0.000039)
- **Stable under perturbation**: weight sensitivity τ = 0.77, cap sensitivity 80.6% with τ > 0.80, PMID-clustered bootstrap τ = 0.76

## Quick Start

```bash
# Install
uv sync --frozen

# Run the canonical pipeline (1,000 permutations, ~7 min)
uv run --frozen --no-sync drugage-skill run \
  --config config/canonical_drugage.yaml \
  --out outputs/canonical

# Verify
uv run --frozen --no-sync drugage-skill verify --run-dir outputs/canonical
```

## Repository Structure

```
├── config/
│   ├── canonical_drugage.yaml      # Pipeline configuration (1,000 permutations)
│   ├── species_normalization.yaml  # Species → taxon mapping
│   └── compound_synonyms.yaml     # Compound name normalization
├── data/
│   ├── drugage_build5_2024-11-29.csv  # DrugAge Build 5 (bundled, SHA256-verified)
│   └── drugage_build4_2021-11-20.csv  # DrugAge Build 4 (for temporal validation)
├── src/drugage_skill/
│   ├── pipeline.py                 # Core scoring pipeline
│   ├── cli.py                      # Command-line interface
│   ├── config.py                   # Configuration loader
│   ├── verify.py                   # Verification logic
│   └── plots.py                    # Visualization
├── scripts/
│   ├── ceiling_analyses.py         # Weight/cap sensitivity, baselines, correlations
│   ├── phase1_analyses.py          # Temporal holdout, leave-mouse-out, neg injection
│   ├── bayesian_robustness.py      # Empirical Bayes beta-binomial model
│   ├── itp_overlap_analysis.py     # Systematic ITP cross-validation
│   ├── dose_response_analysis.py   # Within-compound dose-response
│   └── temporal_confound_analysis.py  # Attention confound controls
├── outputs/canonical/              # Pinned canonical artifacts
│   ├── robustness_rankings.csv     # Full compound rankings
│   ├── compound_evidence_profiles.csv
│   ├── empirical_null_certificate.json
│   ├── itp_overlap_analysis.json
│   ├── itp_holdout_validation.json
│   ├── itp_holdout_auroc.json
│   ├── build4_to_build5_validation.json
│   ├── temporal_confound_analysis.json
│   ├── bayesian_robustness.json
│   ├── ceiling_analyses.json
│   ├── phase1_analyses.json
│   ├── threshold_sensitivity.json
│   ├── cap_sensitivity.json
│   ├── dose_response_analysis.json
│   ├── funnel_analysis.json
│   ├── geroprotectors_overlap_analysis.json
│   ├── compound_class_analysis.json
│   └── figures/*.png
├── paper/
│   ├── main.tex                    # 4-page research note (Claw4S submission)
│   ├── main.pdf                    # Compiled PDF
│   ├── main_full_16page.tex        # Full 16-page version
│   ├── main_full_16page.pdf        # Full compiled PDF
│   └── clawrxiv.md                 # Markdown version
├── tests/                          # 5 tests
├── SKILL.md                        # Agent execution skill
├── pyproject.toml                  # Python project config (uv)
└── uv.lock                         # Locked dependencies
```

## Data Availability

- **DrugAge Build 5** (November 29, 2024): bundled at `data/drugage_build5_2024-11-29.csv`. SHA256: `7ed9771440fa4e1e30be0d3c8e92d919254b572ab40c81e2440ba78c885401d4`. Source: [HAGR DrugAge](https://genomics.senescence.info/drugs/). DrugAge is available under CC BY 3.0 with attribution; please cite Barardo et al. (2017) and the [2024 HAGR resource paper](https://academic.oup.com/nar/article/52/D1/D900/7337614).
- **DrugAge Build 4** (November 20, 2021): bundled at `data/drugage_build4_2021-11-20.csv`. Retrieved from the Wayback Machine archive of the HAGR download endpoint.
- **Geroprotectors.org** data (215 compounds, 259 experiments): scraped for cross-database validation; stored in `outputs/canonical/geroprotectors_overlap_analysis.json`.

## Robustness Score Formula

```
R = 0.35 × breadth + 0.20 × sign_consistency + 0.15 × LOSO + 0.10 × LOTO + 0.15 × aggregation_stability + 0.05 × magnitude
```

Weights prioritize breadth (0.35) and consistency (0.45 combined) over raw effect magnitude (0.05). Sensitivity analysis: median Kendall τ = 0.77 under 200 Dirichlet-sampled weight perturbations.

## Evidence Tiers

| Tier | Criteria | Count |
|------|----------|-------|
| **Robust** | ≥3 experiments, ≥2 species, ≥2 PMIDs, sign consistency ≥0.80, LOSO = 1.0, positive median + trimmed mean | 48 |
| **Promising** | Sign consistency ≥0.67, ≥2 species or ≥3 experiments | 174 |
| **Thin evidence** | Meets consistency but insufficient breadth | 435 |
| **Conflicted** | All remaining | 381 |

## Top 10 Compounds

| Rank | Compound | Score | Species | Experiments |
|------|----------|-------|---------|-------------|
| 1 | Spermidine | 1.000 | 6 | 8 |
| 2 | Apple extract | 0.950 | 3 | 6 |
| 3 | N-acetyl-L-cysteine | 0.941 | 5 | 42 |
| 4 | Minocycline | 0.929 | 3 | 28 |
| 5 | Alpha-ketoglutarate | 0.929 | 4 | 18 |
| 6 | Carnosine | 0.929 | 3 | 9 |
| 7 | Rapamycin | 0.929 | 3 | 37 |
| 8 | Mycophenolic acid | 0.913 | 3 | 7 |
| 9 | EGCG | 0.912 | 3 | 15 |
| 10 | Vitamin E | 0.907 | 4 | 31 |

## Limitations

This index measures **analytical robustness** (cross-species consistency within DrugAge), not biological geroprotective efficacy. Key quantified limitations:

- Data-volume confound: ρ = 0.44, R² = 0.22
- Publication bias: 70.1% positive experiments
- Structure-matched null: z = 1.26, p = 0.13 (not significant)
- 30% hidden-negative injection: 67% robust-tier survival
- ITP holdout tier concordance: p = 0.24 (continuous ranking holds)
- Spermidine heuristic rank 1 → Bayesian rank 17

## License

MIT

## Citation

```bibtex
@article{nguyen2026drugage,
  title={An Evidence-Robustness Index for Longevity Interventions in DrugAge},
  author={Nguyen, Karen and Hughes, Scott and Claw},
  year={2026}
}
```
