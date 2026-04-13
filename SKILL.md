---
name: drugage-evidence-robustness-index
description: Execute the DrugAge evidence-robustness index pipeline, including canonical scoring, ITP holdout validation, temporal Build 4→5 holdout, sensitivity analyses, and paper artifact generation.
allowed-tools: Bash(uv *, python *, ls *, test *, shasum *)
requires_python: "3.12.12"
package_manager: uv (0.9.22)
repo_root: .
canonical_output_dir: outputs/canonical
---

# Evidence-Robustness Index for Longevity Interventions in DrugAge

This skill reproduces every number in the paper from a cold start. No network access is required after cloning — all data is bundled.

## Runtime Expectations

- Platform: CPU-only
- Python: 3.12.12
- Package manager: `uv` 0.9.22
- Canonical input: `data/drugage_build5_2024-11-29.csv` (DrugAge Build 5)
- Temporal input: `data/drugage_build4_2021-11-20.csv` (DrugAge Build 4)
- Offline execution: no network access required
- Wall-clock time: ~7 minutes for canonical pipeline (1,000 permutations)

## Step 1: Confirm Bundled Data

```bash
shasum -a 256 data/drugage_build5_2024-11-29.csv
shasum -a 256 data/drugage_build4_2021-11-20.csv
```

Expected SHA256:

- Build 5: `7ed9771440fa4e1e30be0d3c8e92d919254b572ab40c81e2440ba78c885401d4`
- Build 4: `59a372489690512935df74e893f8ea7876d39f14827180e801c1ac3463bb173f`

## Step 2: Install the Locked Environment

```bash
uv sync --frozen
```

## Step 3: Run the Canonical Pipeline

```bash
uv run --frozen --no-sync drugage-skill run \
  --config config/canonical_drugage.yaml \
  --out outputs/canonical
```

Success: `outputs/canonical/verification.json` exists with status `passed`.

## Step 4: Verify

```bash
uv run --frozen --no-sync drugage-skill verify --run-dir outputs/canonical
```

## Step 5: Run Ceiling Analyses

These scripts generate every sensitivity, validation, and diagnostic artifact referenced in the paper.

```bash
# ITP overlap analysis
uv run --frozen --no-sync python scripts/itp_overlap_analysis.py

# ITP holdout (removes all ITP rows, reruns scoring)
uv run --frozen --no-sync drugage-skill run \
  --config config/itp_holdout_drugage.yaml \
  --out outputs/itp_holdout/run

# Bayesian model + GLMM
uv run --frozen --no-sync python scripts/bayesian_robustness.py

# Dose-response analysis
uv run --frozen --no-sync python scripts/dose_response_analysis.py

# Weight/cap sensitivity, baselines, component correlations, population stats
uv run --frozen --no-sync python scripts/ceiling_analyses.py

# Temporal holdout, leave-mouse-out, negative injection, PMID bootstrap, species weighting
uv run --frozen --no-sync python scripts/phase1_analyses.py
```

## Step 6: Confirm Required Artifacts

Canonical pipeline outputs:

- `outputs/canonical/manifest.json`
- `outputs/canonical/robustness_rankings.csv`
- `outputs/canonical/compound_evidence_profiles.csv`
- `outputs/canonical/empirical_null_certificate.json`
- `outputs/canonical/claim_stability_certificate.json`
- `outputs/canonical/verification.json`

Ceiling analysis outputs:

- `outputs/canonical/itp_overlap_analysis.json`
- `outputs/canonical/itp_holdout_auroc.json`
- `outputs/canonical/itp_holdout_validation.json`
- `outputs/canonical/clean_temporal_holdout.json`
- `outputs/canonical/bayesian_robustness.json`
- `outputs/canonical/ceiling_analyses.json`
- `outputs/canonical/phase1_analyses.json`
- `outputs/canonical/threshold_sensitivity.json`
- `outputs/canonical/cap_sensitivity.json`
- `outputs/canonical/final_analyses.json`
- `outputs/canonical/temporal_confound_analysis.json`
- `outputs/canonical/dose_response_analysis.json`
- `outputs/canonical/funnel_analysis.json`
- `outputs/canonical/compound_class_analysis.json`

## Success Criteria

The pipeline is successful only if:

- the canonical run finishes and verifies
- all ceiling analysis JSONs are written
- every number in the paper can be traced to an artifact in `outputs/canonical/`
