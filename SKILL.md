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
# 5a. ITP overlap analysis → itp_overlap_analysis.json
uv run --frozen --no-sync python scripts/itp_overlap_analysis.py

# 5b. ITP holdout: create filtered CSV, run scoring, compute AUROC
uv run --frozen --no-sync python -c "
import csv
rows = list(csv.DictReader(open('data/drugage_build5_2024-11-29.csv')))
non_itp = [r for r in rows if r.get('ITP','').strip() != 'Yes']
import hashlib, pathlib, io
out = pathlib.Path('outputs/itp_holdout'); out.mkdir(parents=True, exist_ok=True)
buf = io.StringIO(); w = csv.DictWriter(buf, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(non_itp)
(out / 'drugage_no_itp.csv').write_text(buf.getvalue())
"
uv run --frozen --no-sync drugage-skill run \
  --config config/itp_holdout_drugage.yaml \
  --out outputs/itp_holdout/run

# 5c. Bayesian model → bayesian_robustness.json
uv run --frozen --no-sync python scripts/bayesian_robustness.py

# 5d. Dose-response → dose_response_analysis.json
uv run --frozen --no-sync python scripts/dose_response_analysis.py

# 5e. Weight/cap sensitivity, baselines, correlations, population stats
#     → ceiling_analyses.json, figures/
uv run --frozen --no-sync python scripts/ceiling_analyses.py

# 5f. Temporal holdout, leave-mouse-out, negative injection, PMID bootstrap
#     → phase1_analyses.json
uv run --frozen --no-sync python scripts/phase1_analyses.py

# 5g. Clean temporal holdout (Build-4-only scoring)
#     → clean_temporal_holdout.json
uv run --frozen --no-sync python scripts/clean_temporal_holdout.py

# 5h. ITP holdout AUROC → itp_holdout_auroc.json
uv run --frozen --no-sync python scripts/itp_holdout_auroc.py

# 5i. Threshold sensitivity → threshold_sensitivity.json
uv run --frozen --no-sync python scripts/threshold_sensitivity.py

# 5j. Temporal confound controls → temporal_confound_analysis.json
uv run --frozen --no-sync python scripts/temporal_confound_analysis.py

# 5k. Compound class analysis → compound_class_analysis.json
uv run --frozen --no-sync python scripts/compound_class_analysis.py

# 5l. Final analyses (unseen-species, pure-compound, simple baseline)
#     → final_analyses.json
uv run --frozen --no-sync python scripts/final_analyses.py

# 5m. Funnel/positivity analysis → funnel_analysis.json
uv run --frozen --no-sync python scripts/funnel_analysis.py
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
