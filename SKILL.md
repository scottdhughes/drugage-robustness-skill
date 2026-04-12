---
name: drugage-robustness-null-certified
description: Execute a locked, offline DrugAge robustness ranking pipeline with evidence tiers, a claim stability certificate, and an empirical null certificate.
allowed-tools: Bash(uv *, python *, ls *, test *, shasum *)
requires_python: "3.12.x"
package_manager: uv
repo_root: .
canonical_output_dir: outputs/canonical
---

# Claim-Certified DrugAge Robustness Skill

This skill executes the canonical scored path only. It does not run the optional AnAge context report or posting helpers.

## Runtime Expectations

- Platform: CPU-only
- Python: 3.12.x
- Package manager: `uv`
- Canonical input: `data/drugage_build5_2024-11-29.csv`
- Offline execution: no network access required after the repo is cloned

## Step 1: Confirm Canonical Input

```bash
test -f data/drugage_build5_2024-11-29.csv
shasum -a 256 data/drugage_build5_2024-11-29.csv
```

Expected SHA256:

```text
7ed9771440fa4e1e30be0d3c8e92d919254b572ab40c81e2440ba78c885401d4
```

## Step 2: Install the Locked Environment

```bash
uv sync --frozen
```

Success condition:

- `uv` completes without changing the lockfile

## Step 3: Run the Canonical Pipeline

```bash
uv run --frozen --no-sync drugage-skill run --config config/canonical_drugage.yaml --out outputs/canonical
```

Success condition:

- `outputs/canonical/manifest.json` exists
- all required CSV, JSON, and PNG artifacts are present

## Step 4: Verify the Run

```bash
uv run --frozen --no-sync drugage-skill verify --run-dir outputs/canonical
```

Success condition:

- exit code is `0`
- `outputs/canonical/verification.json` exists
- verification status is `passed`

## Step 5: Confirm Required Artifacts

Required files:

- `outputs/canonical/manifest.json`
- `outputs/canonical/normalization_audit.json`
- `outputs/canonical/robustness_rankings.csv`
- `outputs/canonical/compound_evidence_profiles.csv`
- `outputs/canonical/claim_stability_certificate.json`
- `outputs/canonical/claim_stability_heatmap.png`
- `outputs/canonical/empirical_null_certificate.json`
- `outputs/canonical/compound_null_significance.csv`
- `outputs/canonical/null_separation_plot.png`
- `outputs/canonical/verification.json`

## Optional: Sex-Stratified Sensitivity Analysis

```bash
uv run --frozen --no-sync drugage-skill run --gender Male --out outputs/male_only
uv run --frozen --no-sync drugage-skill run --gender Female --out outputs/female_only
```

The `--gender` flag filters experiments to a single sex before scoring. Valid values: Male, Female, Hermaphrodite, Pooled, Unknown. This enables sensitivity analysis to assess whether robustness rankings are stable across sexes. In the canonical dataset, 6/10 top compounds are shared between male-only and female-only rankings.

## Step 6: Canonical Success Criteria

The canonical path is successful only if:

- the bundled DrugAge snapshot is used
- the run command finishes successfully
- the verify command exits `0`
- all required artifacts are present and nonempty
