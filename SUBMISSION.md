# Submission Checklist

This repository is prepared for two related but distinct actions:

1. Claw4S conference submission
2. clawRxiv publication

They are not the same thing.

## Claw4S Conference Submission

Visible live conference requirements:

- `SKILL.md`
- a `Research Note` in LaTeX format
- research note length `1-4 pages`
- Claw must appear as first author or corresponding-author co-author

Current repo status:

- [x] [SKILL.md](/Users/scott/sci_op/drugage_robustness_skill/SKILL.md) exists
- [x] [paper/main.tex](/Users/scott/sci_op/drugage_robustness_skill/paper/main.tex) exists
- [x] [paper/main.pdf](/Users/scott/sci_op/drugage_robustness_skill/paper/main.pdf) exists
- [x] `paper/main.pdf` is `3` pages
- [x] `Claw` is marked as corresponding author in the LaTeX note
- [x] canonical scored path passes with frozen environment

Canonical scored commands:

```bash
uv sync --frozen --extra dev
uv run --frozen --no-sync drugage-skill run --config config/canonical_drugage.yaml --out outputs/canonical
uv run --frozen --no-sync drugage-skill verify --run-dir outputs/canonical
```

Conference bundle to carry into submission:

- `SKILL.md`
- `paper/main.tex`
- `paper/main.pdf`
- `pyproject.toml`
- `uv.lock`
- `config/canonical_drugage.yaml`
- `config/species_normalization.yaml`
- `config/compound_synonyms.yaml`
- `data/drugage_build5_2024-11-29.csv`
- `outputs/canonical/manifest.json`
- `outputs/canonical/verification.json`
- `outputs/canonical/claim_stability_certificate.json`
- `outputs/canonical/empirical_null_certificate.json`
- figure artifacts in `outputs/canonical/`

## Legacy clawRxiv Publication

Prepared files:

- Markdown paper body: [paper/clawrxiv.md](/Users/scott/sci_op/drugage_robustness_skill/paper/clawrxiv.md)
- Generated payload: [submission/clawrxiv_payload.json](/Users/scott/sci_op/drugage_robustness_skill/submission/clawrxiv_payload.json)
- Payload builder: [scripts/build_clawrxiv_payload.py](/Users/scott/sci_op/drugage_robustness_skill/scripts/build_clawrxiv_payload.py)
- Register helper: [scripts/register_clawrxiv_agent.sh](/Users/scott/sci_op/drugage_robustness_skill/scripts/register_clawrxiv_agent.sh)
- Submit helper: [scripts/submit_clawrxiv.sh](/Users/scott/sci_op/drugage_robustness_skill/scripts/submit_clawrxiv.sh)

Legacy payload fields:

- required: `title`, `abstract`, `content`
- optional: `tags`, `human_names`, `skill_md`

Human action required for legacy path:

- choose a unique `claw_name`
- register and save the returned `oc_...` API key, because it is shown only once

Automatable once the key exists:

```bash
python3 scripts/build_clawrxiv_payload.py
export CLAWRXIV_API_KEY='oc_...'
scripts/submit_clawrxiv.sh
```

Live publication:

- agent name: `Claimsmith`
- paper id: `257`
- public URL: [https://www.clawrxiv.io/posts/257](https://www.clawrxiv.io/posts/257)
- human names on record: `Karen Nguyen`, `Scott Hughes`

## New clawrxiv.org API Publication

Prepared files:

- Markdown paper body: [paper/clawrxiv_org.md](/Users/scott/sci_op/drugage_robustness_skill/paper/clawrxiv_org.md)
- Generated payload: [submission/clawrxiv_org_payload.json](/Users/scott/sci_op/drugage_robustness_skill/submission/clawrxiv_org_payload.json)
- Payload builder: [scripts/build_clawrxiv_org_payload.py](/Users/scott/sci_op/drugage_robustness_skill/scripts/build_clawrxiv_org_payload.py)
- Register helper: [scripts/register_clawrxiv_org_agent.sh](/Users/scott/sci_op/drugage_robustness_skill/scripts/register_clawrxiv_org_agent.sh)
- Submit helper: [scripts/submit_clawrxiv_org.sh](/Users/scott/sci_op/drugage_robustness_skill/scripts/submit_clawrxiv_org.sh)

New API payload requirements:

- required: `title`, `abstract`, `content`, `categories`
- categories required: `1-5`
- content must include a top-level `# Title` heading that matches the payload title
- content must include a `## Abstract` section
- no images
- markdown only

Prepared categories:

- `sci.bio`
- `agents.tools`
- `ml.benchmarks`

Human action required for new API path:

1. choose agent `name` and unique `handle`
2. register the agent
3. visit the returned `verification_url`
4. complete GitHub verification to activate the API key

Automatable once the key exists:

```bash
python3 scripts/build_clawrxiv_org_payload.py
export CLAWRXIV_ORG_API_KEY='clrx_...'
scripts/submit_clawrxiv_org.sh
```

Live publication:

- agent name: `Claimsmith`
- agent handle: `claimsmith`
- paper id: `clawrxiv:2026.00014`
- public URL: [https://clawrxiv.org/papers/2026.00014](https://clawrxiv.org/papers/2026.00014)
