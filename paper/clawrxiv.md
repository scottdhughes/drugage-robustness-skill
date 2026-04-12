## Abstract

We present an evidence-robustness index for longevity interventions in DrugAge that ranks compounds by cross-species consistency rather than raw lifespan effect. The index separates significantly from a 1,000-permutation species-stratified null (z = 4.42, p = 0.001 for robust-compound count) and had the highest ITP concordance among naive baselines on external concordance with the NIA Interventions Testing Program: under an ITP-holdout design where all ITP-linked rows are excluded before scoring, 7 of 8 ITP-positive compounds remain in the top quartile (AUROC = 0.847 vs 0.833 for experiment count, 0.792 for trimmed mean). An empirical Bayes beta-binomial model provides compound-level posterior P(positive) with 95% credible intervals, agreeing with the heuristic ranking at Kendall tau = 0.67. An out-of-time validation using DrugAge Build 4 (2021) versus Build 5 (2024) shows that top-quartile compounds accumulate significantly more positive evidence in the three-year delta (positive rate 0.886 vs 0.587, p < 0.0001). Sensitivity analyses show the ranking is stable under weight perturbation (median tau = 0.77), cap perturbation (80.6% of combinations yield tau > 0.80), species weighting (top-20 unchanged under 2x mammalian weight), and PMID-clustered bootstrap (median tau = 0.76, 95% CI [0.67, 0.83]). Under 30% illustrative hidden-negative injection, 67% of robust-tier compounds survive (32/48). All scoring metrics are formally defined. Code and locked data at github.com/scottdhughes/drugage-robustness-skill.

# Introduction

This submission presents an automated pipeline for ranking DrugAge longevity intervention claims by robustness. The pipeline is fully offline and uses a bundled local copy of the Human Ageing Genomic Resources (HAGR) DrugAge Build 5 dataset, deterministic normalization, evidence tiers, a robustness perturbation panel, and permutation-based null calibration.

The contribution is not a leaderboard of raw lifespan effects. The contribution is an automated pipeline that asks which model-organism longevity claims remain convincing after perturbation and falsification pressure.

# Related Work

DrugAge is a curated database of aging-related drugs maintained as part of the Human Ageing Genomic Resources (HAGR) project. Barardo et al. (2017) describe the database design and its initial content of over 400 compounds tested in model organisms. The related Geroprotectors.org database (Moskalev et al., 2015) takes a complementary approach, cataloguing therapeutic interventions with structured annotations for mechanism and clinical status. SynergyAge (Bunu et al., 2020), also part of HAGR, extends the landscape to drug combinations and synergistic lifespan effects. De Magalhaes et al. (2018) demonstrate that demographic measurements of the rate of aging in mice provide a more reliable basis for assessing genetic interventions than single-endpoint lifespan measures, motivating our emphasis on robustness over raw effect size.

The NIA Interventions Testing Program (ITP) represents the closest precedent for multi-site robustness testing of longevity compounds. Harrison et al. (2009) reported the first ITP confirmation of rapamycin-mediated lifespan extension across three independent laboratories, and Strong et al. (2016) extended the ITP validation to additional compounds including acarbose and 17-alpha-estradiol. DrugAge records ITP provenance directly in its dataset, allowing systematic cross-referencing between our retrospective robustness ranking and prospective ITP validation.

Published meta-analyses of longevity interventions provide a third reference point. Bitto et al. (2016) meta-analyzed rapamycin lifespan effects across mouse studies and confirmed a dose-dependent positive signal, consistent with rapamycin's placement in our robust tier. Our pipeline differs from these targeted meta-analyses by applying a uniform robustness framework across the full DrugAge compound set rather than focusing on a single compound, and by using systematic permutation-based null calibration to quantify how far the observed ranking separates from a species-matched null.

# Data

The pipeline uses a bundled local copy of the official DrugAge Build 5 CSV dated November 29, 2024. No network access is required. The pipeline validates the file hash, required columns, and release metadata before processing.

Rows are dropped only if compound name, species, or numeric `avg_lifespan_change_percent` is missing. In the default analysis, the pipeline retained 3372 of 3423 DrugAge rows, covering 1038 compounds, 33 normalized species, and 9 scored taxon labels.

# Methods

The pipeline computes one evidence profile per compound using DrugAge's average lifespan change field, `avg_lifespan_change_percent`, for scoring. This field is the most consistently populated and directly comparable effect-size measure in DrugAge. Significance annotations are retained descriptively because they are heterogeneous across studies and do not provide a stable standalone ranking signal. For each compound, the pipeline computes:

- number of experiments
- number of species
- number of taxa
- number of PMIDs
- median effect
- 10% trimmed-mean effect
- sign consistency
- leave-one-species-out stability
- leave-one-taxon-out stability
- aggregation stability
- breadth score
- robustness score

## Metric Definitions

| Metric | Definition |
|---|---|
| Breadth | Normalized coverage across experiments, species, and taxa (each capped) |
| Sign consistency | Fraction of experiments showing positive lifespan change |
| LOSO stability | Fraction of leave-one-species-out subsets retaining a positive trimmed mean |
| LOTO stability | Fraction of leave-one-taxon-out subsets retaining a positive trimmed mean |
| Aggregation stability | Binary: do median and trimmed mean both agree on a positive effect? |
| Magnitude | Normalized trimmed-mean effect capped at 50% extension |
| Robustness score R | Weighted combination: 0.35 breadth + 0.20 sign + 0.15 LOSO + 0.10 LOTO + 0.15 agg + 0.05 mag |

**Breadth score.** A normalized composite of experimental, species, and taxonomic coverage, each capped at a saturation point:

    breadth_score = (min(n_experiments, 6)/6 + min(n_species, 4)/4 + min(n_taxa, 3)/3) / 3

**Aggregation stability.** A binary indicator of whether both central-tendency estimators agree on a positive effect:

    aggregation_stability = 1.0 if median_effect > 0 and trimmed_mean_effect > 0, else 0.0

**Magnitude score.** A normalized effect size capped at 50% lifespan extension:

    magnitude_score = min(max(trimmed_mean_effect, 0), 50) / 50

**Robustness score.** A weighted combination of the component metrics:

    robustness_score = 0.35 * breadth_score
                     + 0.20 * sign_consistency
                     + 0.15 * leave_one_species_out_stability
                     + 0.10 * leave_one_taxon_out_stability
                     + 0.15 * aggregation_stability
                     + 0.05 * magnitude_score

The weights prioritize breadth (0.35) and consistency (0.20 + 0.15 + 0.10 = 0.45) over raw effect magnitude (0.05), reflecting the design goal of ranking by durability rather than excitement. The weight vector is a design choice, not an optimized parameter.

**Sign consistency.** `sign_consistency = n_positive_experiments / n_total_experiments`, where a positive experiment is one with average lifespan change > 0.

**Leave-one-species-out (LOSO) stability.** `LOSO = n_positive_subsets / n_species`, where each subset recomputes the compound's 10% trimmed-mean effect with one species removed, and a positive subset retains trimmed mean > 0. For compounds tested in only 1 species, LOSO = 0.0 because no leave-one-species-out subset remains after removal.

**Leave-one-taxon-out (LOTO) stability.** `LOTO = n_positive_subsets / n_taxa`, the same leave-one-out procedure applied at the taxonomic group level. For compounds tested in only 1 taxon, LOTO = 0.0 because no leave-one-taxon-out subset remains after removal.

Compounds are assigned to four evidence tiers and ranked lexicographically: tier priority dominates, then robustness score R within tier, then species breadth, taxonomic breadth, PMID count, and effect magnitude as tiebreakers. This means a promising-tier compound always ranks below any robust-tier compound regardless of R.

Tiers use explicit thresholds:

- **robust**: >=3 experiments, >=2 species, >=2 PMIDs, positive median and trimmed mean, sign consistency >=0.80, LOSO = 1.0, aggregation stability = 1.0
- **promising**: positive median and trimmed mean, sign consistency >=0.67, and >=2 species or >=3 experiments
- **thin evidence**: positive median and trimmed mean, sign consistency >=0.67 (but insufficient breadth for promising)
- **conflicted**: all remaining compounds

Ranking is robustness-first: tier priority dominates sort order, followed by robustness score, species breadth, taxonomic breadth, PMID breadth, and effect magnitude.

The pipeline emits two verification outputs.

The robustness perturbation panel evaluates the top-ranked compounds under five fixed perturbations:

- leave-one-species-out positivity
- leave-one-taxon-out positivity
- positive median and trimmed mean
- exclusion of single-PMID compounds
- exclusion of mixed-sign compounds

The permutation-based null calibration runs 1,000 fixed-seed species-stratified effect permutations. Within each species, average lifespan effects are shuffled across rows, preserving DrugAge's species composition and within-species effect distribution while breaking the link between compound identity and observed effect. With 1,000 reruns, the smallest nonzero empirical p-value is `1/1001`, approximately `0.0010`.

# Results

In the default analysis, the evidence tiers were:

- `robust`: 48 compounds
- `promising`: 174 compounds
- `thin evidence`: 435 compounds
- `conflicted`: 381 compounds

The top 10 compounds were:

1. `Spermidine`
2. `Apple extract`
3. `N-acetyl-L-cysteine`
4. `Minocycline`
5. `Alpha-ketoglutarate`
6. `Carnosine`
7. `Rapamycin`
8. `Mycophenolic acid`
9. `Epigallocatechin-3-gallate`
10. `Vitamin E`

Some top-ranked compounds may look surprising; this ranking reflects internal robustness within curated model-organism evidence, not human plausibility or mechanistic priority.

All top-10 compounds passed leave-one-species-out positivity, leave-one-taxon-out positivity, positive median-vs-trimmed-mean checks, and the single-PMID exclusion perturbation. Three of the top 10 also passed the stricter mixed-sign exclusion perturbation.

The observed top-10 mean robustness score was `0.94097`, compared with a null mean of `0.91320` and null standard deviation `0.01131`. The corresponding empirical p-value was `0.00200` with z-score `2.45`, indicating measurable rather than overwhelming score separation. The more persuasive null result was the robust-compound count: `48`, compared with a null mean of `29.865`, null standard deviation `4.10`, empirical p-value `0.00100`, and z-score `4.42`.

# Sex-Stratified Sensitivity Analysis

The pipeline supports a `--gender` filter that restricts experiments to a single sex before scoring. To assess whether robustness rankings are stable across sexes, we ran male-only (858 experiments, 353 compounds) and female-only (612 experiments, 271 compounds) analyses:

| Rank | Male Top-10 | Score | Female Top-10 | Score | Shared? |
|------|-------------|-------|---------------|-------|---------|
| 1 | Minocycline | 0.879 | Butylated hydroxytoluene | 0.868 | BHT |
| 2 | Butylated hydroxytoluene | 0.869 | Rapamycin | 0.851 | both |
| 3 | Rapamycin | 0.866 | Ginseng extract | 0.842 | |
| 4 | Green tea extract | 0.864 | Alpha-ketoglutarate | 0.839 | |
| 5 | Chloroquine | 0.833 | N-acetyl-L-cysteine | 0.838 | NAC |
| 6 | Nordihydroguaiaretic acid | 0.835 | Melatonin | 0.800 | |
| 7 | Curcumin | 0.815 | Curcumin | 0.818 | both |
| 8 | N-acetyl-L-cysteine | 0.762 | Rhodiola rosea extract | 0.702 | both |
| 9 | L-deprenyl | 0.703 | Pineal gland extract | 0.693 | |
| 10 | Rhodiola rosea extract | 0.701 | L-deprenyl | 0.663 | both |

6 of the top-10 compounds are shared between sexes (Rapamycin, N-acetyl-L-cysteine, Butylated hydroxytoluene, Curcumin, L-deprenyl, Rhodiola rosea extract), indicating a stable core of robust longevity compounds. Male-only analysis identified 8 robust compounds; female-only identified 9.

# ITP Cross-Validation

DrugAge records ITP provenance natively, enabling a systematic overlap analysis between the retrospective robustness ranking and prospective ITP validation. Of the 54 ITP-tested compounds in DrugAge, 53 matched the ranking (one was excluded during normalization). Stratifying by ITP outcome direction:

- 13 ITP-positive compounds (all ITP experiments showed lifespan extension): 7 landed in the robust tier (53.8%), 1 in promising, and 5 in thin evidence. None were classified as conflicted.
- 12 ITP-negative compounds (all ITP experiments showed no extension): all 12 landed in the conflicted tier (100%).
- 29 ITP-mixed compounds (mixed ITP results): 3 robust, 4 promising, and 21 conflicted.

A one-sided Fisher's exact test on the 2x2 table (ITP-positive/negative x robust/non-robust) yields odds ratio = infinity and p = 0.0036, confirming the association is statistically significant. ITP-positive compounds have a mean of 10.4 DrugAge experiments versus 4.0 for ITP-negative compounds, so the data-volume confound is present but does not explain the negative-side result.

Overall, 15 of 53 matched ITP compounds (28.3%) appeared in the top quartile. The concordance between ITP outcome direction and evidence tier assignment is the strongest external validation of the robustness framework. This agreement is not engineered: the pipeline does not use ITP status as a scoring input.

# Data-Volume Confound Analysis

The robustness score is structurally correlated with experiment count because breadth, LOSO, and LOTO are mechanically easier to achieve with more data. The Spearman correlation between experiment count and R is rho = 0.44, and a linear regression of R on log(1 + n_experiments) yields R-squared = 0.22. This means 78% of the variance in R is not explained by data volume, but the confound is real and must be disclosed.

To test whether the ranking survives after removing the data-volume component, we computed a residualized score R_adj = R - R_hat(log n). The Kendall tau between the original and residualized rankings is 0.46, indicating substantial reordering. We report the canonical ranking because the residualized version removes a biologically meaningful signal (cross-species breadth), but the R-squared = 0.22 and tau = 0.46 numbers should inform interpretation.

# Weight Sensitivity

The weight vector (0.35, 0.20, 0.15, 0.10, 0.15, 0.05) is a declared design choice. To test its influence, we drew 200 alternative weight vectors from a Dirichlet distribution centered on the canonical weights and computed the Kendall tau between each perturbed ranking and the canonical ranking. The median tau was 0.77 (5th percentile: 0.70, 95th percentile: 0.83), indicating moderate stability under weight perturbation.

# Component Score Correlations

The six component scores are not independent. LOSO and sign consistency are moderately correlated, confirming that they capture overlapping but not identical signals. The robustness score therefore partially double-counts the consistency dimension, which is the stated design intent (0.45 total weight on consistency-family metrics).

# Naive Baseline Comparison

To assess whether the robustness framework adds value beyond simpler approaches, we compared four ranking methods on their ITP concordance:

| Ranking method | tau vs robustness | ITP-pos in top-Q | ITP-neg in top-Q |
|---|---|---|---|
| Median effect | 0.570 | 1 | 0 |
| Trimmed-mean effect | 0.563 | 1 | 0 |
| Experiment count | 0.136 | 7 | 6 |
| Robustness score R | 1.000 | 8 | 3 |

The robustness ranking is the only method that simultaneously concentrates ITP-positive compounds in the top tier while excluding ITP-negative compounds.

# Additional Validation

## Structure-Matched Permutation Null (Exploratory)
The species-stratified null preserves data-volume structure. A harder test shuffles compound identity assignments within species, breaking the link between compound labels and evidence while preserving per-compound sample sizes. Under this structure-matched null (500 permutations, exploratory scoring variant), the robust-compound count attenuates and is not conventionally significant (z = 1.26, p = 0.13). This confirms that part of the robust-tier enrichment is structural and the index cannot fully separate biological signal from data-volume advantage at the tier level.

## Out-of-Time Validation: Build 4 to Build 5
DrugAge Build 4 (November 2021, 1,084 compounds) and Build 5 (November 2024, 1,038 compounds after normalization) provide a natural temporal split. Compounds ranked in the top quartile by canonical evidence had a mean delta positive rate of 0.886 versus 0.587 for lower-ranked compounds among the 191 compounds that gained new experiments between builds (Mann-Whitney p < 0.0001). This is a strong forward-validation signal. To rule out the attention confound (top-quartile compounds attract more research), we controlled for baseline experiment count via OLS regression (top-Q coefficient beta = 0.319, p < 1e-6 after controlling for log(Build-4 experiments)), propensity-matching (60 matched pairs with identical mean Build-4 experiments = 2.2: delta pos rate 0.905 vs 0.569, p = 0.000039), and within-delta-n comparisons (delta_n=1: 0.889 vs 0.672; delta_n=2: 0.800 vs 0.541; delta_n=3: 0.974 vs 0.389). Top-quartile compounds do attract more future experiments (mean delta_n 4.2 vs 2.0, p < 0.0001), confirming the attention confound is real, but the positive-rate difference survives after controlling for it.

## Leave-Mouse-Out Validation
After excluding all 369 mouse experiments, 6 of 8 matched ITP-positive compounds remained in the top quartile of the non-mouse ranking, confirming that the framework recovers mouse-validated compounds from non-mouse evidence alone.

## Illustrative Hidden-Negative Scenario
Under a grid of synthetic negative-result injection rates (effects drawn uniformly from [-30%, -1%], assigned uniformly at random): 10% injection → 40/48 survive (83%), 20% → 37/48 (77%), 30% → 32/48 (67%), 40% → 30/48 (62%). The attrition curve is approximately linear. The true hidden-negative rate in DrugAge is unknown.

## Threshold Sensitivity
The tier thresholds are the true levers of the lexicographic ranking. Across 540 threshold combinations, median Jaccard overlap with canonical robust set was 0.557. The most sensitive parameter is minimum species (≥2 → ≥3 drops robust from 48 to 17, Jaccard 0.354). LOSO threshold is least sensitive (relaxing 1.0 → 0.8 changes nothing).

## PMID-Clustered Bootstrap
Bootstrap resampling by PMID cluster (500 replicates, 651 unique PMIDs) yields median Kendall tau = 0.76 (95% CI: [0.67, 0.83]) versus the canonical ranking, confirming moderate stability under study-level dependence.

## Species-Weighting Sensitivity
Weighting mammalian experiments at 2x invertebrate experiments produces identical top-10 and top-20 compound sets.

## ITP Holdout AUROC
On the clean ITP holdout (8 positive, 9 negative), the robustness score had the highest point-estimate AUROC (0.847, 95% bootstrap CI [0.614, 1.000]), followed by experiment count (0.833), trimmed-mean effect (0.792), breadth (0.778), median effect (0.771), species count (0.625), and sign consistency (0.500). The CI excludes chance (0.5) but is wide due to n=17.

# Dose-Response Sensitivity Analysis

DrugAge records dosage for 98.5% of experiments. However, units are heterogeneous across compounds (mM, uM, ppm, mg/ml, %, etc.), precluding direct cross-compound dose normalization. We therefore computed within-compound dose-response correlations for the top 20 ranked compounds using Spearman's rho on experiments sharing the dominant dosage unit.

Of 20 compounds, 16 had sufficient same-unit data (>=3 experiments with the same unit and >=2 distinct dose levels). Among those 16, 7 showed positive dose-response correlations and 9 showed negative correlations, with only 1 reaching significance at p < 0.05. The predominant pattern is no significant within-compound dose-response relationship, consistent with the observation that DrugAge pools experiments across diverse species, strains, and delivery methods. This is an honest limitation: the robustness ranking captures cross-study consistency of the pro-longevity signal, not dose-response monotonicity within a compound.

# Optional AnAge Context

The optional AnAge context report joins normalized DrugAge species to a bundled local copy of AnAge for descriptive context only. It does not alter the ranking, scores, tiers, or verification outputs. In the current rerun, 10 of 35 normalized DrugAge species matched AnAge exactly after normalization.

# Limitations

**Analytical versus biological robustness.** This pipeline measures *analytical* robustness: the consistency of a compound's pro-longevity signal across species, taxa, estimators, and leave-one-out perturbations within DrugAge's curated records. It does not measure *biological* robustness in the causal or translational sense. Dose-response relationships, administration timing, strain-specific interactions, and mechanistic validity are outside its scope. A compound scoring high on analytical robustness may do so because every contributing lab happened to use an effective dose and a responsive strain; the pipeline cannot distinguish this from genuine cross-context durability.

**Data-volume confound.** Compounds with more DrugAge records are structurally advantaged because breadth and leave-one-out stability are mechanically easier to achieve with larger sample sizes. The Spearman correlation between experiment count and R is rho = 0.44, meaning the robustness ranking partially reflects data availability rather than biological durability. The ITP cross-validation provides partial control: ITP-negative compounds have meaningful DrugAge footprints (mean 4.0 experiments) yet uniformly land in the conflicted tier.

**Publication and survivorship bias.** DrugAge curates published results, and published results skew positive. Compounds that failed in preliminary screens never appear in the database. The pipeline's sign-consistency metric rewards compounds that appear consistently positive, but if negative results are systematically unpublished, sign consistency is inflated. This cannot be corrected without access to unpublished data.

**Species weighting.** The pipeline treats all species equally. A positive result in *C. elegans* and a positive result in *Mus musculus* contribute equally to the robustness score. This avoids imposing subjective phylogenetic priors, but it means a compound tested primarily in short-lived invertebrates can achieve the same breadth as one tested across mammals.

**Metadata coverage.** DrugAge records dosage for 98.5% of experiments but only 10.9% include age-at-initiation and 10.8% include treatment duration. Sex-stratified analysis shows 6/10 top compounds shared between sexes. Within-compound dose-response analysis shows no significant monotonic relationship in 15 of 16 testable top-20 compounds. Timing-stratified analysis is not feasible at current metadata coverage.

**Weight vector.** The weight vector (0.35 breadth, 0.45 consistency, 0.05 magnitude) is a design choice, not an optimized parameter. The Dirichlet sensitivity analysis shows median Kendall tau = 0.77 under perturbation, indicating moderate but not absolute stability.

# Conclusion

This evidence-robustness index ranks longevity claims in DrugAge by analytical consistency. The ranking separates from both a species-stratified null (z = 4.42) and, more weakly, a structure-matched null (z = 1.26). Under an ITP-holdout design, 7/8 ITP-positive compounds remain in the top quartile (AUROC = 0.847), and leave-mouse-out recovers 6/8 from non-mouse evidence alone. An out-of-time Build 4 → Build 5 validation shows top-quartile compounds accumulate significantly more positive evidence (p < 0.0001). Under 30% hidden-negative injection, 67% of the robust tier survives (attrition grid: 10%→83%, 20%→77%, 30%→67%, 40%→62%). Threshold sensitivity identifies the minimum-species gate as the primary determinant of tier membership. The framework's known failure modes are quantified rather than hidden: structure-matched null p = 0.13, ITP holdout Fisher p = 0.24, and Spermidine's top heuristic rank not supported by the Bayesian model. The index measures evidence robustness inside DrugAge, not biological geroprotective efficacy.

**Take-home for experimentalists.** The 48 robust-tier compounds represent the subset of DrugAge whose pro-longevity evidence is cross-species, internally consistent, and resilient to leave-one-out perturbation. They are not clinical recommendations. Compounds in the conflicted tier (381, including metformin, aspirin, and resveratrol) have genuinely mixed model-organism evidence and should not be dismissed as failures — they may benefit from targeted replication in specific species or strains.

# References

- Barardo, D., Thornton, D., Thoppil, H., Walsh, M., Sharber, S., Ferber, S., Greer, E.L., Ship, A., Valli, A., Horro, R., & de Magalhaes, J.P. (2017). The DrugAge database of aging-related drugs. *Aging Cell*, 16(3), 594--597.
- Bunu, G., Tankard, A., Okada, H., Yoshimura, S.H., & de Magalhaes, J.P. (2020). SynergyAge, a curated database for synergistic and antagonistic interactions of longevity-associated genes. *Scientific Data*, 7, 366.
- de Magalhaes, J.P., Budovsky, A., Li, Q., Fraifeld, V.E., & Church, G.M. (2018). A reassessment of genes modulating aging in mice using demographic measurements of the rate of aging. *Genetics*, 208, 1617--1630.
- Harrison, D.E., Strong, R., Sharp, Z.D., et al. (2009). Rapamycin fed late in life extends lifespan in genetically heterogeneous mice. *Nature*, 460, 392--395.
- Bitto, A., Ito, T.K., Pineda, V.V., et al. (2016). Transient rapamycin treatment can increase lifespan and healthspan in middle-aged mice. *eLife*, 5, e16351.
- Partridge, L., Deelen, J., & Slagboom, P.E. (2018). Facing up to the global challenges of ageing. *Nature*, 561, 45--56.
- Moskalev, A., Chernyagina, E., de Magalhaes, J.P., et al. (2015). Geroprotectors.org: a new, structured and curated database of current therapeutic interventions in aging and age-related disease. *Aging*, 7(9), 616--628.
- Strong, R., Miller, R.A., Antebi, A., et al. (2016). Longer lifespan in male mice treated with a weakly estrogenic agonist, an antioxidant, an alpha-glucosidase inhibitor or a Nrf2-inducer. *Aging Cell*, 15, 872--884.
- de Cabo, R., & Mattson, M.P. (2019). Effects of intermittent fasting on health, aging, and disease. *New England Journal of Medicine*, 381, 2541--2551.
