"""Sex-stratified sensitivity analysis: run pipeline on Male-only and Female-only subsets."""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from drugage_skill.pipeline import run_pipeline
import yaml

with open('config/canonical_drugage.yaml') as f:
    cfg = yaml.safe_load(f)

# Read full data
with open(cfg['dataset']['path']) as f:
    all_rows = list(csv.DictReader(f))
    headers = list(all_rows[0].keys())

results = {}
for sex in ['Male', 'Female']:
    subset = [r for r in all_rows if r['gender'] == sex]
    
    # Write subset to temp file
    temp_path = f'data/drugage_{sex.lower()}_only.csv'
    with open(temp_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        w.writerows(subset)
    
    # Update config for subset
    cfg_sex = dict(cfg)
    cfg_sex['dataset'] = dict(cfg['dataset'])
    cfg_sex['dataset']['path'] = temp_path
    # Remove SHA check for subset
    cfg_sex['dataset'].pop('sha256', None)
    cfg_sex['pipeline']['null_reruns'] = 32  # fewer perms for speed
    
    config_path = f'config/sex_{sex.lower()}.yaml'
    with open(config_path, 'w') as f:
        yaml.dump(cfg_sex, f)
    
    outdir = f'outputs/sex_{sex.lower()}'
    os.makedirs(outdir, exist_ok=True)
    
    try:
        import subprocess
        result = subprocess.run(
            [sys.executable, '-m', 'drugage_skill.pipeline', '--config', config_path, '--outdir', outdir],
            capture_output=True, text=True, timeout=120
        )
        
        # Read top 10
        ranking_path = os.path.join(outdir, 'robustness_rankings.csv')
        if os.path.exists(ranking_path):
            with open(ranking_path) as f:
                ranked = list(csv.DictReader(f))
            top10 = [r['compound_name'] for r in ranked[:10]]
            n_robust = sum(1 for r in ranked if r.get('evidence_tier') == 'robust')
            results[sex] = {'n_rows': len(subset), 'top10': top10, 'n_robust': n_robust, 'total': len(ranked)}
            print(f"\n{sex}: {len(subset)} experiments, {len(ranked)} compounds, {n_robust} robust")
            print(f"  Top 10: {', '.join(top10)}")
        else:
            print(f"{sex}: pipeline failed")
            results[sex] = {'error': result.stderr[:200]}
    except Exception as e:
        print(f"{sex}: {e}")
        results[sex] = {'error': str(e)}

# Compare
if 'top10' in results.get('Male', {}) and 'top10' in results.get('Female', {}):
    m_set = set(results['Male']['top10'])
    f_set = set(results['Female']['top10'])
    overlap = m_set & f_set
    print(f"\n=== COMPARISON ===")
    print(f"Male top-10: {results['Male']['top10']}")
    print(f"Female top-10: {results['Female']['top10']}")
    print(f"Overlap: {len(overlap)}/10 — {overlap}")
    print(f"Male-only: {m_set - f_set}")
    print(f"Female-only: {f_set - m_set}")
    print(f"Male robust: {results['Male']['n_robust']}, Female robust: {results['Female']['n_robust']}")

with open('outputs/sex_stratified_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to outputs/sex_stratified_results.json")
