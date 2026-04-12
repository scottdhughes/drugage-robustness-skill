"""Run the pipeline with 1000 permutations instead of 128."""
import yaml, subprocess, sys, json

# Load config, change null_reruns
with open('config/canonical_drugage.yaml') as f:
    cfg = yaml.safe_load(f)

cfg['pipeline']['null_reruns'] = 1000
with open('config/null_1000.yaml', 'w') as f:
    yaml.dump(cfg, f)

# Run pipeline with new config
result = subprocess.run(
    [sys.executable, '-m', 'drugage_skill.pipeline', '--config', 'config/null_1000.yaml', '--outdir', 'outputs/null_1000'],
    capture_output=True, text=True
)
print(result.stdout[-500:] if result.stdout else "")
if result.returncode != 0:
    print("STDERR:", result.stderr[-500:])

# Read results
with open('outputs/null_1000/empirical_null_certificate.json') as f:
    cert = json.load(f)
print(f"\n=== 1000 PERMUTATION RESULTS ===")
print(f"Observed top-10 mean: {cert.get('observed_top10_mean', '?')}")
print(f"Null mean: {cert.get('null_mean_top10', '?')}")
print(f"p-value (top-10): {cert.get('p_value_top10', '?')}")
print(f"Observed robust count: {cert.get('observed_robust_count', '?')}")
print(f"Null mean robust: {cert.get('null_mean_robust_count', '?')}")
print(f"p-value (robust): {cert.get('p_value_robust_count', '?')}")
