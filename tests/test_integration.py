from pathlib import Path

from drugage_skill.anage import run_anage_context_report
from drugage_skill.config import load_config
from drugage_skill.constants import OPTIONAL_ANAGE_OUTPUTS, REQUIRED_CANONICAL_OUTPUTS
from drugage_skill.pipeline import run_pipeline
from drugage_skill.verify import verify_run


def test_full_canonical_run_and_verify(tmp_path: Path) -> None:
    config = load_config("config/canonical_drugage.yaml")
    run_dir = tmp_path / "canonical"
    manifest = run_pipeline(config, run_dir)
    assert manifest["evidence_tiers"]["robust"] > 0
    for name in REQUIRED_CANONICAL_OUTPUTS:
        path = run_dir / config.outputs[name]
        assert path.exists()
        assert path.stat().st_size > 0
    verification = verify_run(config, run_dir)
    assert verification["status"] == "passed"
    run_anage_context_report(config, run_dir)
    for name in OPTIONAL_ANAGE_OUTPUTS:
        path = run_dir / config.outputs[name]
        assert path.exists()
        assert path.stat().st_size > 0
