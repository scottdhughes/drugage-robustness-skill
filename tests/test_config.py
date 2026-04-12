from drugage_skill.config import load_config


def test_load_config_has_expected_contract() -> None:
    config = load_config("config/canonical_drugage.yaml")
    assert config.title.startswith("Claim-Certified")
    assert config.dataset["build"] == "Build 5"
    assert config.pipeline["null_reruns"] == 1000
    assert "robustness_rankings" in config.outputs

