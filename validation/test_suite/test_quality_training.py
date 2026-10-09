"""Training setup checks that do not load or download model weights."""

import json

import pytest

from training import acmg_training, finetune_acmg


@pytest.mark.parametrize("method,extra", [("mlx", "mlx"), ("transformers", "training")])
def test_missing_training_dependencies_require_explicit_install(
    monkeypatch, tmp_path, method, extra
):
    monkeypatch.setattr(finetune_acmg, "find_spec", lambda package: None)
    commands = []
    monkeypatch.setattr(
        finetune_acmg.subprocess, "run", lambda *args, **kwargs: commands.append(args)
    )
    with pytest.raises(RuntimeError, match=f"uv sync --extra {extra}"):
        if method == "mlx":
            finetune_acmg.finetune_mlx(tmp_path, tmp_path, "test-model", 10)
        else:
            finetune_acmg.finetune_transformers(tmp_path / "training.json", tmp_path, "test-model")
    assert commands == []


def test_mlx_training_and_fusion_keep_explicit_model_paths(monkeypatch, tmp_path):
    monkeypatch.setattr(finetune_acmg, "find_spec", lambda package: object())
    commands = []
    monkeypatch.setattr(
        finetune_acmg.subprocess,
        "run",
        lambda command, **kwargs: commands.append((command, kwargs)),
    )
    model_path = finetune_acmg.finetune_mlx(
        tmp_path / "data", tmp_path / "output", "base/model", 12
    )
    assert model_path == tmp_path / "output" / "model"
    training, fusion = [command for command, _ in commands]
    assert training[training.index("--model") + 1] == "base/model"
    assert training[training.index("--iters") + 1] == "12"
    assert fusion[fusion.index("--save-path") + 1] == str(model_path)
    assert all(options == {"check": True} for _, options in commands)


def test_mlx_data_split_preserves_examples_and_chat_roles(tmp_path):
    source = tmp_path / "input.json"
    examples = [
        {"instruction": "Classify", "input": f"variant {number}", "output": f"class {number}"}
        for number in range(10)
    ]
    source.write_text(json.dumps(examples))
    destination = finetune_acmg.prepare_mlx_data(source, tmp_path / "prepared")
    training = [json.loads(line) for line in (destination / "train.jsonl").read_text().splitlines()]
    validation = [
        json.loads(line) for line in (destination / "valid.jsonl").read_text().splitlines()
    ]
    assert len(training) == 9
    assert len(validation) == 1
    messages = validation[0]["messages"]
    assert messages == [
        {"role": "user", "content": "Classify\n\nvariant 9"},
        {"role": "assistant", "content": "class 9"},
    ]


def test_training_output_preserves_uncertain_class_and_missing_evidence():
    output = acmg_training.format_acmg_output("Uncertain_significance", [], "BRCA1")
    assert output.startswith("ACMG Classification: Uncertain significance\n")
    assert "- Limited computational evidence available" in output
    assert output.endswith(
        "Conclusion: This variant in BRCA1 has insufficient or conflicting evidence for definitive classification."
    )
