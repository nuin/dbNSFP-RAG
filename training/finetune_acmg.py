#!/usr/bin/env python
"""
Fine-tune LLM for ACMG variant classification.

Supports:
- Apple Silicon via mlx-lm (recommended for M1/M2/M3/M4)
- NVIDIA GPU via transformers + LoRA

Usage:
    # Generate training data first
    python training/acmg_training.py

    # Fine-tune with mlx-lm (Apple Silicon)
    python training/finetune_acmg.py --method mlx

    # Fine-tune with transformers (NVIDIA)
    python training/finetune_acmg.py --method transformers
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def prepare_mlx_data(input_path: Path, output_dir: Path):
    """Prepare data in mlx-lm format (train.jsonl, valid.jsonl)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    with open(input_path) as f:
        data = json.load(f)

    # Split 90/10
    split_idx = int(len(data) * 0.9)
    train_data = data[:split_idx]
    valid_data = data[split_idx:]

    # Convert to chat format for mlx-lm
    def to_chat(example):
        return {
            "messages": [
                {"role": "user", "content": f"{example['instruction']}\n\n{example['input']}"},
                {"role": "assistant", "content": example["output"]}
            ]
        }

    # Save
    train_path = output_dir / "train.jsonl"
    valid_path = output_dir / "valid.jsonl"

    with open(train_path, "w") as f:
        for ex in train_data:
            f.write(json.dumps(to_chat(ex)) + "\n")

    with open(valid_path, "w") as f:
        for ex in valid_data:
            f.write(json.dumps(to_chat(ex)) + "\n")

    print(f"Prepared {len(train_data)} training, {len(valid_data)} validation examples")
    return output_dir


def finetune_mlx(data_dir: Path, output_dir: Path, base_model: str, iters: int):
    """Fine-tune using mlx-lm on Apple Silicon."""
    print(f"\n{'='*60}")
    print("Fine-tuning with mlx-lm (Apple Silicon)")
    print(f"{'='*60}\n")

    # Check mlx-lm is installed
    try:
        import mlx_lm
    except ImportError:
        print("Installing mlx-lm...")
        subprocess.run([sys.executable, "-m", "pip", "install", "mlx-lm"], check=True)

    # Run LoRA training (new mlx-lm API)
    cmd = [
        sys.executable, "-m", "mlx_lm",
        "lora",
        "--model", base_model,
        "--data", str(data_dir),
        "--train",
        "--batch-size", "4",
        "--num-layers", "16",
        "--iters", str(iters),
        "--adapter-path", str(output_dir / "adapters"),
        "--steps-per-report", "50",
    ]

    print(f"Running: {' '.join(cmd)}\n")
    subprocess.run(cmd, check=True)

    # Fuse adapters into model
    print("\nFusing adapters into model...")
    fuse_cmd = [
        sys.executable, "-m", "mlx_lm",
        "fuse",
        "--model", base_model,
        "--adapter-path", str(output_dir / "adapters"),
        "--save-path", str(output_dir / "model"),
    ]
    subprocess.run(fuse_cmd, check=True)

    print(f"\nModel saved to: {output_dir / 'model'}")
    return output_dir / "model"


def finetune_transformers(data_path: Path, output_dir: Path, base_model: str):
    """Fine-tune using transformers + LoRA."""
    print(f"\n{'='*60}")
    print("Fine-tuning with Transformers + LoRA")
    print(f"{'='*60}\n")

    try:
        from datasets import load_dataset
        from transformers import (
            AutoModelForCausalLM,
            AutoTokenizer,
            TrainingArguments,
            Trainer,
        )
        from peft import LoraConfig, get_peft_model
        import torch
    except ImportError:
        print("Installing required packages...")
        subprocess.run([
            sys.executable, "-m", "pip", "install",
            "transformers", "datasets", "peft", "accelerate", "bitsandbytes"
        ], check=True)
        from datasets import load_dataset
        from transformers import (
            AutoModelForCausalLM,
            AutoTokenizer,
            TrainingArguments,
            Trainer,
        )
        from peft import LoraConfig, get_peft_model
        import torch

    # Load tokenizer and model
    print(f"Loading {base_model}...")
    tokenizer = AutoTokenizer.from_pretrained(base_model)
    tokenizer.pad_token = tokenizer.eos_token

    model = AutoModelForCausalLM.from_pretrained(
        base_model,
        torch_dtype=torch.float16,
        device_map="auto",
    )

    # Configure LoRA
    lora_config = LoraConfig(
        r=16,
        lora_alpha=32,
        target_modules=["q_proj", "k_proj", "v_proj", "o_proj"],
        lora_dropout=0.05,
        bias="none",
        task_type="CAUSAL_LM",
    )

    model = get_peft_model(model, lora_config)
    model.print_trainable_parameters()

    # Load dataset
    dataset = load_dataset("json", data_files=str(data_path))

    def format_example(example):
        text = f"### Instruction:\n{example['instruction']}\n\n### Input:\n{example['input']}\n\n### Response:\n{example['output']}"
        return {"text": text}

    dataset = dataset.map(format_example)

    def tokenize(example):
        return tokenizer(
            example["text"],
            truncation=True,
            max_length=1024,
            padding="max_length",
        )

    tokenized = dataset.map(tokenize, remove_columns=dataset["train"].column_names)

    # Training arguments
    training_args = TrainingArguments(
        output_dir=str(output_dir),
        num_train_epochs=3,
        per_device_train_batch_size=4,
        gradient_accumulation_steps=4,
        learning_rate=2e-4,
        fp16=True,
        save_steps=500,
        logging_steps=100,
        warmup_steps=100,
    )

    # Train
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized["train"],
    )

    trainer.train()

    # Save
    model.save_pretrained(output_dir / "model")
    tokenizer.save_pretrained(output_dir / "model")

    print(f"\nModel saved to: {output_dir / 'model'}")
    return output_dir / "model"


def test_model(model_path: Path, method: str):
    """Test the fine-tuned model."""
    print(f"\n{'='*60}")
    print("Testing fine-tuned model")
    print(f"{'='*60}\n")

    test_input = """Variant: chr17:41197801 T>A
Gene: BRCA1
CADD phred: 19.09"""

    prompt = f"Classify this variant according to ACMG/AMP guidelines and provide the evidence criteria.\n\n{test_input}"

    if method == "mlx":
        from mlx_lm import load, generate
        model, tokenizer = load(str(model_path))
        response = generate(model, tokenizer, prompt=prompt, max_tokens=300)
    else:
        from transformers import AutoModelForCausalLM, AutoTokenizer
        import torch

        tokenizer = AutoTokenizer.from_pretrained(model_path)
        model = AutoModelForCausalLM.from_pretrained(model_path, device_map="auto")

        inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
        outputs = model.generate(**inputs, max_new_tokens=300, temperature=0.7)
        response = tokenizer.decode(outputs[0], skip_special_tokens=True)

    print(f"Input:\n{test_input}\n")
    print(f"Response:\n{response}")


def main():
    parser = argparse.ArgumentParser(description="Fine-tune LLM for ACMG classification")
    parser.add_argument(
        "--method", "-m",
        choices=["mlx", "transformers"],
        default="mlx",
        help="Fine-tuning method (mlx for Apple Silicon, transformers for NVIDIA)"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        default=Path("data/training/acmg_training.json"),
        help="Input training data (Alpaca JSON format)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("models/acmg-classifier"),
        help="Output directory"
    )
    parser.add_argument(
        "--base-model",
        default="mlx-community/Llama-3.2-3B-Instruct-4bit",
        help="Base model to fine-tune"
    )
    parser.add_argument(
        "--iters",
        type=int,
        default=1000,
        help="Training iterations (mlx only)"
    )
    parser.add_argument(
        "--test-only",
        action="store_true",
        help="Only test existing model"
    )
    args = parser.parse_args()

    if args.test_only:
        test_model(args.output / "model", args.method)
        return

    if not args.input.exists():
        print(f"Training data not found: {args.input}")
        print("Run: python training/acmg_training.py")
        return

    args.output.mkdir(parents=True, exist_ok=True)

    if args.method == "mlx":
        # Prepare data for mlx-lm
        data_dir = prepare_mlx_data(args.input, args.output / "data")
        model_path = finetune_mlx(data_dir, args.output, args.base_model, args.iters)
    else:
        model_path = finetune_transformers(args.input, args.output, args.base_model)

    # Test the model
    test_model(model_path, args.method)

    print(f"\n{'='*60}")
    print("Fine-tuning complete!")
    print(f"{'='*60}")
    print(f"\nModel saved to: {model_path}")
    print(f"\nTo use with API:")
    print(f"  export ACMG_MODEL_PATH={model_path}")
    print(f"  uvicorn api.server:app --host 0.0.0.0 --port 8000")


if __name__ == "__main__":
    main()
