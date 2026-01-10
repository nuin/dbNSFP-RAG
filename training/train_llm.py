#!/usr/bin/env python
"""
Fine-tune LLM for variant interpretation.

Supports:
- Llama 3.2 with LoRA (recommended for M1/M2/M3 Macs)
- Full fine-tuning on larger hardware

Usage:
    python training/train_llm.py --panel hereditary_cancer --method lora
    python training/train_llm.py --input data/exports/custom_llm_training.jsonl
"""

import json
import argparse
from pathlib import Path


def prepare_dataset(filepath: Path, output_path: Path):
    """Convert to Alpaca format for training."""
    samples = []
    with open(filepath) as f:
        for line in f:
            data = json.loads(line)
            samples.append({
                "instruction": data["instruction"],
                "input": data.get("input", ""),
                "output": data["output"],
            })

    with open(output_path, "w") as f:
        json.dump(samples, f, indent=2)

    print(f"Prepared {len(samples)} samples -> {output_path}")
    return output_path


def train_lora(data_path: Path, base_model: str = "meta-llama/Llama-3.2-3B-Instruct"):
    """
    Fine-tune with LoRA using mlx-lm (Apple Silicon optimized).

    Install: pip install mlx-lm
    """
    print(f"""
=== LoRA Fine-tuning on Apple Silicon ===

1. Install mlx-lm:
   pip install mlx-lm

2. Convert training data:
   python training/train_llm.py --panel your_panel --prepare-only

3. Run LoRA training:
   mlx_lm.lora \\
       --model {base_model} \\
       --data {data_path.parent} \\
       --train \\
       --batch-size 4 \\
       --lora-layers 16 \\
       --iters 1000

4. Merge adapters:
   mlx_lm.fuse \\
       --model {base_model} \\
       --adapter-path adapters \\
       --save-path models/variant-llm

5. Test:
   mlx_lm.generate \\
       --model models/variant-llm \\
       --prompt "Interpret variant BRCA1 c.5266dupC"
""")


def train_full_finetune(data_path: Path):
    """
    Full fine-tuning using Hugging Face Transformers.

    For larger hardware (multi-GPU).
    """
    print(f"""
=== Full Fine-tuning with Transformers ===

Install:
    pip install transformers datasets accelerate bitsandbytes

Training script:
```python
from datasets import load_dataset
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    TrainingArguments,
    Trainer,
)

# Load model
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-3.2-3B-Instruct",
    torch_dtype="auto",
    device_map="auto",
)
tokenizer = AutoTokenizer.from_pretrained("meta-llama/Llama-3.2-3B-Instruct")

# Load data
dataset = load_dataset("json", data_files="{data_path}")

# Training
trainer = Trainer(
    model=model,
    train_dataset=dataset["train"],
    args=TrainingArguments(
        output_dir="models/variant-llm",
        num_train_epochs=3,
        per_device_train_batch_size=4,
        gradient_accumulation_steps=4,
        learning_rate=2e-5,
        fp16=True,
        save_steps=500,
    ),
)
trainer.train()
```
""")


def main():
    parser = argparse.ArgumentParser(description="Fine-tune LLM for variant interpretation")
    parser.add_argument("--panel", "-p", help="Panel name")
    parser.add_argument("--input", "-i", help="Pre-exported JSONL file")
    parser.add_argument("--method", "-m", choices=["lora", "full"], default="lora")
    parser.add_argument("--prepare-only", action="store_true", help="Only prepare data, don't train")
    args = parser.parse_args()

    # Get training data
    if args.input:
        data_path = Path(args.input)
    elif args.panel:
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from src.export import export_for_llm_finetuning
        data_path = export_for_llm_finetuning(args.panel)
    else:
        print("Specify --panel or --input")
        return

    # Prepare Alpaca format
    alpaca_path = data_path.with_suffix(".alpaca.json")
    prepare_dataset(data_path, alpaca_path)

    if args.prepare_only:
        return

    # Show training instructions
    if args.method == "lora":
        train_lora(alpaca_path)
    else:
        train_full_finetune(alpaca_path)


if __name__ == "__main__":
    main()
