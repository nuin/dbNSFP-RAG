"""Load optional language-model backends for legacy inference callers."""

import logging
import os
from pathlib import Path

from api.inference_generation import generate_classification, parse_classification_response

__all__ = ["generate_classification", "load_model", "parse_classification_response"]

logger = logging.getLogger(__name__)


def _load_gguf(model_path: str) -> dict | None:
    try:
        from llama_cpp import Llama

        n_ctx = int(os.environ.get("LLM_N_CTX", "4096"))
        n_threads = int(os.environ.get("LLM_N_THREADS", "4"))
        n_gpu_layers = int(os.environ.get("LLM_N_GPU_LAYERS", "0"))
        llm = Llama(
            model_path=model_path,
            n_ctx=n_ctx,
            n_threads=n_threads,
            n_gpu_layers=n_gpu_layers,
            verbose=False,
        )
    except (ImportError, OSError, RuntimeError, ValueError, TypeError) as exc:
        logger.warning("GGUF backend unavailable; trying the next backend: %s", exc)
        return None
    logger.info("Loaded llama.cpp model (n_ctx=%s, threads=%s)", n_ctx, n_threads)
    return {"type": "llama_cpp", "model": llm}


def _load_mlx(model_path: str) -> dict | None:
    try:
        from mlx_lm import load

        model, tokenizer = load(model_path)
    except (ImportError, OSError, RuntimeError, ValueError, TypeError) as exc:
        logger.warning("MLX backend unavailable; trying Transformers: %s", exc)
        return None
    logger.info("Loaded MLX model from %s", model_path)
    return {"type": "mlx", "model": model, "tokenizer": tokenizer}


def _load_transformers(model_path: str) -> dict | None:
    try:
        import torch
        from transformers import AutoModelForCausalLM, AutoTokenizer

        device = "cuda" if torch.cuda.is_available() else "cpu"
        tokenizer = AutoTokenizer.from_pretrained(model_path)
        model = AutoModelForCausalLM.from_pretrained(
            model_path,
            torch_dtype=torch.float16 if device == "cuda" else torch.float32,
            device_map="auto" if device == "cuda" else None,
        )
    except (ImportError, OSError, RuntimeError, ValueError, TypeError) as exc:
        logger.warning("Transformers backend unavailable; using Ollama: %s", exc)
        return None
    logger.info("Loaded Transformers model from %s", model_path)
    return {"type": "transformers", "model": model, "tokenizer": tokenizer, "device": device}


def load_model() -> dict:
    """Try the configured local backends before the external Ollama service."""
    model_path = os.environ.get("ACMG_MODEL_PATH", "models/acmg-classifier/model")
    gguf_path = os.environ.get("LLM_MODEL_PATH", "/app/models/model.gguf")
    use_bundled = os.environ.get("USE_BUNDLED_LLM", "false").lower() == "true"
    if use_bundled or Path(gguf_path).exists():
        model = _load_gguf(gguf_path)
        if model is not None:
            return model
    for loader in (_load_mlx, _load_transformers):
        model = loader(model_path)
        if model is not None:
            return model
    return {"type": "ollama", "model": "llama3.2:3b"}
