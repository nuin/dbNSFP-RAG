"""Generate and parse responses from optional inference backends."""

import os


def generate_classification(variant_input: str, model_info: dict) -> str:
    """Generate ACMG classification using the model."""
    # Instruction matching training format
    instruction = (
        "Classify this variant according to ACMG/AMP guidelines and provide the evidence criteria."
    )

    if model_info["type"] == "llama_cpp":
        # llama-cpp-python backend (GGUF model for Docker deployment)
        llm = model_info["model"]

        # System prompt for ACMG classification
        system_prompt = """You are an expert clinical geneticist performing ACMG/AMP variant classification.
Analyze the variant data provided and classify it according to the ACMG/AMP 2015 guidelines.
Provide the classification (Pathogenic, Likely_pathogenic, Uncertain_significance, Likely_benign, or Benign)
followed by the specific ACMG criteria that support your classification."""

        # Generate using chat format
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"{instruction}\n\n{variant_input}"},
        ]

        response = llm.create_chat_completion(
            messages=messages,
            max_tokens=512,
            temperature=0.1,
            top_p=0.9,
            repeat_penalty=1.1,
        )

        return response["choices"][0]["message"]["content"]

    elif model_info["type"] == "ollama":
        # Ollama uses raw prompt (not fine-tuned)
        prompt = f"{instruction}\n\n{variant_input}\n\nACMG Classification:"
        import requests

        ollama_host = os.environ.get("OLLAMA_HOST", "localhost:11434")
        if not ollama_host.startswith("http"):
            ollama_host = f"http://{ollama_host}"
        response = requests.post(
            f"{ollama_host}/api/generate",
            json={
                "model": model_info["model"],
                "prompt": prompt,
                "stream": False,
            },
            timeout=120,
        )
        return response.json().get("response", "")

    elif model_info["type"] == "mlx":
        from mlx_lm import generate
        from mlx_lm.sample_utils import make_sampler

        model = model_info["model"]
        tokenizer = model_info["tokenizer"]

        # Apply chat template to match training format
        messages = [{"role": "user", "content": f"{instruction}\n\n{variant_input}"}]
        prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
        # Use low temperature for consistent classification
        sampler = make_sampler(temp=0.1)
        return generate(model, tokenizer, prompt=prompt, max_tokens=500, sampler=sampler)

    elif model_info["type"] == "transformers":
        model = model_info["model"]
        tokenizer = model_info["tokenizer"]
        device = model_info["device"]

        # Apply chat template to match training format
        messages = [{"role": "user", "content": f"{instruction}\n\n{variant_input}"}]
        prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)

        inputs = tokenizer(prompt, return_tensors="pt")
        if device == "cuda":
            inputs = {k: v.cuda() for k, v in inputs.items()}

        outputs = model.generate(
            **inputs,
            max_new_tokens=500,
            temperature=0.1,  # Low temperature for consistent output
            do_sample=True,
        )
        # Extract only the generated response (after the prompt)
        full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)
        if "assistant" in full_output.lower():
            return full_output.split("assistant")[-1].strip()
        return full_output

    return "Classification unavailable"


def parse_classification_response(response: str) -> tuple[str, list[dict], float]:
    """Parse model response into structured data."""
    acmg_class = "Uncertain_significance"
    criteria = []
    confidence = 0.5

    response_lower = response.lower()

    if "pathogenic" in response_lower:
        if "likely pathogenic" in response_lower or "likely_pathogenic" in response_lower:
            acmg_class = "Likely_pathogenic"
            confidence = 0.75
        elif "benign" not in response_lower:
            acmg_class = "Pathogenic"
            confidence = 0.9

    if "benign" in response_lower:
        if "likely benign" in response_lower or "likely_benign" in response_lower:
            acmg_class = "Likely_benign"
            confidence = 0.75
        elif "pathogenic" not in response_lower:
            acmg_class = "Benign"
            confidence = 0.9

    import re

    criteria_pattern = r"(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])[:.]?\s*([^.\n]+)"
    for match in re.finditer(criteria_pattern, response, re.IGNORECASE):
        criteria.append({"code": match.group(1).upper(), "description": match.group(2).strip()})

    return acmg_class, criteria, confidence
