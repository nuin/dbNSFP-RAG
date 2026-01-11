# ACMG Variant Classification API
#
# Build:
#   docker build -t acmg-api .
#
# Run (with Ollama):
#   docker run -p 8000:8000 -e OLLAMA_HOST=host.docker.internal:11434 acmg-api
#
# Run (with bundled model - requires NVIDIA GPU):
#   docker build --build-arg INCLUDE_MODEL=true -t acmg-api-gpu .
#   docker run --gpus all -p 8000:8000 acmg-api-gpu

FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY pyproject.toml .
RUN pip install --no-cache-dir \
    fastapi \
    uvicorn \
    pydantic \
    sentence-transformers \
    faiss-cpu \
    numpy \
    torch --index-url https://download.pytorch.org/whl/cpu \
    requests

# Copy application code
COPY src/ src/
COPY api/ api/

# Copy database (required for variant lookup)
COPY data/vectordb/grch37-ngsgenes/ data/vectordb/grch37-ngsgenes/

# Optional: Copy model for GPU deployments
ARG INCLUDE_MODEL=false
RUN if [ "$INCLUDE_MODEL" = "true" ]; then \
    pip install --no-cache-dir transformers accelerate; \
    fi
COPY models/acmg-classifier/model/ models/acmg-classifier/model/

# Environment
ENV PYTHONUNBUFFERED=1
ENV DBNSFP_DEVICE=cpu
ENV ACMG_DB_PATH=/app/data/vectordb/grch37-ngsgenes

# Default to Ollama backend (set ACMG_MODEL_PATH to use local model)
ENV OLLAMA_HOST=host.docker.internal:11434

EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

CMD ["uvicorn", "api.server:app", "--host", "0.0.0.0", "--port", "8000"]
