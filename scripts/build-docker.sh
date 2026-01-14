#!/bin/bash
#
# Build the bundled Docker image for ACMG Variant Classification API
#
# This script builds a self-contained Docker image that includes:
# - The FastAPI application
# - Llama 3.2 3B Instruct GGUF model for CPU inference
# - GRCh37 vector database for variant lookups
#
# Usage:
#   ./scripts/build-docker.sh [OPTIONS]
#
# Options:
#   --tag TAG       Image tag (default: acmg-api:bundled)
#   --no-cache      Build without cache
#   --push          Push to registry after build
#   --gpu           Build GPU-enabled image (requires NVIDIA base)
#
# Requirements:
#   - Docker installed and running
#   - ~20GB disk space for build
#   - data/vectordb/grch37-ngsgenes/ must exist
#
# The resulting image will be ~8-10GB containing:
#   - Model (~2GB)
#   - Vector database (~2GB compressed)
#   - Python dependencies (~4GB)

set -e

# Default values
TAG="acmg-api:bundled"
DOCKERFILE="Dockerfile.bundled"
NO_CACHE=""
PUSH=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --tag)
            TAG="$2"
            shift 2
            ;;
        --no-cache)
            NO_CACHE="--no-cache"
            shift
            ;;
        --push)
            PUSH=true
            shift
            ;;
        --gpu)
            DOCKERFILE="Dockerfile.gpu"
            echo "GPU build not yet implemented"
            exit 1
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=== ACMG API Docker Build ==="
echo "Tag: $TAG"
echo "Dockerfile: $DOCKERFILE"
echo ""

# Check prerequisites
if [ ! -d "data/vectordb/grch37-ngsgenes" ]; then
    echo "ERROR: Vector database not found at data/vectordb/grch37-ngsgenes/"
    echo "Please build the database first:"
    echo "  uv run python -m src.main build-panel --panel NGSgenes --build grch37"
    exit 1
fi

# Check database size
DB_SIZE=$(du -sh data/vectordb/grch37-ngsgenes/ | cut -f1)
echo "Vector database size: $DB_SIZE"

# Build the image
echo ""
echo "Building Docker image..."
docker build \
    $NO_CACHE \
    -f "$DOCKERFILE" \
    -t "$TAG" \
    .

# Show image size
echo ""
echo "Build complete!"
docker images "$TAG"

# Push if requested
if [ "$PUSH" = true ]; then
    echo ""
    echo "Pushing to registry..."
    docker push "$TAG"
fi

echo ""
echo "To run the container:"
echo "  docker run -p 8000:8000 $TAG"
echo ""
echo "Or use docker-compose:"
echo "  docker compose -f docker-compose.bundled.yml up"
