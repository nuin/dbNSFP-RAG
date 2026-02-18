// PM2 Ecosystem - ACMG Variant Classification API
//
// Prerequisites (RHEL 8):
//   sudo dnf install -y nodejs npm python3.11 python3.11-pip
//   sudo npm install -g pm2
//   curl -LsSf https://astral.sh/uv/install.sh | sh
//   cd /apps/data/src/dbNSFP-RAG && uv venv --python python3.11 && uv pip install -e ".[api]"
//
// Usage:
//   pm2 start ecosystem.config.js          # start all
//   pm2 start ecosystem.config.js --only acmg-api   # API only
//   pm2 start ecosystem.config.js --only ollama      # Ollama only
//   pm2 logs acmg-api                      # tail logs
//   pm2 monit                              # dashboard
//   pm2 save && pm2 startup                # persist across reboots
//
// Database:
//   The SQLite DB must exist at data/sqlite/grch37-all-panels.db (~397 MB).
//   Copy it from your build machine before starting.

module.exports = {
  apps: [
    // ---------------------------------------------------------------
    // Main API server (FastAPI + uvicorn)
    // ---------------------------------------------------------------
    {
      name: "acmg-api",
      script: ".venv/bin/uvicorn",
      args: "api.server:app --host 0.0.0.0 --port 8029 --workers 2",
      cwd: "/apps/data/src/dbNSFP-RAG",
      interpreter: "none",        // uvicorn is already a Python entry point
      env: {
        ACMG_DB_PATH: "/apps/data/src/dbNSFP-RAG/data/sqlite/grch37-all-panels.db",
        PYTHONUNBUFFERED: "1",
      },
      // Restart policy
      max_restarts: 10,
      min_uptime: "10s",
      restart_delay: 3000,
      autorestart: true,

      // Logging
      log_date_format: "YYYY-MM-DD HH:mm:ss Z",
      error_file: "/apps/data/src/dbNSFP-RAG/logs/acmg-api-error.log",
      out_file: "/apps/data/src/dbNSFP-RAG/logs/acmg-api-out.log",
      merge_logs: true,
      max_size: "50M",           // rotate at 50 MB
      retain: 5,                 // keep 5 rotated files

      // Health
      listen_timeout: 10000,
      kill_timeout: 5000,
    },

    // ---------------------------------------------------------------
    // Ollama LLM server (optional — only needed for LLM inference)
    //
    // Install first: curl -fsSL https://ollama.com/install.sh | sh
    // Then pull the model: ollama pull llama3.2:3b
    // ---------------------------------------------------------------
    {
      name: "ollama",
      script: "/usr/local/bin/ollama",
      args: "serve",
      interpreter: "none",
      env: {
        OLLAMA_HOST: "127.0.0.1:11434",
        OLLAMA_NUM_PARALLEL: "2",
        OLLAMA_MAX_LOADED_MODELS: "1",
      },
      // Restart policy
      max_restarts: 5,
      min_uptime: "10s",
      restart_delay: 5000,
      autorestart: true,

      // Logging
      log_date_format: "YYYY-MM-DD HH:mm:ss Z",
      error_file: "/apps/data/src/dbNSFP-RAG/logs/ollama-error.log",
      out_file: "/apps/data/src/dbNSFP-RAG/logs/ollama-out.log",
      merge_logs: true,
      max_size: "50M",
      retain: 3,

      kill_timeout: 10000,
    },
  ],
};
