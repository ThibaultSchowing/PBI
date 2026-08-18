#!/bin/bash
# Dynamic entrypoint for the analysis container.
#
# 1. Patches /etc/passwd so that Python's getpass.getuser() works when the
#    container runs as an arbitrary host UID/GID (Red Hat UBI / OpenShift pattern).
#
# 2. Builds the Jupyter Lab command with optional token authentication:
#    - JUPYTER_TOKEN env var not set / empty: no authentication (default)
#    - JUPYTER_TOKEN env var set: token-based auth required
#
# The CMD from Dockerfile is ignored — this entrypoint owns the full command
# so that auth flags can be injected based on environment variables.

# ---------------------------------------------------------------------------
# /etc/passwd patching for arbitrary UID support
# ---------------------------------------------------------------------------
CURRENT_UID=$(id -u)
CURRENT_GID=$(id -g)

export USER="${USER:-jupyter}"
export LOGNAME="${LOGNAME:-${USER}}"

if ! getent passwd "${CURRENT_UID}" > /dev/null 2>&1 && [ -w /etc/passwd ]; then
    echo "jupyter:x:${CURRENT_UID}:${CURRENT_GID}:Jupyter user:${HOME:-/workspace}:/bin/sh" \
        >> /etc/passwd 2>/dev/null || true
fi

# ---------------------------------------------------------------------------
# Jupyter Lab command construction
# ---------------------------------------------------------------------------
CMD_ARGS=(
    "--ip=0.0.0.0"
    "--port=8888"
    "--no-browser"
    "--notebook-dir=/workspace"
    "--ServerApp.disable_check_xsrf=True"
)

if [ -n "${JUPYTER_TOKEN}" ]; then
    CMD_ARGS+=("--ServerApp.token=${JUPYTER_TOKEN}")
    echo "=== Jupyter Lab: token authentication enabled ==="
    echo "=== Open http://localhost:8888/lab?token=${JUPYTER_TOKEN} ==="
else
    CMD_ARGS+=("--ServerApp.token=")
    CMD_ARGS+=("--ServerApp.password=")
    CMD_ARGS+=("--ServerApp.allow_remote_access=True")
    echo "=== Jupyter Lab: running without authentication ==="
fi

# If arguments were passed to `docker compose run`, exec them directly.
# Otherwise start Jupyter Lab (the default service).
if [ $# -gt 0 ]; then
    exec "$@"
else
    exec jupyter lab "${CMD_ARGS[@]}"
fi
