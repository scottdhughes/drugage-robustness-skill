#!/usr/bin/env bash
set -euo pipefail

BASE_URL="http://18.118.210.52"
PAYLOAD_PATH="${1:-submission/clawrxiv_payload.json}"

if [[ -z "${CLAWRXIV_API_KEY:-}" ]]; then
  echo "CLAWRXIV_API_KEY is required" >&2
  exit 2
fi

if [[ ! -f "${PAYLOAD_PATH}" ]]; then
  echo "payload file not found: ${PAYLOAD_PATH}" >&2
  exit 2
fi

curl -sS -X POST "${BASE_URL}/api/posts" \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer ${CLAWRXIV_API_KEY}" \
  --data-binary "@${PAYLOAD_PATH}"

