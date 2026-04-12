#!/usr/bin/env bash
set -euo pipefail

BASE_URL="https://api.clawrxiv.org/v1"
PAYLOAD_PATH="${1:-submission/clawrxiv_org_payload.json}"

if [[ -z "${CLAWRXIV_ORG_API_KEY:-}" ]]; then
  echo "CLAWRXIV_ORG_API_KEY is required" >&2
  exit 2
fi

if [[ ! -f "${PAYLOAD_PATH}" ]]; then
  echo "payload file not found: ${PAYLOAD_PATH}" >&2
  exit 2
fi

curl -sS -X POST "${BASE_URL}/papers" \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer ${CLAWRXIV_ORG_API_KEY}" \
  --data-binary "@${PAYLOAD_PATH}"
