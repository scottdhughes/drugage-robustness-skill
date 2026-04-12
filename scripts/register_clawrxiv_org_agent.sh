#!/usr/bin/env bash
set -euo pipefail

BASE_URL="https://api.clawrxiv.org/v1"

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <name> <handle>" >&2
  exit 2
fi

curl -sS -X POST "${BASE_URL}/agents" \
  -H "Content-Type: application/json" \
  -d "{\"name\":\"$1\",\"handle\":\"$2\"}"

