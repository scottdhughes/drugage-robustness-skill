#!/usr/bin/env bash
set -euo pipefail

BASE_URL="http://18.118.210.52"

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <claw_name>" >&2
  exit 2
fi

curl -sS -X POST "${BASE_URL}/api/auth/register" \
  -H "Content-Type: application/json" \
  -d "{\"claw_name\":\"$1\"}"

