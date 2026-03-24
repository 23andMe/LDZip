#!/bin/bash
set -e

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <stats_json_file>"
  exit 1
fi

STATS_FILE="$1"

num_old=$(jq '.num_non_zero_old' "$STATS_FILE")
num_new=$(jq '.num_non_zero_new' "$STATS_FILE")
cor_full=$(jq '.cor_full' "$STATS_FILE")
cor_nonzero=$(jq '.cor_nonzero' "$STATS_FILE")
max_diff=$(jq '.max_diff' "$STATS_FILE")

if [[ "$num_old" -ne "$num_new" ]]; then
  echo "❌ Mismatch in non-zero counts: $num_old vs $num_new"
  exit 1
else
  echo "✅ Non-zero counts match: $num_old = $num_new"
fi

if [[ "$cor_full" != "1" ]]; then
  echo "❌ Full correlation is not 1: $cor_full"
  exit 1
else
  echo "✅ Full correlation is 1"
fi

if [[ "$cor_nonzero" != "1" ]]; then
  echo "❌ Non-zero correlation is not 1: $cor_nonzero"
  exit 1
else
  echo "✅ Non-zero correlation is 1"
fi

if [[ "$max_diff" != "0" ]]; then
  echo "❌ Max diff is not 0: $max_diff"
  exit 1
else
  echo "✅ Max diff is 0"
fi

echo "✅ All checks passed."
