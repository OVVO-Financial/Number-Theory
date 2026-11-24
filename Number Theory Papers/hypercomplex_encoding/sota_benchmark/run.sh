#!/bin/sh
set -eu
# Run from Number Theory Papers/hypercomplex_encoding
MODE=${1:-full}
CONFIG=${CONFIG:-/tmp/hypercomplex-sota-frozen.json}
OUT=${OUT:-/tmp/hypercomplex-sota-results}
python -m unittest -v sota_benchmark.test_codec sota_benchmark.test_protocol
if [ "$MODE" = smoke ]; then
  python -m sota_benchmark.benchmark calibrate --smoke --config "$CONFIG"
  python -m sota_benchmark.benchmark evaluate --smoke --config "$CONFIG" --output "$OUT" --synthetic --repeats 2
elif [ "$MODE" = full ]; then
  python -m sota_benchmark.benchmark calibrate --config "$CONFIG"
  python -m sota_benchmark.benchmark evaluate --config "$CONFIG" --output "$OUT" --synthetic --repeats 3
else
  echo 'Usage: sh sota_benchmark/run.sh [smoke|full]' >&2
  exit 2
fi
