#!/usr/bin/env bash
set -euo pipefail
# Usage: bash scripts/generate_carbon.sh /path/to/SLICES C10000.csv
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLICES_DIR="$(cd "${1:?Provide the SLICES checkout path}" && pwd)"
OUTPUT_NAME="${2:-C10000.csv}"
if [[ "$OUTPUT_NAME" == */* || "$OUTPUT_NAME" != *.csv ]]; then
  echo 'Output must be a CSV filename without a directory.' >&2
  exit 1
fi
mkdir -p "$REPO_DIR/results/generated"
OUTPUT_PATH="$REPO_DIR/results/generated/$OUTPUT_NAME"
[[ ! -e "$OUTPUT_PATH" ]] || { echo "Output already exists: $OUTPUT_PATH" >&2; exit 1; }
cd "$SLICES_DIR/MatterGPT/1_train_generate"
for file in model/8-31.pt model/8-31.ini model/8-31_vocab.json; do
  [[ -f "$file" ]] || { echo "Missing model file: $file" >&2; exit 1; }
done
python generate.py \
  --model_weight 8-31.pt \
  --output_csv "$OUTPUT_PATH" \
  --prop_targets '[[40],[50],[60],[70],[80],[90],[100]]' \
  --gen_size 5000 --batch_size 8 \
  --train_dataset "$REPO_DIR/data/train_data_c.csv"
