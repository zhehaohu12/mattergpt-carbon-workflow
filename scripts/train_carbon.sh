#!/usr/bin/env bash
set -euo pipefail
# Usage: bash scripts/train_carbon.sh /path/to/SLICES
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLICES_DIR="$(cd "${1:?Provide the SLICES checkout path}" && pwd)"
cd "$SLICES_DIR/MatterGPT/1_train_generate"
if [[ -e model/8-31.pt || -e model/8-31.ini || -e model/8-31_vocab.json ]]; then
  echo 'Existing 8-31 model files found. Use a separate checkout for a new training run.' >&2
  exit 1
fi
python train.py \
  --run_name 8-31 \
  --batch_size 8 \
  --max_epochs 100 \
  --n_embd 128 \
  --n_layer 4 \
  --n_head 4 \
  --learning_rate 3e-4 \
  --train_dataset "$REPO_DIR/data/train_data_c.csv" \
  --val_dataset "$REPO_DIR/data/val_data_c.csv" \
  --slices_column_index 0 \
  --prop_column_index_list 1
