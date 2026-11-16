#!/usr/bin/env bash
set -euo pipefail
usage() {
  cat <<'HELP'
Ensemble_gene_annotation.sh (placeholder)
Future steps: prepare_input, analysis, summary, all
This file is a scaffold. It does not run analyses or generate results.
HELP
}
case "${1:-help}" in
  help|-h|--help) usage ;;
  all|prepare_input|analysis|summary)
    echo "[NOT IMPLEMENTED] Ensemble_gene_annotation.sh: this step is a reserved template." >&2
    exit 2 ;;
  *) usage; exit 2 ;;
esac
