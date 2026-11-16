#!/usr/bin/env bash
set -euo pipefail
usage() {
  cat <<'HELP'
2.Ensemble_gene_annotation.sh (placeholder)
Future steps: prepare_input, analysis, summary, all
This file is a scaffold. It does not run analyses or generate results.
HELP
}
case "${1:-help}" in
  help|-h|--help) usage ;;
  all|prepare_input|analysis|summary)
    echo "[NOT IMPLEMENTED] 2.Ensemble_gene_annotation.sh: this step is a reserved template." >&2
    exit 2 ;;
  *) usage; exit 2 ;;
esac
