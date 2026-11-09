#!/usr/bin/env bash
# Tek komutla yeniden üretim: veri indirme -> EDA -> rapor sayıları -> PDF
set -euo pipefail
cd "$(dirname "$0")"
python scripts/download_data.py
cd src
python eda_golub.py
python eda_gse13159.py
python plot_gse_qc.py
python make_macros.py
python build_bib.py
cd ../report
for d in tr; do
  (cd $d && f=$(ls *.tex | head -1) && n=${f%.tex} && pdflatex -interaction=nonstopmode $n >/dev/null && bibtex $n >/dev/null \
   && pdflatex -interaction=nonstopmode $n >/dev/null && pdflatex -interaction=nonstopmode $n >/dev/null && echo "PDF: report/$d/$n.pdf")
done
