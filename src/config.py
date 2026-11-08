"""Proje genelinde kullanılan sabitler / project-wide constants."""
from pathlib import Path

SEED = 42                       # tüm rastgele işlemler için sabit seed
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"

# Golub ve ark. (1999) / Dudoit ve ark. (2002) ön işleme eşikleri
CLIP_LOW, CLIP_HIGH = 100, 16000  # alt-üst eşik (Dudoit ve ark. 2002); yönergedeki 20 de raporlanır
FOLD_MIN = 5                      # max/min > 5
DIFF_MIN = 500                    # max - min > 500

LABEL_MAP = {"ALL": 0, "AML": 1}  # ALL = 0, AML = 1 (modüldeki kodlama)

# Renkler (renk-körü dostu, doğrulanmış kategorik palet)
C_ALL = "#2a78d6"   # mavi
C_AML = "#eb6834"   # turuncu
C_3 = "#1baf7a"     # aqua (üçüncü kategori)
GREY = "#8a8985"

for p in (PROC, FIG / "en", FIG / "tr", TAB):
    p.mkdir(parents=True, exist_ok=True)
