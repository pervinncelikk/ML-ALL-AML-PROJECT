"""Ham verileri indirir (ham veri depoya KONMAZ; bkz. .gitignore).

Kullanım:  python scripts/download_data.py            # hepsi
           python scripts/download_data.py --golub    # yalnız Golub
           python scripts/download_data.py --gse      # yalnız GSE13159 + GPL anotasyonları

Golub: Kaggle API (kaggle.json gerekir) denenir; olmazsa Kaggle sürümünün birebir aynısı olan
genel GitHub aynasından indirilir (sayısal içerik Kaggle ile karşılaştırılarak doğrulanmıştır).
"""
import argparse
import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

GOLUB_FILES = ["data_set_ALL_AML_train.csv", "data_set_ALL_AML_independent.csv", "actual.csv"]
GOLUB_MIRROR = "https://raw.githubusercontent.com/tkhan11/Gene-Expression-Classification/master/"
GOLUB_INFO = "https://raw.githubusercontent.com/cran/hddplot/master/data/golubInfo.rda"   # örnek meta verisi (BM/PB, merkez)

GEO = {
    "GSE13159_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13159/matrix/GSE13159_series_matrix.txt.gz",
    "GPL570.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "GPL80.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL80/annot/GPL80.annot.gz",
}


def fetch(url, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print("var:", dest.name); return
    print("indiriliyor:", url)
    urllib.request.urlretrieve(url, dest)


def golub():
    d = RAW / "golub"; d.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(["kaggle", "datasets", "download", "-d", "crawford/gene-expression", "-p", str(d), "--unzip"],
                       check=True, capture_output=True)
        print("Kaggle API ile indirildi")
    except Exception:
        print("Kaggle API yok/başarısız -> GitHub aynası")
        for f in GOLUB_FILES:
            fetch(GOLUB_MIRROR + f, d / f)
    info = d / "golub_sample_info.csv"
    if not info.exists():
        import pyreadr
        tmp = d / "golubInfo.rda"; fetch(GOLUB_INFO, tmp)
        m = pyreadr.read_r(str(tmp))["golubInfo"].rename(columns={"BM.PB": "BM_PB", "tissue.mf": "tissue_mf"})
        m["Samples"] = m["Samples"].astype(int); m.to_csv(info, index=False)
    for f in GOLUB_FILES:
        print(f, hashlib.md5((d / f).read_bytes()).hexdigest())


def gse():
    for name, url in GEO.items():
        try:
            fetch(url, RAW / "gse13159" / name)
        except Exception as e:  # anotasyon dosyası yoksa GEO sayfasından 'Download full table' ile alınabilir
            print(f"UYARI: {name} indirilemedi ({e}). GEO sayfasından elle indirip data/raw/gse13159/ içine koyun.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--golub", action="store_true"); ap.add_argument("--gse", action="store_true")
    a = ap.parse_args()
    if a.golub or not (a.golub or a.gse):
        golub()
    if a.gse or not (a.golub or a.gse):
        gse()
