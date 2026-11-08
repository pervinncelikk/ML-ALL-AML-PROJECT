"""GSE13159 (MILE Stage 1, Haferlach ve ark. 2010) okuma ve ALL vs AML alt kümesi.

Beklenen dosyalar (data/raw/gse13159/):
  GSE13159_series_matrix.txt.gz   (GEO > Series Matrix File(s) > TXT)
  GPL570 anotasyonu  (GPL570.annot.gz  veya  'Download full table' çıktısı)
  GPL80  anotasyonu  (Golub çipi Hu6800; aynı biçimlerden biri)

Çıktılar (data/processed/):
  gse13159_allaml_expr.parquet  : örnek x prob (float32), yalnız ALL+AML
  gse13159_allaml_meta.csv      : GSM, alt tip (MILE sınıfı), ikili etiket
"""
import gzip
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

from config import RAW, PROC, LABEL_MAP

GDIR = RAW / "gse13159"


def _open(p: Path):
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.suffix == ".gz" else open(p, encoding="utf-8", errors="replace")


def find(pattern: str):
    hits = sorted(GDIR.glob(pattern))
    return hits[0] if hits else None


def read_series_matrix(path: Path):
    """Series matrix: '!' ile başlayan meta satırları + ifade tablosu (bellek dostu iki geçiş)."""
    meta, begin, end = {}, None, None
    with _open(path) as f:
        for i, line in enumerate(f):
            if line.startswith("!series_matrix_table_begin"):
                begin = i; continue
            if line.startswith("!series_matrix_table_end"):
                end = i; break
            if begin is None and line.startswith(("!Sample_", "!Series_")):
                key, *vals = line.rstrip("\n").split("\t")
                meta.setdefault(key, []).append([v.strip('"') for v in vals])
    gsm = meta["!Sample_geo_accession"][0]
    sm = pd.DataFrame(index=gsm)
    for key, rows in meta.items():
        if not key.startswith("!Sample_") or len(rows[0]) != len(gsm):
            continue
        for j, r in enumerate(rows):
            name = key.replace("!Sample_", "") + (f"_{j}" if len(rows) > 1 else "")
            sm[name] = r
    nrows = end - begin - 2                      # başlık satırı hariç
    hdr = pd.read_csv(path, sep="\t", skiprows=begin + 1, nrows=0).columns
    dtypes = {h: np.float32 for h in hdr[1:]}
    expr = pd.read_csv(path, sep="\t", skiprows=begin + 1, nrows=nrows, index_col=0,
                       dtype=dtypes, na_values=["null", "NA", ""], low_memory=False)
    expr.index = expr.index.astype(str)
    series = {k: v for k, v in meta.items() if k.startswith("!Series_")}
    return expr, sm, series


def leukemia_class(sm: pd.DataFrame) -> pd.Series:
    """characteristics sütunlarından 'leukemia class' bilgisini bul."""
    for c in sm.columns:
        if c.startswith("characteristics"):
            v = sm[c].astype(str)
            if v.str.contains("leukemia class", case=False).mean() > .5:
                return v.str.split(":", n=1).str[1].str.strip()
    # yedek: başlık/kaynak alanlarında arama
    for c in sm.columns:
        v = sm[c].astype(str)
        if v.str.contains("ALL|AML|CLL|CML", regex=True).mean() > .5:
            return v
    raise ValueError("leukemia class alanı bulunamadı")


def binary_label(cls: str):
    """MILE alt tipini ALL / AML / None (diğer) olarak eşle."""
    c = cls.upper()
    if "AML" in c:
        return "AML"
    if "ALL" in c and "CLL" not in c and "NON" not in c:   # 'Non-leukemia' içinde ALL yok ama emniyet
        return "ALL"
    return None


def read_gpl_annotation(path: Path) -> pd.DataFrame:
    """GPL anotasyonu: .annot.gz veya 'full table' – 'ID' başlık satırından itibaren oku."""
    with _open(path) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if re.match(r'^"?ID"?\t', l))
    end = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
    df = pd.read_csv(io.StringIO("".join(lines[start:end])), sep="\t", dtype=str)
    df.columns = [c.strip('"') for c in df.columns]
    sym_col = next(c for c in df.columns if c.lower() in ("gene symbol", "gene_symbol", "symbol", "identifier"))
    out = df[["ID", sym_col]].rename(columns={sym_col: "symbol"}).dropna()
    out["symbol"] = out["symbol"].str.split(r"\s*///\s*").str[0].str.strip()
    out = out[out.symbol.str.len() > 0]
    return out.set_index("ID")["symbol"]


def collapse_to_genes(X: pd.DataFrame, probe2sym: pd.Series) -> pd.DataFrame:
    """Örnek x prob -> örnek x gen: her gen için ortalama ifadesi en yüksek prob (yaygın 'maxMean' kuralı)."""
    common = X.columns.intersection(probe2sym.index)
    m = X[common].mean()
    df = pd.DataFrame({"probe": common, "sym": probe2sym[common].values, "mean": m.values})
    best = df.sort_values("mean", ascending=False).drop_duplicates("sym")
    G = X[best.probe.values].copy()
    G.columns = best.sym.values
    return G


def build_subset():
    sm_path = find("GSE13159_series_matrix*.txt*")
    assert sm_path is not None, f"{GDIR} içinde GSE13159 series matrix yok"
    expr, sm, series = read_series_matrix(sm_path)
    cls = leukemia_class(sm)
    sm = sm.assign(mile_class=cls.values)
    sm["label"] = sm.mile_class.map(binary_label)
    keep = sm.index[sm.label.notna()]
    meta = sm.loc[keep, ["mile_class", "label"] + [c for c in sm.columns if c.startswith(("source_name", "characteristics", "data_processing", "platform_id"))]]
    meta["y"] = meta.label.map(LABEL_MAP)
    X = expr[keep].T.astype(np.float32)
    X.to_parquet(PROC / "gse13159_allaml_expr.parquet")
    meta.to_csv(PROC / "gse13159_allaml_meta.csv")
    sm[["mile_class", "label"]].to_csv(PROC / "gse13159_all_samples_classes.csv")
    info = {"n_total_samples": int(expr.shape[1]), "n_probes": int(expr.shape[0]),
            "platform": sorted(set(sm.get("platform_id", pd.Series(["?"])))),
            "processing": sorted(set(sm.get("data_processing", pd.Series(["?"]))))[:3],
            "class_counts_all": sm.mile_class.value_counts().to_dict(),
            "n_ALL": int((meta.label == "ALL").sum()), "n_AML": int((meta.label == "AML").sum()),
            "series_title": series.get("!Series_title", [[""]])[0][0]}
    return X, meta, info


if __name__ == "__main__":
    import json
    X, meta, info = build_subset()
    print(json.dumps(info, indent=1, default=str))
    print(X.shape, float(np.nanmin(X.values)), float(np.nanmax(X.values)))
