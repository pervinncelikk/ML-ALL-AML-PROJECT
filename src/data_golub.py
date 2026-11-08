"""Golub ALL/AML verisinin yüklenmesi ve etiket eşleştirmesi.

Önemli: CSV dosyalarındaki hasta sütunları numara sırasında DEĞİLDİR.
Etiketler sütun SIRASINA göre değil, hasta NUMARASINA göre eşlenir ve doğrulanır.
"""
import numpy as np
import pandas as pd

from config import RAW, LABEL_MAP

META_COLS = ("Gene Description", "Gene Accession Number")


def _split(df: pd.DataFrame):
    """Sayısal sütunları ve her birine ait 'call' (P/A/M) sütununu ayır."""
    cols = list(df.columns)
    expr_cols, call_cols = [], []
    for i, c in enumerate(cols):
        if c in META_COLS or c.startswith("call"):
            continue
        expr_cols.append(c)
        call_cols.append(cols[i + 1])          # her hastanın hemen sağındaki sütun
        assert cols[i + 1].startswith("call"), f"{c} sonrası call sütunu yok"
    return expr_cols, call_cols


def load_golub():
    """72 örnek x 7129 prob matrisini, etiketleri, set bilgisini ve call'ları döndür."""
    d = RAW / "golub"
    tr = pd.read_csv(d / "data_set_ALL_AML_train.csv")
    te = pd.read_csv(d / "data_set_ALL_AML_independent.csv")
    lab = pd.read_csv(d / "actual.csv").set_index("patient")["cancer"]

    assert (tr["Gene Accession Number"] == te["Gene Accession Number"]).all()
    genes = tr[list(META_COLS)].rename(
        columns={"Gene Description": "description", "Gene Accession Number": "probe"})

    Xs, Cs, sets, order = [], [], [], {}
    for name, df in (("train", tr), ("test", te)):
        ec, cc = _split(df)
        order[name] = [int(c) for c in ec]
        X = df[ec].astype(float).T
        X.index = [int(c) for c in ec]
        C = df[cc].T
        C.index = X.index
        Xs.append(X); Cs.append(C); sets += [name] * len(ec)

    X = pd.concat(Xs); X.columns = genes["probe"].values
    C = pd.concat(Cs); C.columns = genes["probe"].values
    s = pd.Series(sets, index=X.index, name="set")

    # --- etiket eşleştirme: hasta numarası üzerinden ---
    y_lab = lab.loc[X.index]                       # numaraya göre seç
    y = y_lab.map(LABEL_MAP).astype(int).rename("y")

    # --- doğrulama ---
    assert X.index.is_unique and len(X) == 72
    assert set(X.index) == set(lab.index)
    counts = pd.crosstab(s, y_lab)
    assert counts.loc["train", "ALL"] == 27 and counts.loc["train", "AML"] == 11
    assert counts.loc["test", "ALL"] == 20 and counts.loc["test", "AML"] == 14
    # Sütun sırasına göre (yanlış) eşleme yapılsaydı kaç etiket hatalı olurdu?
    naive = lab.sort_index().values                # 1..72 sırası
    wrong_if_naive = int((naive != y_lab.values).sum())

    info = {
        "train_order": order["train"],
        "test_order": order["test"],
        "train_sorted": order["train"] == sorted(order["train"]),
        "test_sorted": order["test"] == sorted(order["test"]),
        "n_mislabeled_if_positional": wrong_if_naive,
        "class_counts": counts.to_dict(),
    }
    return X, y, s, C, genes.set_index("probe"), info


def preprocess_log(X: pd.DataFrame, low=20, high=16000):
    """Alt-üst eşik kırpma + log2 (EDA amaçlı; modelde katlar içinde uygulanacak)."""
    return np.log2(X.clip(lower=low, upper=high))


def variation_filter(Xc: pd.DataFrame, fold=5, diff=500):
    """Dudoit ve ark. (2002) filtresi: max/min > fold ve max-min > diff (kırpılmış ölçekte)."""
    mx, mn = Xc.max(), Xc.min()
    return Xc.columns[(mx / mn > fold) & (mx - mn > diff)]
