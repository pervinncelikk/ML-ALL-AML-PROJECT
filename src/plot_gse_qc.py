"""GSE13159 korelasyon tabanlı kalite kontrol şekli (tek panel). Girdi: results/tables/gse_sample_qc.csv"""
import numpy as np
import pandas as pd
from config import TAB, C_ALL, C_AML
from plotstyle import plt, save, T, COL_W

qc = pd.read_csv(TAB / "gse_sample_qc.csv", index_col=0)
mr = qc.mean_corr
q1, q3 = mr.quantile(.25), mr.quantile(.75)
thr = q1 - 1.5 * (q3 - q1)          # ders (Modül 3.2) IQR kuralı
for lang, L in T.items():
    fig, ax = plt.subplots(figsize=(2.3, 1.95))
    bins = np.linspace(mr.min(), mr.max(), 45)
    for c_, col in (("ALL", C_ALL), ("AML", C_AML)):
        ax.hist(mr[qc.label == c_], bins=bins, color=col, alpha=.75, label=c_)
    ax.axvline(thr, color="#0b0b0b", ls="--", lw=.7)
    ax.set_xlabel("mean correlation to other samples" if lang == "en" else "diğer örneklerle ort. korelasyon")
    ax.set_ylabel(L["n"]); ax.legend(frameon=False, loc="upper left"); ax.grid(axis="y")
    n_out = int((mr < thr).sum())
    ax.set_title(f"(d) QC: {n_out} " + ("outliers" if lang == "en" else "aykırı") + " (IQR)")
    save(fig, "gse_qc_corr", lang)
print("aykırı:", int((mr < thr).sum()))
