"""Golub ALL/AML verisi için EDA ve veri analitiği (Ara Rapor, A1).

Çalıştırma:  python src/eda_golub.py
Çıktılar  :  results/figures/{en,tr}/golub_*.pdf|png, results/tables/golub_*.csv,
             results/tables/golub_numbers.json (rapordaki her sayı buradan gelir)

Not (sızıntı): Bu betik yalnızca betimseldir; model eğitmez. Ayırt edici gen
sıralaması YALNIZCA eğitim setinde (38) hesaplanır, bağımsız test seti sadece
bu sıralamanın tutarlılığını görmek için kullanılır. Modelleme aşamasında tüm
ön işleme / öznitelik seçimi adımları CV katlarının içinde yeniden yapılacaktır.
"""
import json
import warnings

import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from config import SEED, RAW, TAB, CLIP_LOW, CLIP_HIGH, FOLD_MIN, DIFF_MIN, C_ALL, C_AML, C_3, GREY
from data_golub import load_golub, preprocess_log, variation_filter
from plotstyle import plt, save, T, COL_W, FULL_W

warnings.filterwarnings("ignore", category=RuntimeWarning)
np.random.seed(SEED)
NUM = {}  # rapora giren sayılar

# ------------------------------------------------------------------ 1. veri
X, y, s, C, genes, info = load_golub()
cls = y.map({0: "ALL", 1: "AML"})
meta = pd.read_csv(RAW / "golub" / "golub_sample_info.csv").set_index("Samples").loc[X.index]
assert (meta.cancer.str.startswith("all") == (cls == "ALL")).all(), "meta etiket uyuşmazlığı"
NUM["n_samples"], NUM["n_probes"] = X.shape
NUM["n_affx_control"] = int(genes.index.str.startswith("AFFX").sum())
NUM["label_check"] = {k: info[k] for k in ("train_sorted", "test_sorted", "n_mislabeled_if_positional")}
NUM["train_col_order_head"] = info["train_order"][:30]
NUM["test_col_order_head"] = info["test_order"][:12]
NUM["class_counts"] = pd.crosstab(s, cls).to_dict()
mt = pd.crosstab([s, cls], [meta.BM_PB]).join(pd.crosstab([s, cls], meta.Source))
mt.to_csv(TAB / "golub_meta_crosstab.csv")
NUM["meta_crosstab"] = {f"{a}|{b}": r.to_dict() for (a, b), r in mt.iterrows()}
NUM["ALL_B_T"] = meta.cancer.value_counts().to_dict()

# ------------------------------------------------------------------ 2. veri kalitesi
q_rows = []
for name in ("train", "test"):
    v = X[s == name].values
    c = C[s == name].values.ravel()
    q_rows.append({
        "set": name, "n": int((s == name).sum()),
        "missing_NaN": int(np.isnan(v).sum()), "INF": int(np.isinf(v).sum()),
        "neg_or_zero_%": round(100 * (v <= 0).mean(), 2),
        "min": v.min(), "max": v.max(), "median": float(np.median(v)),
        f"<{CLIP_LOW}_%": round(100 * (v < CLIP_LOW).mean(), 2),
        "<20_%": round(100 * (v < 20).mean(), 2),
        f">{CLIP_HIGH}_%": round(100 * (v > CLIP_HIGH).mean(), 3),
        "P_%": round(100 * (c == "P").mean(), 2), "A_%": round(100 * (c == "A").mean(), 2),
        "M_%": round(100 * (c == "M").mean(), 2),
    })
quality = pd.DataFrame(q_rows)
quality.to_csv(TAB / "golub_quality.csv", index=False)
NUM["quality"] = quality.to_dict(orient="records")

# Ham veriye doğrudan log alınırsa oluşacak INF/NaN (modüldeki INF'in kaynağı)
with np.errstate(divide="ignore", invalid="ignore"):
    L0 = np.log2(X.values)
NUM["log_raw_neg_inf"] = int(np.isneginf(L0).sum())
NUM["log_raw_nan"] = int(np.isnan(L0).sum())
NUM["log_raw_bad_%"] = round(100 * (~np.isfinite(L0)).mean(), 2)

# Call istatistikleri: hiç P almamış probe sayısı
nP = (C == "P").sum(axis=0)
NUM["probes_never_P"] = int((nP == 0).sum())
NUM["probes_P_in_all"] = int((nP == 72).sum())
NUM["calls_P_per_sample_min_max"] = [round(100 * (C == "P").mean(axis=1).min(), 1),
                                     round(100 * (C == "P").mean(axis=1).max(), 1)]
# Negatif değerlerin call ile ilişkisi
neg = (X.values <= 0)
NUM["neg_values_with_A_call_%"] = round(100 * (C.values[neg] == "A").mean(), 2)

# Uç değerler: log ölçeğinde gen bazında |robust z| > 5
Xl = preprocess_log(X, CLIP_LOW, CLIP_HIGH)
med = Xl.median(); mad = (Xl - med).abs().median() * 1.4826
rz = ((Xl - med) / mad.replace(0, np.nan)).abs()
NUM["robust_z_gt5_%"] = round(100 * float((rz > 5).values.mean()), 3)
# Ders (Modül 3.2) IQR kuralı: gen bazında Q1-1.5*IQR / Q3+1.5*IQR dışındaki değerler
q1, q3 = Xl.quantile(.25), Xl.quantile(.75); iq = q3 - q1
iqr_out = (Xl < q1 - 1.5 * iq) | (Xl > q3 + 1.5 * iq)
NUM["iqr_outlier_%"] = round(100 * float(iqr_out.values.mean()), 2)
NUM["iqr_outlier_genes_any_%"] = round(100 * float(iqr_out.any().mean()), 1)
NUM["iqr_outlier_sample_max"] = {"patient": int(iqr_out.mean(axis=1).idxmax()), "pct": round(100 * float(iqr_out.mean(axis=1).max()), 2)}
smp_out = (rz > 5).mean(axis=1)
NUM["sample_most_outlier"] = {"patient": int(smp_out.idxmax()), "pct": round(100 * smp_out.max(), 2)}

# ------------------------------------------------------------------ 3. filtre
NUM["filter_genes_72"] = int(len(variation_filter(X.clip(CLIP_LOW, CLIP_HIGH), FOLD_MIN, DIFF_MIN)))
NUM["filter_genes_72_clip20"] = int(len(variation_filter(X.clip(20, CLIP_HIGH), FOLD_MIN, DIFF_MIN)))
keep_tr = variation_filter(X[s == "train"].clip(CLIP_LOW, CLIP_HIGH), FOLD_MIN, DIFF_MIN)
NUM["filter_genes_train"] = int(len(keep_tr))
keep72 = variation_filter(X.clip(CLIP_LOW, CLIP_HIGH), FOLD_MIN, DIFF_MIN)
Xf = Xl[keep72]

# ------------------------------------------------------------------ 4. betimsel
keep72 = variation_filter(X.clip(CLIP_LOW, CLIP_HIGH), FOLD_MIN, DIFF_MIN)
desc = pd.DataFrame({
    "set": s, "class": cls,
    "raw_median": X.median(axis=1), "raw_mean": X.mean(axis=1), "raw_sd": X.std(axis=1),
    "log_median": Xl[keep72].median(axis=1),
    "log_iqr": Xl[keep72].quantile(.75, axis=1) - Xl[keep72].quantile(.25, axis=1),
    "P_call_%": 100 * (C == "P").mean(axis=1),
})
desc.to_csv(TAB / "golub_sample_descriptives.csv")
g = desc.groupby(["set", "class"])[["raw_median", "raw_mean", "log_median", "log_iqr", "P_call_%"]].mean().round(2)
g.to_csv(TAB / "golub_descriptives_by_group.csv")
NUM["descriptives_by_group"] = {f"{a}|{b}": r.to_dict() for (a, b), r in g.iterrows()}
NUM["skew_raw_filtered_median"] = round(float(np.median(stats.skew(X[keep72].values, axis=0))), 2)
NUM["skew_log_filtered_median"] = round(float(np.median(stats.skew(Xl[keep72].values, axis=0))), 2)

# ------------------------------------------------------------------ 5. boyut indirgeme
Z = StandardScaler().fit_transform(Xf)
pca = PCA(n_components=10, random_state=SEED).fit(Z)
P = pca.transform(Z)
ev = 100 * pca.explained_variance_ratio_
NUM["pca_explained_%"] = [round(v, 2) for v in ev]
NUM["pca_cum10_%"] = round(float(ev.sum()), 2)
pca_full = PCA(random_state=SEED).fit(Z)
NUM["pca_n_for_80%"] = int(np.searchsorted(np.cumsum(pca_full.explained_variance_ratio_), .80) + 1)
# PC1 ile sınıf/set ilişkisi
NUM["pc1_vs_class_MWU_p"] = float(stats.mannwhitneyu(P[y == 0, 0], P[y == 1, 0]).pvalue)
NUM["pc1_vs_set_MWU_p"] = float(stats.mannwhitneyu(P[s == "train", 0], P[s == "test", 0]).pvalue)
NUM["pc2_vs_class_MWU_p"] = float(stats.mannwhitneyu(P[y == 0, 1], P[y == 1, 1]).pvalue)
NUM["pc2_vs_set_MWU_p"] = float(stats.mannwhitneyu(P[s == "train", 1], P[s == "test", 1]).pvalue)
NUM["pc1_vs_source_KW_p"] = float(stats.kruskal(*[P[(meta.Source == k).values, 0] for k in meta.Source.unique()]).pvalue)
NUM["pc1_vs_BMPB_MWU_p"] = float(stats.mannwhitneyu(P[(meta.BM_PB == "BM").values, 0], P[(meta.BM_PB == "PB").values, 0]).pvalue)
NUM["pc1_median_by_source"] = pd.Series(P[:, 0], index=X.index).groupby(meta.Source).median().round(2).to_dict()
# PCA uzayında silhouette (sınıf ayrışabilirliği, betimsel)
from sklearn.metrics import silhouette_score
NUM["silhouette_pc12_class"] = round(float(silhouette_score(P[:, :2], y)), 3)
NUM["silhouette_pc12_set"] = round(float(silhouette_score(P[:, :2], (s == "test").astype(int))), 3)

tsne = TSNE(n_components=2, perplexity=15, random_state=SEED, init="pca").fit_transform(Z)
import umap
um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=SEED).fit_transform(Z)

# ------------------------------------------------------------------ 6. ayırt edici genler (yalnız eğitim)
Xtr = Xl.loc[s == "train", keep_tr]; ytr = y[s == "train"]
a, b = Xtr[ytr == 0], Xtr[ytr == 1]
t, p = stats.ttest_ind(b, a, equal_var=False)            # Welch t-testi (AML vs ALL)
pu = stats.mannwhitneyu(b, a, axis=0).pvalue              # Mann–Whitney U
lfc = b.mean() - a.mean()
snr = (b.mean() - a.mean()) / (b.std() + a.std())          # Golub sinyal/gürültü
de = pd.DataFrame({"description": genes.loc[keep_tr, "description"],
                   "log2FC_AML_vs_ALL": lfc, "welch_t": t, "p_welch": p,
                   "fdr_welch": multipletests(p, method="fdr_bh")[1],
                   "p_mwu": pu, "fdr_mwu": multipletests(pu, method="fdr_bh")[1],
                   "golub_SNR": snr}, index=keep_tr)
de = de.sort_values("p_welch")
de.to_csv(TAB / "golub_de_train.csv")
NUM["de_n_fdr05_welch"] = int((de.fdr_welch < .05).sum())
NUM["de_n_fdr01_welch"] = int((de.fdr_welch < .01).sum())
NUM["de_n_fdr05_mwu"] = int((de.fdr_mwu < .05).sum())
NUM["de_n_fdr05_and_absLFC1"] = int(((de.fdr_welch < .05) & (de.log2FC_AML_vs_ALL.abs() > 1)).sum())
NUM["de_up_in_AML_fdr05"] = int(((de.fdr_welch < .05) & (de.log2FC_AML_vs_ALL > 0)).sum())
NUM["de_up_in_ALL_fdr05"] = int(((de.fdr_welch < .05) & (de.log2FC_AML_vs_ALL < 0)).sum())
top50 = de.index[:50]
top20 = de.head(20).copy()
top20["gene"] = top20.description.str.split().str[0]
top20[["description", "log2FC_AML_vs_ALL", "welch_t", "fdr_welch", "golub_SNR"]].round(4).to_csv(TAB / "golub_top20_train.csv")
NUM["top10_probes"] = list(de.index[:10])
NUM["top10_desc"] = list(de.description[:10])

# Bağımsız testte aynı genlerin tutarlılığı (sıralama testte YAPILMAZ)
Xte = Xl.loc[s == "test", keep_tr]; yte = y[s == "test"]
lfc_te = Xte[yte == 1].mean() - Xte[yte == 0].mean()
NUM["lfc_corr_train_test_all_genes_pearson"] = round(float(np.corrcoef(lfc, lfc_te[lfc.index])[0, 1]), 3)
NUM["lfc_sign_agree_top50_%"] = round(100 * float((np.sign(lfc[top50]) == np.sign(lfc_te[top50])).mean()), 1)
pt = stats.ttest_ind(Xte.loc[yte == 1, top50], Xte.loc[yte == 0, top50], equal_var=False).pvalue
NUM["top50_p<0.05_in_test"] = int((pt < .05).sum())

# ------------------------------------------------------------------ 7. batch: eğitim vs test
NUM["median_log_train_vs_test_MWU_p"] = float(stats.mannwhitneyu(desc.log_median[s == "train"], desc.log_median[s == "test"]).pvalue)
NUM["iqr_log_train_vs_test_MWU_p"] = float(stats.mannwhitneyu(desc.log_iqr[s == "train"], desc.log_iqr[s == "test"]).pvalue)
NUM["Pcall_train_vs_test_MWU_p"] = float(stats.mannwhitneyu(desc["P_call_%"][s == "train"], desc["P_call_%"][s == "test"]).pvalue)
ks_res = {}
for k, lab_ in ((0, "ALL"), (1, "AML")):
    A_ = Xl.loc[(s == "train") & (y == k), keep72]; B_ = Xl.loc[(s == "test") & (y == k), keep72]
    pk = np.array([stats.ks_2samp(A_[g_], B_[g_]).pvalue for g_ in keep72])
    ks_res[lab_] = pk
    NUM[f"ks_{lab_}_fdr05_n"] = int((multipletests(pk, method="fdr_bh")[1] < .05).sum())
    NUM[f"ks_{lab_}_p05_%"] = round(100 * float((pk < .05).mean()), 1)
# AFFX kontrol problarının set farkı
affx = [c for c in X.columns if c.startswith("AFFX")]
pa = np.array([stats.mannwhitneyu(X.loc[s == "train", c], X.loc[s == "test", c]).pvalue for c in affx])
NUM["affx_p05_n"] = int((pa < .05).sum())

# ------------------------------------------------------------------ 8. kaydet
pd.DataFrame(P[:, :3], index=X.index, columns=["PC1", "PC2", "PC3"]).assign(set=s, cls=cls).to_csv(TAB / "golub_pca_scores.csv")
with open(TAB / "golub_numbers.json", "w") as f:
    json.dump(NUM, f, indent=1, default=float)

# ------------------------------------------------------------------ 9. şekiller (iki dil)
pal = {"ALL": C_ALL, "AML": C_AML}
mk = {"train": "o", "test": "^"}


import re
def gene_label(desc_, probe):
    """Açıklamadan gen sembolü çıkar; yoksa kısaltılmış açıklama."""
    if "HoxA9" in desc_: return "HOXA9"
    tok = desc_.split()
    if tok and re.fullmatch(r"[A-Z0-9-]{2,10}", tok[0]) and tok[0] not in ("GB", "DEF", "HKR", "T-COMPLEX", "PROTEIN"):
        return tok[0]
    m = re.search(r"\(([A-Z0-9-]{2,10})\)", desc_)
    if m: return m.group(1)
    d = re.sub(r"^GB DEF = |mRNA.*$|gene.*$|, complete.*$", "", desc_).strip()
    return (d[:16] + "…") if len(d) > 16 else d


def scatter_emb(ax, E, title, L, xl, yl):
    for st in ("train", "test"):
        for c_ in ("ALL", "AML"):
            m = (s == st).values & (cls == c_).values
            ax.scatter(E[m, 0], E[m, 1], s=16, marker=mk[st], c=pal[c_], edgecolors="white",
                       linewidths=0.5, label=f"{c_} · {L[st + '_s']}")
    ax.set_title(title); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(True)


for lang, L in T.items():
    # F1 sınıf dağılımı
    fig, ax = plt.subplots(figsize=(COL_W, 1.9))
    ct = pd.crosstab(s, cls).loc[["train", "test"]]
    xs = np.arange(2); w = 0.36
    for i, c_ in enumerate(("ALL", "AML")):
        bars = ax.bar(xs + (i - .5) * w, ct[c_], w - .04, color=pal[c_], label=c_)
        for bb, v in zip(bars, ct[c_]):
            pct = 100 * v / ct.sum(axis=1).iloc[list(bars).index(bb)]
            ax.text(bb.get_x() + bb.get_width() / 2, v + .4, f"{v} ({pct:.0f}%)", ha="center", fontsize=6.5)
    ax.set_xticks(xs, [L["train"], L["test"]]); ax.set_ylabel(L["n"]); ax.set_ylim(0, 31)
    ax.legend(frameon=False, ncols=2, loc="upper right"); ax.grid(axis="y")
    save(fig, "golub_class_distribution", lang)

    # F2 dağılımlar: ham vs log
    order = desc.sort_values(["set", "class"], ascending=[False, True]).index
    fig, axs = plt.subplots(1, 2, figsize=(FULL_W, 1.5))
    for ax_, data, ttl, yl in ((axs[0], X, L["a_raw"], L["raw"]), (axs[1], Xl[keep72], L["b_log"] + f" · {len(keep72)} " + L["filt"], L["log"])):
        bp = ax_.boxplot([data.loc[i].values for i in order], showfliers=False, widths=.6,
                         patch_artist=True, medianprops=dict(color="#0b0b0b", lw=.8),
                         whiskerprops=dict(lw=.5), capprops=dict(lw=.5))
        for patch, i in zip(bp["boxes"], order):
            patch.set(facecolor=pal[cls[i]], edgecolor="white", lw=.3)
        ax_.set_xticks([]); ax_.set_xlabel(L["sample"]); ax_.set_ylabel(yl); ax_.set_title(ttl)
        ntr = (s == "train").sum(); ax_.axvline(ntr + .5, color=GREY, ls="--", lw=.6)
        ax_.text(ntr / 2, ax_.get_ylim()[1], L["train_s"], ha="center", va="top", fontsize=6.5, color="#52514e")
        ax_.text(ntr + (72 - ntr) / 2, ax_.get_ylim()[1], L["test_s"], ha="center", va="top", fontsize=6.5, color="#52514e")
    from matplotlib.patches import Patch
    axs[0].legend(handles=[Patch(color=C_ALL, label="ALL"), Patch(color=C_AML, label="AML")],
                  frameon=False, loc="lower left", ncols=2)
    save(fig, "golub_distributions", lang)

    # F2b yoğunluk (ham değer histogramı + log yoğunluk)
    fig, axs = plt.subplots(1, 2, figsize=(FULL_W, 1.8))
    rv = X.values.ravel()
    axs[0].hist(np.clip(rv, -2000, 6000), bins=120, color=GREY)
    axs[0].axvline(0, color="#0b0b0b", lw=.6); axs[0].axvline(CLIP_LOW, color=C_AML, lw=.8, ls="--")
    axs[0].set_xlabel(L["raw"] + " [−2000, 6000]"); axs[0].set_ylabel(L["n"].split()[0] if lang == "en" else "Frekans")
    axs[0].set_title(L["a_raw"])
    for i in X.index:
        sns.kdeplot(Xl.loc[i, keep72], ax=axs[1], color=pal[cls[i]], lw=.5, alpha=.6)
    axs[1].set_xlabel(L["log"]); axs[1].set_ylabel(L["density"]); axs[1].set_title(L["b_log"] + f" · {len(keep72)} " + L["filt"])
    save(fig, "golub_hist_density", lang)

    # F3 PCA
    fig, axs = plt.subplots(1, 2, figsize=(FULL_W, 2.3), gridspec_kw=dict(width_ratios=[1, 1.35]))
    axs[0].bar(np.arange(1, 11), ev, color=C_ALL, width=.7)
    axs[0].plot(np.arange(1, 11), np.cumsum(ev), color="#0b0b0b", marker="o", ms=3, lw=1, label=L["cum"])
    for i_ in range(3):
        axs[0].text(i_ + 1, ev[i_] + 1, f"{ev[i_]:.1f}", ha="center", fontsize=6)
    axs[0].set_xlabel(L["ncomp"]); axs[0].set_ylabel(L["expl"]); axs[0].set_xticks(range(1, 11))
    axs[0].legend(frameon=False); axs[0].grid(axis="y"); axs[0].set_title("(a) Scree")
    scatter_emb(axs[1], P, "(b) PCA", L, f"{L['pc']}1 ({ev[0]:.1f}%)", f"{L['pc']}2 ({ev[1]:.1f}%)")
    axs[1].legend(frameon=False, fontsize=6, loc="best")
    save(fig, "golub_pca", lang)

    # F3b PCA + UMAP (ana metin için birleşik)
    fig, axs = plt.subplots(1, 3, figsize=(FULL_W, 1.95), gridspec_kw=dict(width_ratios=[.8, 1, 1]))
    axs[0].bar(np.arange(1, 11), ev, color=C_ALL, width=.7)
    axs[0].plot(np.arange(1, 11), np.cumsum(ev), color="#0b0b0b", marker="o", ms=2.5, lw=1, label=L["cum"])
    axs[0].set_xlabel(L["ncomp"]); axs[0].set_ylabel(L["expl"]); axs[0].set_xticks(range(1, 11))
    axs[0].legend(frameon=False); axs[0].grid(axis="y"); axs[0].set_title("(a) Scree")
    scatter_emb(axs[1], P, "(b) PCA", L, f"{L['pc']}1 ({ev[0]:.1f}%)", f"{L['pc']}2 ({ev[1]:.1f}%)")
    scatter_emb(axs[2], um, "(c) UMAP", L, "UMAP 1", "UMAP 2")
    axs[2].legend(frameon=False, fontsize=5.5, loc="best")
    save(fig, "golub_pca_umap", lang)

    # F4 t-SNE / UMAP
    fig, axs = plt.subplots(1, 2, figsize=(FULL_W, 2.3))
    scatter_emb(axs[0], tsne, "(a) t-SNE (perplexity=15)", L, "t-SNE 1", "t-SNE 2")
    scatter_emb(axs[1], um, "(b) UMAP (n_neighbors=15)", L, "UMAP 1", "UMAP 2")
    axs[1].legend(frameon=False, fontsize=6, loc="best")
    save(fig, "golub_tsne_umap", lang)

    # F5 volcano
    fig, ax = plt.subplots(figsize=(COL_W, 2.4))
    sig = de.fdr_welch < .05
    ax.scatter(de.log2FC_AML_vs_ALL[~sig], -np.log10(de.fdr_welch[~sig]), s=3, c="#c9c8c3", lw=0)
    up = sig & (de.log2FC_AML_vs_ALL > 0); dn = sig & (de.log2FC_AML_vs_ALL < 0)
    ax.scatter(de.log2FC_AML_vs_ALL[up], -np.log10(de.fdr_welch[up]), s=4, c=C_AML, lw=0, label=f"AML↑ ({up.sum()})")
    ax.scatter(de.log2FC_AML_vs_ALL[dn], -np.log10(de.fdr_welch[dn]), s=4, c=C_ALL, lw=0, label=f"ALL↑ ({dn.sum()})")
    ax.axhline(-np.log10(.05), color=GREY, ls="--", lw=.6)
    from adjustText import adjust_text
    txts = [ax.text(de.log2FC_AML_vs_ALL[pr], -np.log10(de.fdr_welch[pr]), gene_label(de.description[pr], pr), fontsize=5.5)
            for pr in de.index[:8]]
    adjust_text(txts, ax=ax, arrowprops=dict(arrowstyle="-", lw=.4, color=GREY))
    ax.set_xlabel(L["lfc"]); ax.set_ylabel(L["fdr"]); ax.legend(frameon=False, loc="upper left"); ax.grid(True)
    save(fig, "golub_volcano", lang)

    # F6 ısı haritası + hiyerarşik kümeleme (top-50, eğitimde seçilmiş)
    H = Xl[top50]
    Hz = (H - H.mean()) / H.std()
    col_colors = pd.DataFrame({L["class_"]: cls.map(pal), L["set_"]: s.map({"train": "#52514e", "test": "#d9d8d3"})})
    rowlab = [gene_label(de.description[p_], p_) for p_ in top50]
    cg = sns.clustermap(Hz.T.set_axis(rowlab), method="average", metric="correlation", cmap="RdBu_r",
                        vmin=-2.5, vmax=2.5, col_colors=col_colors, figsize=(5.9, 3.55),
                        xticklabels=True, yticklabels=True, dendrogram_ratio=(.08, .1),
                        cbar_kws=dict(), cbar_pos=(0.02, 0.84, 0.012, 0.1))
    cg.ax_cbar.set_title("z", fontsize=6)
    cg.ax_heatmap.tick_params(labelsize=4.6)
    if lang == "en":
        cut = hierarchy.fcluster(cg.dendrogram_col.linkage, 2, "maxclust")
        cl_ = pd.Series(cut, index=Hz.index)
        maj = cl_.groupby(cls).agg(lambda v: v.mode()[0])
        mis = [int(i) for i in Hz.index if cl_[i] != maj[cls[i]]]
        NUM["heatmap_2cluster_misplaced"] = {str(i): f"{cls[i]}/{meta.cancer[i]}/{meta.BM_PB[i]}/{meta.Source[i]}/{s[i]}" for i in mis}
        with open(TAB / "golub_numbers.json", "w") as f:
            json.dump(NUM, f, indent=1, default=float)
    cg.ax_heatmap.set_xlabel("")
    save(cg.fig, "golub_heatmap_top50", lang)

    # F7 batch: eğitim vs test
    fig, axs = plt.subplots(1, 4, figsize=(FULL_W, 1.85), gridspec_kw=dict(wspace=.55))
    for ax_, col, yl in ((axs[0], "log_median", L["median"]), (axs[1], "P_call_%", L["calls_p"])):
        for j, st in enumerate(("train", "test")):
            for k, c_ in enumerate(("ALL", "AML")):
                v = desc[(desc.set == st) & (desc["class"] == c_)][col]
                xpos = j * 2.4 + k
                ax_.boxplot([v], positions=[xpos], widths=.7, showfliers=False, patch_artist=True,
                            boxprops=dict(facecolor=pal[c_], edgecolor="white"), medianprops=dict(color="#0b0b0b"))
                ax_.scatter(np.full(len(v), xpos) + np.random.uniform(-.15, .15, len(v)), v, s=4, c="#0b0b0b", zorder=3, lw=0)
        ax_.set_xticks([.5, 2.9], [L["train_s"], L["test_s"]]); ax_.set_ylabel(yl); ax_.grid(axis="y")
    axs[0].set_title(f"(a) p(MWU)={NUM['median_log_train_vs_test_MWU_p']:.2f}")
    axs[1].set_title(f"(b) p(MWU)={NUM['Pcall_train_vs_test_MWU_p']:.2f}")
    srcs = ["DFCI", "CALGB", "St-Jude", "CCG"]
    for j, src in enumerate(srcs):
        for c_ in ("ALL", "AML"):
            m = ((meta.Source == src) & (cls == c_)).values
            if m.sum() == 0: continue
            ax_ = axs[2]
            ax_.scatter(np.full(m.sum(), j) + (0.15 if c_ == "AML" else -0.15) + np.random.uniform(-.08, .08, m.sum()),
                        P[m, 0], s=9, c=pal[c_], marker="o", edgecolors="white", lw=.4)
    axs[2].set_xticks(range(4), srcs, rotation=30); axs[2].set_ylabel(f"{L['pc']}1"); axs[2].grid(axis="y")
    axs[2].set_title(f"(c) p(KW)={NUM['pc1_vs_source_KW_p']:.2f}")
    axs[3].scatter(lfc, lfc_te[lfc.index], s=2, c="#c9c8c3", lw=0)
    axs[3].scatter(lfc[top50], lfc_te[top50], s=6, c=C_AML, lw=0, label="top-50")
    lim = [-5, 5]; axs[3].plot(lim, lim, color=GREY, ls="--", lw=.6)
    axs[3].set_xlim(lim); axs[3].set_ylim(lim)
    axs[3].set_xlabel(L["lfc_tr"]); axs[3].set_ylabel(L["lfc_te"])
    axs[3].set_title(f"(d) r={NUM['lfc_corr_train_test_all_genes_pearson']:.2f}")
    axs[3].legend(frameon=False, loc="upper left"); axs[3].grid(True)
    save(fig, "golub_batch", lang)

print(json.dumps({k: NUM[k] for k in ['meta_crosstab','pc1_vs_source_KW_p','pc1_vs_BMPB_MWU_p','pc1_median_by_source','skew_raw_filtered_median','skew_log_filtered_median','median_log_train_vs_test_MWU_p','iqr_log_train_vs_test_MWU_p']}, indent=1, default=float))
