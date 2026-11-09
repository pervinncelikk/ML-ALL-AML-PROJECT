"""Rapordaki sayıları JSON çıktılarından LaTeX makrolarına dönüştürür.

Raporda elle yazılmış sayı yerine \\Gxxx / \\Sxxx makroları kullanılır; böylece rapordaki
her sayı kod çıktısıyla birebir eşleşir (Yönerge 6.9).  Çalıştırma: python src/make_macros.py
"""
import json

from config import TAB, ROOT

out = []


def fmt(v, nd=2, tr=False):
    if isinstance(v, bool):
        v = "evet" if (v and tr) else ("hayır" if tr else str(v).lower())
        return v
    if isinstance(v, int):
        s = f"{v:,}"
        return s.replace(",", ".") if tr else s.replace(",", "{,}")
    if isinstance(v, float):
        if v != 0 and abs(v) < 1e-3:
            m, e = f"{v:.1e}".split("e")
            m = m.replace(".", "{,}") if tr else m
            return f"{m}\\times10^{{{int(e)}}}"
        s = f"{v:.{nd}f}"
        return s.replace(".", "{,}") if tr else s
    return str(v)


def add(name, v, nd=2):
    out.append(f"\\newcommand{{\\{name}}}{{{fmt(v, nd)}}}")
    out.append(f"\\newcommand{{\\{name}TR}}{{{fmt(v, nd, tr=True)}}}")


g = json.load(open(TAB / "golub_numbers.json"))
q = {r["set"]: r for r in g["quality"]}
add("GnegTr", q["train"]["neg_or_zero_%"]); add("GnegTe", q["test"]["neg_or_zero_%"])
add("GminTr", int(q["train"]["min"])); add("GmaxTr", int(q["train"]["max"]))
add("GminTe", int(q["test"]["min"])); add("GmaxTe", int(q["test"]["max"]))
add("GltHundTr", q["train"]["<100_%"]); add("GltHundTe", q["test"]["<100_%"])
add("GltTwentyTr", q["train"]["<20_%"]); add("GgtCeilTr", q["train"][">16000_%"], 3); add("GgtCeilTe", q["test"][">16000_%"], 3)
add("GPpct", round((q["train"]["P_%"] * 38 + q["test"]["P_%"] * 34) / 72, 1), 1)
add("GApct", round((q["train"]["A_%"] * 38 + q["test"]["A_%"] * 34) / 72, 1), 1)
add("GMpct", round((q["train"]["M_%"] * 38 + q["test"]["M_%"] * 34) / 72, 1), 1)
add("GlogNegInf", g["log_raw_neg_inf"]); add("GlogNaN", g["log_raw_nan"]); add("GlogBad", g["log_raw_bad_%"])
add("GneverP", g["probes_never_P"]); add("GallP", g["probes_P_in_all"])
add("GPminS", g["calls_P_per_sample_min_max"][0], 1); add("GPmaxS", g["calls_P_per_sample_min_max"][1], 1)
add("GnegA", g["neg_values_with_A_call_%"]); add("GrobustZ", g["robust_z_gt5_%"], 2); add("GiqrOut", g["iqr_outlier_%"]); add("GiqrGenes", g["iqr_outlier_genes_any_%"], 1)
add("GiqrPat", g["iqr_outlier_sample_max"]["patient"]); add("GiqrPatPct", g["iqr_outlier_sample_max"]["pct"])
add("GoutPat", g["sample_most_outlier"]["patient"]); add("GoutPct", g["sample_most_outlier"]["pct"])
add("GmisPos", g["label_check"]["n_mislabeled_if_positional"])
add("GfiltAll", g["filter_genes_72"]); add("GfiltTwenty", g["filter_genes_72_clip20"]); add("GfiltTr", g["filter_genes_train"])
add("GskewRaw", g["skew_raw_filtered_median"]); add("GskewLog", g["skew_log_filtered_median"])
add("GbelowClip", g["raw_below_clip_%"], 0)
ev = g["pca_explained_%"]
add("GpcA", ev[0], 1); add("GpcB", ev[1], 1); add("GpcC", ev[2], 1); add("GpcTen", g["pca_cum10_%"], 1)
add("GpcEighty", g["pca_n_for_80%"])
add("GpcAclass", g["pc1_vs_class_MWU_p"]); add("GpcAset", g["pc1_vs_set_MWU_p"], 3)
add("GpcBclass", g["pc2_vs_class_MWU_p"]); add("GpcBset", g["pc2_vs_set_MWU_p"])
add("GpcAsource", g["pc1_vs_source_KW_p"]); add("GpcAbmpb", g["pc1_vs_BMPB_MWU_p"])
add("GsilClass", g["silhouette_pc12_class"], 3); add("GsilSet", g["silhouette_pc12_set"], 3)
add("GdeFive", g["de_n_fdr05_welch"]); add("GdeOne", g["de_n_fdr01_welch"]); add("GdeMWU", g["de_n_fdr05_mwu"])
add("GdeLfc", g["de_n_fdr05_and_absLFC1"]); add("GdeUpAML", g["de_up_in_AML_fdr05"]); add("GdeUpALL", g["de_up_in_ALL_fdr05"])
add("GlfcR", g["lfc_corr_train_test_all_genes_pearson"]); add("GsignAgree", g["lfc_sign_agree_top50_%"], 0)
add("GtopTest", g["top50_p<0.05_in_test"])
add("GmedSetP", g["median_log_train_vs_test_MWU_p"], 3); add("GiqrSetP", g["iqr_log_train_vs_test_MWU_p"], 3)
add("GPcallSetP", g["Pcall_train_vs_test_MWU_p"], 2)
add("GksALL", g["ks_ALL_fdr05_n"]); add("GksAML", g["ks_AML_fdr05_n"])
add("GksALLp", g["ks_ALL_p05_%"], 1); add("GksAMLp", g["ks_AML_p05_%"], 1); add("Gaffx", g["affx_p05_n"])
mis = g.get("heatmap_2cluster_misplaced", {})
out.append("\\newcommand{\\GmisClust}{" + ", ".join(sorted(mis, key=int)) + "}")
add("GmisN", len(mis))
mc = g["meta_crosstab"]
add("GtePB", mc["test|ALL"]["PB"] + mc["test|AML"]["PB"])
add("GteSJ", mc["test|ALL"]["St-Jude"] + mc["test|AML"]["St-Jude"]); add("GteCCG", mc["test|AML"]["CCG"])
add("GteCALGB", mc["test|AML"]["CALGB"]); add("GteDFCI", mc["test|ALL"]["DFCI"])
add("GallT", g["ALL_B_T"].get("allT", 0)); add("GallB", g["ALL_B_T"].get("allB", 0))

try:
    s = json.load(open(TAB / "gse_numbers.json"))
    add("Sn", s["n"]); add("Sp", s["p"]); add("SnALL", s["class_counts"].get("ALL", 0)); add("SnAML", s["class_counts"].get("AML", 0))
    add("Stotal", s["build"]["n_total_samples"]); add("Snsub", len(s["subtype_counts"]))
    add("SallPct", round(100 * s["class_counts"].get("ALL", 0) / s["n"], 1), 1)
    add("SamlPct", round(100 * s["class_counts"].get("AML", 0) / s["n"], 1), 1)
    add("Snan", s["nan"]); add("Smin", s["min"]); add("Smax", s["max"]); add("Smedian", s["median"])
    add("Sneg", s["neg_or_zero_%"], 3)
    add("SoutMed", s["outlier_samples_median"]); add("SoutIqr", s["outlier_samples_iqr"]); add("SoutCor", s["outlier_samples_corr"])
    add("SoutAny", s["outlier_samples_any"]); add("Sexpr", s["probes_expressed_20pct"]); add("Siqr", s["probes_iqr_gt_0.5"])
    add("Sskew", s["skew_median_raw_or_log"])
    add("SpcA", s["pca_explained_%"][0], 1); add("SpcB", s["pca_explained_%"][1], 1)
    add("SpcAclass", max(s["pc1_vs_class_MWU_p"], 1e-300)); add("SsilClass", s["silhouette_pc12_class"], 3)
    add("SdeFive", s["de_fdr05"]); add("SdeLfc", s["de_fdr05_absLFC1"])
    add("ScommonG", s["common_genes"]); add("SgolubG", s["golub_genes"]); add("SgseG", s["gse_genes"])
    add("SlfcR", s["lfc_corr_pearson"]); add("SlfcRho", s["lfc_corr_spearman"])
    add("SovFifty", s["top50_overlap"]); add("SovHund", s["top100_overlap"]); add("SovHundExp", s["top100_overlap_expected"])
    add("SovHundP", s["top100_overlap_hypergeom_p"]); add("SovTwoH", s["top200_overlap"])
    add("SmedG", s["golub_common_median"]); add("SmedS", s["gse_common_median"])
    add("SiqrG", s["golub_common_iqr"]); add("SiqrS", s["gse_common_iqr"])
    add("SmeanRho", s["gene_mean_spearman"]); add("SsortDiff", s["sorted_dist_max_absdiff"], 4)
    add("SaffxDrop", s["affx_all_nan_dropped"]); add("ScorMin", s["mean_corr_range"][0], 3); add("ScorMed", s["mean_corr_range"][1], 3)
    tb = s["tissue_by_class"]
    add("SbmALL", tb["ALL"]["BM"]); add("SpbALL", tb["ALL"]["PB"]); add("SbmAML", tb["AML"]["BM"]); add("SpbAML", tb["AML"]["PB"])
    add("SbmN", tb["ALL"]["BM"] + tb["AML"]["BM"]); add("SpbN", tb["ALL"]["PB"] + tb["AML"]["PB"])
    add("SpcBclass", s["pc2_vs_class_MWU_p"]); add("SpcAtissue", s["pc1_vs_tissue_MWU_p"], 3); add("SpcBtissue", s["pc2_vs_tissue_MWU_p"])
    add("StisALL", s["tissue_DE_ALL_fdr05_top5k"]); add("StisAML", s["tissue_DE_AML_fdr05_top5k"])
    add("SpcC", s["pca_explained_%"][2], 1)
    add("SbetaDiff", s["beta_1.2_3_maxdiff"], 3); add("SbetaMed", s["beta_1.2_3_median"], 4)
    add("SsilDraw", s["comb_raw_pc1_dataset_silhouette"], 3); add("SsilDz", s["comb_z_pc12_dataset_silhouette"], 3)
    add("SsilCz", s["comb_z_pc12_class_silhouette"], 3)
    add("SrawPCone", s["comb_raw_explained_%"][0], 1)
    add("SheatAgree", s["heatmap_2cluster_agree_%"], 1)
    _ct = s["heatmap_2cluster_crosstab"]; _aml = max(_ct, key=lambda c: _ct[c]["AML"]); _all = [c for c in _ct if c != _aml][0]
    add("SheatMisALL", _ct[_aml]["ALL"]); add("SheatMisAML", _ct[_all]["AML"])
    out.append("\\newcommand{\\SovGenes}{" + ", ".join(s["top100_overlap_genes"][:12]) + "}")
    out.append("\\newcommand{\\StopTen}{" + ", ".join(str(x) for x in s["top10"][:6]) + "}")
except FileNotFoundError:
    print("UYARI: gse_numbers.json yok – GSE makroları üretilmedi")

(ROOT / "report").mkdir(exist_ok=True)
open(ROOT / "report" / "numbers.tex", "w").write("% OTOMATİK ÜRETİLDİ: src/make_macros.py\n" + "\n".join(out) + "\n")
print(len(out), "makro yazıldı")
