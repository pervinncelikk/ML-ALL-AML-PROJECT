"""IEEE sütun genişliğine uygun ortak çizim stili ve iki dilli etiketler."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import FIG

COL_W = 3.5      # IEEE tek sütun (inç)
FULL_W = 7.16    # IEEE çift sütun (inç)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5,
    "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.linewidth": 0.6, "grid.color": "#e4e3df",
    "grid.linewidth": 0.5, "lines.linewidth": 1.2, "savefig.dpi": 300,
    "pdf.fonttype": 42, "figure.dpi": 150,
})

T = {
    "en": dict(
        n="Number of samples", train="Train (38)", test="Independent test (34)",
        train_s="train", test_s="test", raw="Raw intensity (AU)",
        log="log$_2$ intensity (clipped 100–16,000)", density="Density",
        sample="Samples (ordered by set and class)", pc="PC", expl="Explained variance (%)",
        cum="Cumulative", ncomp="Principal component", fdr="$-\\log_{10}$ FDR",
        lfc="log$_2$ fold change (AML / ALL)", class_="Class", set_="Set",
        median="Per-sample median log$_2$", iqr="Per-sample IQR log$_2$",
        lfc_tr="log$_2$FC in train", lfc_te="log$_2$FC in test",
        frac="Fraction of genes", pvals="Train-vs-test KS p-value",
        calls_p="Present (P) calls per sample (%)",
        a_raw="(a) Raw values", b_log="(b) After clipping + log$_2$",
        genes_k="Number of probes", nP="Number of samples with P call",
        dataset="Dataset", common="common genes", filt="filtered probes",
    ),
    "tr": dict(
        n="Örnek sayısı", train="Eğitim (38)", test="Bağımsız test (34)",
        train_s="eğitim", test_s="test", raw="Ham yoğunluk (AU)",
        log="log$_2$ yoğunluk (100–16.000 kırpılmış)", density="Yoğunluk",
        sample="Örnekler (set ve sınıfa göre sıralı)", pc="BB", expl="Açıklanan varyans (%)",
        cum="Kümülatif", ncomp="Temel bileşen", fdr="$-\\log_{10}$ FDR",
        lfc="log$_2$ kat değişimi (AML / ALL)", class_="Sınıf", set_="Set",
        median="Örnek başına medyan log$_2$", iqr="Örnek başına IQR log$_2$",
        lfc_tr="Eğitimde log$_2$FC", lfc_te="Testte log$_2$FC",
        frac="Gen oranı", pvals="Eğitim-test KS p-değeri",
        calls_p="Örnek başına Present (P) çağrısı (%)",
        a_raw="(a) Ham değerler", b_log="(b) Kırpma + log$_2$ sonrası",
        genes_k="Prob sayısı", nP="P çağrısı alan örnek sayısı",
        dataset="Veri seti", common="ortak gen", filt="filtrelenmiş prob",
    ),
}


def save(fig, name, lang):
    out = FIG / lang
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out / f"{name}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
