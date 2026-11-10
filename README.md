# ALL/AML Lösemi Alt Tip Sınıflandırması

**Ders:** Makine Öğrenmesi ve Anomali Tespiti, Dönem Projesi

**Öğrenci:** Pervin Çelik (503020250010) · **Aşama:** Ara Rapor (EDA, veri setleri, literatür)

**Tanıtım videosu (ara rapor):** https://drive.google.com/file/d/11Pfx1a9HtEmbLNI-eDoWHuMRzq_28QuK/view?usp=sharing

Gen ifadesi (mikrodizi) verisinden **ALL vs. AML alt tip sınıflandırması** için uçtan uca, tekrarlanabilir
ve sızıntısız bir ML çalışması. Ana veri: Golub ve ark. (1999); benzer açık veri: GEO **GSE13159** (MILE).

> Not: Ders modülü ALL'yi "non-cancer" olarak adlandırır; bu yanlıştır. ALL ve AML'nin ikisi de lösemidir.
> Kodlama: ALL = 0, AML = 1.

## Hızlı başlangıç (tek komut)
```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
bash run_all.sh            # veri indir -> EDA -> sayı makroları -> PDF
```
Adım adım:
```bash
python scripts/download_data.py      # data/raw/ altına indirir (ham veri depoya konmaz)
python src/plot_gse_qc.py            # GSE13159 korelasyon tabanlı kalite kontrol şekli
python src/eda_golub.py              # Golub EDA: results/figures/{tr,en}/golub_*, results/tables/golub_*
python src/eda_gse13159.py           # GSE13159 EDA + iki veri seti karşılaştırması
python src/make_macros.py            # rapordaki tüm sayılar -> report/numbers.tex
python src/build_bib.py              # DOI'li kaynakça
```

## Tekrarlanabilirlik
- Sabit tohum: `SEED = 42` (`src/config.py`); kütüphane sürümleri `requirements.txt`.
- Rapordaki **her sayı** `results/tables/*_numbers.json` → `report/numbers.tex` makrolarıyla basılır; elle sayı yazılmaz.
- Etiketler sütun sırasına göre değil, **hasta numarasına göre** eşlenir ve `assert` ile doğrulanır (`src/data_golub.py`).
- Doğrulama kontrolü: Dudoit ve ark. (2002) filtresi literatürdeki 3571 (72 örnek) ve 3051 (38 eğitim) gen sayılarını üretir.

## Depo yapısı
```
scripts/download_data.py   veri indirme (Kaggle API veya ayna; GEO FTP)
src/config.py              sabitler, seed, renkler
src/data_golub.py          Golub yükleme, etiket eşleme, kırpma/log, varyans filtresi
src/data_gse13159.py       GSE13159 series matrix okuma, ALL/AML alt kümesi, prob->gen eşleme
src/eda_golub.py           Golub EDA (kalite, dağılımlar, PCA/t-SNE/UMAP, DE genler, ısı haritası, batch)
src/eda_gse13159.py        GSE13159 EDA ve platformlar arası karşılaştırma
src/make_macros.py         JSON -> LaTeX sayı makroları
report/tr                  IEEE konferans şablonunda ara rapor (TR)
results/figures, tables    tüm şekiller (PDF vektör + 300 dpi PNG) ve tablolar
```

## Veri kaynakları ve lisans
| Veri | Kaynak | Not |
|---|---|---|
| Golub ALL/AML | Golub et al., *Science* 1999, doi:10.1126/science.286.5439.531; Kaggle `crawford/gene-expression` (lisans: CC0 Public Domain) | 72 örnek, 7129 prob, Hu6800 (GPL80) |
| Golub örnek meta verisi | CRAN `hddplot::golubInfo` (GPL) | BM/PB, merkez, B/T-ALL |
| GSE13159 (MILE Stage 1) | GEO GSE13159; Haferlach et al., *J Clin Oncol* 2010, doi:10.1200/JCO.2009.23.4732 | HG-U133 Plus 2.0 (GPL570) |

## GenAI beyanı
Kod ve metin taslaklarında Anthropic Claude (Claude Opus 5.5) yardımcı araç olarak kullanılmıştır; ayrıntılı beyan raporun Ek A'sındadır.
