# Sabitlenmiş prob → gen sembolü eşlemeleri

| Dosya | Platform | Kaynak (GEO "full table") | Satır |
|---|---|---|---|
| `gpl80_probe2symbol.csv` | Affymetrix Hu6800 (Golub) | `GPL80-30376.txt` | 6.428 |
| `gpl570_probe2symbol.csv` | Affymetrix HG-U133 Plus 2.0 (GSE13159) | `GPL570-55999.txt` | 45.782 |

Tablolar, GEO platform sayfalarından Kasım 2026'da indirilen tam tablolardan `src/data_gse13159.py`
içindeki `read_gpl_annotation` ile üretilmiştir (çok genli problarda ilk sembol; sembolsüz problar çıkarılmıştır).

GEO, platform anotasyonlarını zaman zaman günceller; güncel `.annot` dosyalarıyla bazı probların gen sembolü
değişebilir ve bu, ortak gen sayısı gibi sonuçları küçük ölçüde etkiler. Sonuçların her ortamda raporla birebir
aynı olması için kod, bu klasördeki eşlemeler varsa öncelikle onları kullanır.
