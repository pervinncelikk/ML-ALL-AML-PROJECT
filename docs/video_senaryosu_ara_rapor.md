# Ara Rapor Tanıtım Videosu — Senaryo (hedef 8–9 dk, en fazla 10 dk)

**Kurallar (Yönerge Bölüm 10):** Kendi sesin (yapay ses/avatar yok) · yüzün açılışta ve kod anlatımı sırasında köşede görünür ·
YouTube / Google Drive / OneDrive · bağlantı en az 6 ay erişilebilir (YouTube'da "Liste dışı" önerilir) · URL raporun ilk sayfasına.

**Kayıt önerisi:** OBS Studio (ekran + web kamerası köşede) veya Zoom/Teams'te kendinle toplantı kaydı.
Ekranda sırayla: rapor PDF'i → VS Code/Jupyter → `results/figures/tr/` şekilleri.

| Süre | Bölüm | Ekranda | Anlatılacak ana noktalar |
|---|---|---|---|
| 0:00–0:45 | Açılış (yüz tam ekran) | Kamera | Ad, numara, ders; proje: ALL vs AML alt tip sınıflandırması. "Modül ALL'yi non-cancer diyor, bu yanlış; ikisi de lösemi." |
| 0:45–1:45 | Problem ve motivasyon | Rapor Bölüm I | Klinik önem; mikrodizi ve ortalama fark/negatif değerler; p≫n ve dengesizlik; 3 araştırma sorusu. |
| 1:45–3:15 | Literatür | Tablo I | Arama sorguları + dahil/dışlama; Golub, SVM-RFE, elastik ağ; Ambroise'in seçim yanlılığı; son dönem %100 iddiaları ve neden şüpheyle bakılmalı; boşluk paragrafı. |
| 3:15–4:30 | Veri setleri | Tablo II | Golub (Hu6800, 38/34, BM/PB, merkezler) ve GSE13159 (GPL570, alt tipler); neden seçtim; ortak gen eşlemesi. |
| 4:30–6:45 | **Kod canlı** (yüz köşede) | `src/data_golub.py`, terminal | `load_golub()` içinde etiket eşlemesini göster, `python src/eda_golub.py` çalıştır, `golub_numbers.json`'da "n_mislabeled_if_positional: 14" ve 3571/3051 gen sayısını göster. Neden bu adımlar katlar içinde yapılmalı (sızıntı) — kendi cümlelerinle. |
| 6:45–8:15 | EDA bulguları | Şekil 1–5 | Ham vs log dağılım; PCA'da PC1=set, PC2=sınıf; ısı haritası (12, 25, 67 numaralı örnekler); merkez–sınıf karışıklığı; GSE13159 ile karşılaştırma. |
| 8:15–9:00 | Plan ve kapanış | Bölüm V | 3 yöntem + gerekçe, iç içe CV, bonus hedefleri, takvim, riskler; GenAI beyanı. |

**Kendi sesinle anlatabilmen gereken kritik sorular (hoca sorabilir):**
1. Neden INF'i ortalama ile doldurmak yerine log öncesi kırpma yaptık?
2. Etiketler neden sütun sırasıyla değil hasta numarasıyla eşlenmeli? Kaç örnek yanlış etiketlenirdi?
3. EDA'da gen sıralamasını neden yalnızca eğitim setinde yaptık?
4. PC1'in set ile ilişkili olması modellemede neyi değiştirir?
5. İki platform (Hu6800 vs HG-U133 Plus 2.0) nasıl eşlendi, kaç ortak gen kaldı?
