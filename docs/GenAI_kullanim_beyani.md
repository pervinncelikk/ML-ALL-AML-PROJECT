# Ek A: Üretken Yapay Zekâ (GenAI) Kullanım Beyanı

Makine Öğrenmesi ve Anomali Tespiti Dersi, Dönem Projesi Ara Raporu · Pervin Çelik (503020250010) · Eskişehir Osmangazi Üniversitesi, Bilgisayar Mühendisliği (YL)

## 1. Beyan

Bu çalışmanın hazırlanmasında Anthropic'in Claude Opus 5.5 modeli; literatür taraması, keşifsel veri analizi kodlarının geliştirilmesi ve raporun dizgisi aşamalarında yardımcı araç olarak kullanılmıştır. Veri setlerinin seçimi, örneklem tasarımı, ön işleme yaklaşımı ve bulguların yorumuna ilişkin kararlar tarafımca verilmiş; yapay zekâ tarafından üretilen tüm içerik birincil kaynaklar ve kod çıktılarıyla doğrulanmıştır. Çalışmanın içeriğine ilişkin tüm sorumluluk bana aittir.

## 2. Kullanım tablosu

| Araç / sürüm | Amaç                               | Bölüm / dosya        | Öğrencinin katkısı                                                                                                                         | Doğrulama                                                                                 |
|------------------|----------------------------------------|--------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| Claude Opus 5.5  | Literatür taraması ve kaynak doğrulama | Bölüm II; `references.bib` | Arama protokolünü (5 sorgu, 1999–2026) Google Scholar'da yürüttü; kaynak ölçütlerini yönergeye göre belirledi                                  | DOI'ler Crossref kayıtlarıyla eşlendi; literatür tablosu birincil kaynaklarla karşılaştırıldı |
| Claude Opus 5.5  | Keşifsel veri analizi kodu             | Bölüm III–IV; `src/*.py`  | GSE13159'u seçti ve temin etti; örneklem kapsamını ve 0–1 ölçeğin işlenme biçimini belirledi                                                   | Dudoit filtresi literatürdeki 3.571/3.051 gen sayısını yeniden üretti                         |
| Claude Opus 5.5  | Rapor metni ve LaTeX dizgisi           | `ara_rapor_tr.tex`         | Yönerge isterlerine göre bölüm bölüm inceledi; hatalı ve eksik içeriğin düzeltilmesini sağladı                                                 | Rapordaki tüm sayılar kod çıktılarından otomatik aktarıldı                                    |
| Claude Opus 5.5  | Tekrarlanabilirlik defterleri          | `notebooks/*.ipynb`       | Analizleri Google Colab ortamında bağımsız olarak yeniden çalıştırdı; anotasyon sürüm farkını tespit edip prob–gen eşlemesini depoda sabitledi | Golub 149 ve GSE13159 134 rapor sayısı birebir eşleşti                                        |

## 3. Öğrencinin Katkıları

1. **Benzer veri setinin seçimi ve temini.** GSE13159 (MILE) veri seti; ana veri setiyle aynı klinik soruyu ele alması, GEO erişim numarası, çok merkezli yapısı ve alt tip çeşitliliği gözetilerek seçilmiştir. Seri matris ile GPL80 ve GPL570 anotasyon dosyaları GEO'dan temin edilmiş; 886 MB'lık seri matris parçalara bölünerek çalışma ortamına aktarılmıştır.

2. **Örneklem tasarımı.** Analize tüm ALL (8 alt tip, 750 örnek) ve AML (6 alt tip, 542 örnek) örneklerinin ve iki doku tipinin (kemik iliği, periferik kan) birlikte alınmasına karar verilmiştir. Bu tercih, B- ve T-ALL'ı birlikte içeren Golub ALL tanımıyla örtüşmeyi sağlamış; doku tipinin etkisi EDA'da ayrıca analiz edilmiştir.

3. **Veri ölçeğinin incelenmesi.** GSE13159 değerlerinin 0–1 aralığında olmasının nedeni araştırılmış, yöntemin özgün kaynağı (Liu ve ark., 2006, J. Theor. Biol.) temin edilerek değerlerin kuramsal Beta dağılımına göre kantil normalize edildiği doğrulanmıştır. Bu bulguya dayanarak değerler yayımlandığı haliyle kullanılmış, log dönüşümü uygulanmamış ve veri setleri arası karşılaştırmada veri seti içi standartlaştırma tercih edilmiştir.

4. **Lisans ve kullanım koşulları.** Golub veri setinin Kaggle üzerindeki lisansı (CC0) ve GEO verilerinin kullanım koşulları birincil kaynaklardan kontrol edilmiştir.

5. **Literatür arama protokolünün yürütülmesi.** Rapordaki beş arama sorgusu Google Scholar'da 1999–2026 yayın aralığıyla çalıştırılmış ve belgelenmiştir; sonuçların raporda kullanılan temel kaynaklarla (Ambroise ve McLachlan 2002, Varma ve Simon 2006, Johnson ve ark. 2007, Haferlach ve ark. 2010) tutarlı olduğu teyit edilmiştir.

6. **Raporun yönerge isterlerine göre doğrulanması.** Rapor; problem tanımı (ALL/AML alt tip sınıflandırması), etiket–sütun eşleştirmesinin hasta numarası üzerinden doğrulanması, p ≫ n ve sınıf dengesizliği isterleri açısından madde madde incelenmiştir. Bu inceleme sonucunda literatür tablosundaki bir çalışmanın (Nagra ve ark., 2024) yöntem ve doğrulama bilgileri birincil kaynakla düzeltilmiş, koddan üretilmeyen bir sayısal ifade kod çıktısına (%47) bağlanmış ve bir şekilde açıklama kutusunun veriyi örtmesi giderilmiştir.

7. **Ayırt edici gen analizinin tamamlanması.** GSE13159 için ısı haritası ve hiyerarşik kümeleme analizinin eksik olduğu tespit edilerek tamamlanmıştır. Bu süreçte p-değerlerinin büyük örneklemde sayısal alt sınıra düşmesi nedeniyle sıralamanın \|t\| istatistiğine göre yapılması gerektiği belirlenmiş; en iyi 30 genle yapılan kümeleme örneklerin %97,3'ünü kendi sınıfıyla gruplamıştır.

8. **Tekrarlanabilirliğin bağımsız ortamda doğrulanması.** Analiz betikleri bölüm bölüm Jupyter defterlerine aktarılmış ve Google Colab ortamında sabit kütüphane sürümleriyle yeniden çalıştırılmıştır. Defterlerin ürettiği değerler rapordaki sayılarla otomatik olarak karşılaştırılmış; Golub için 149, GSE13159 için 134 değerin tamamı eşleşmiştir. Rastgelelik içeren UMAP gösteriminde koordinatların ortama göre değiştiği, ancak kümelenme yapısının korunduğu gözlenmiştir. Bağımsız ortamda gözlenen 21 değer farkının GEO anotasyon dosyalarının sürüm farkından kaynaklandığı belirlenmiş; prob–gen eşlemesi depoda sabitlenerek sonuçların ortamdan bağımsız olarak yeniden üretilebilmesi sağlanmıştır.

9. **Sürüm kontrolü ve depo düzeni.** Kod deposu iş adımlarına (veri temini, Golub EDA, GSE13159 EDA ve karşılaştırma, rapor, tekrarlanabilirlik defterleri) göre ayrı commit'lerle oluşturulmuş; ham veri `.gitignore` ile depo dışında tutulmuş ve güncelliğini yitiren çıktılar kaldırılmıştır.

## 4. Önemli istemlerin özeti

Yönerge isterlerinin çıkarılması ve ara rapor kapsamının belirlenmesi; ana veri setine benzer açık veri setinin belirlenmesi; iki veri seti için veri kalitesi, boyut indirgeme, ayırt edici gen ve batch etkisi analizleri; veri setleri arası karşılaştırma; raporun IEEE şablonunda dizilmesi ve kod deposunun tekrarlanabilirlik isterlerine göre hazırlanması.

**Pervin Çelik**
