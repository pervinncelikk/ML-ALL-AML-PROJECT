# Ek A — GenAI Kullanım Beyanı (Ara Rapor)

**Öğrenci:** Pervin Çelik (503020250010) · **Ders:** Makine Öğrenmesi (YL), Eskişehir Osmangazi Üniversitesi
**Kod deposu:** https://github.com/pervinncelikk/ML-ALL-AML-PROJECT

## Etik beyan

Bu çalışmanın hazırlanmasında Anthropic'in Claude Opus 5.5 modelinden yararlanılmıştır. Model; literatür
taraması, keşifsel veri analizi kodlarının geliştirilmesi ve raporun biçimlendirilmesi aşamalarında yardımcı
araç olarak kullanılmıştır. Veri setlerinin seçimi ve analiz kapsamına ilişkin kararlar tarafımca verilmiş;
yapay zekâ tarafından üretilen tüm içerik tarafımca gözden geçirilmiş ve doğrulanmıştır. Çalışmanın içeriğine
ilişkin tüm sorumluluk bana aittir.

## Kullanım tablosu

| Araç / sürüm | Amaç | Bölüm / dosya | Öğrencinin yaptığı değişiklik | Doğrulama |
|---|---|---|---|---|
| Claude Opus 5.5 (Kasım 2026) | Literatür taraması | Bölüm II, `report/references.bib` | Kaynak ölçütlerini yönergeye göre belirledi; tabloyu 6 sayfa sınırına göre kısalttı | DOI'ler Crossref kayıtlarıyla kontrol edildi |
| Claude Opus 5.5 | Keşifsel veri analizi | Bölüm III–IV, `src/*.py` | GSE13159'u seçip GEO'dan indirdi; örneklem kapsamını ve 0–1 ölçeğin kullanımını belirledi | Kodlar yeniden çalıştırıldı; filtre sonucu literatürdeki 3571/3051 gen sayısıyla eşleşti |
| Claude Opus 5.5 | Rapor metni ve LaTeX | `report/tr/ara_rapor_tr.tex` | İçeriği yönerge isterleriyle sınırladı; gereksiz kısımları çıkardı | Rapordaki sayılar `results/` çıktılarıyla karşılaştırıldı |
| Claude Opus 5.5 | Kod deposu düzeni | GitHub deposu | Depoyu oluşturdu; dosyaları anlamlı mesajlı commit'lerle yükledi | Depo yapısı ve commit geçmişi kontrol edildi |

## Önemli istemlerin özeti

1. Yönergenin incelenmesi ve ara rapor isterlerinin hazırlanması.
2. Ana veri setine benzer açık veri setinin belirlenmesi ve iki veri setinin keşifsel analizi.
3. Raporun IEEE şablonunda 6 sayfaya uyarlanması ve kod deposunun hazırlanması.
