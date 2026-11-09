"""references.bib -> references_doi.bib: DOI'leri IEEE stilinde kaynakçada görünür kılar."""
import re
from config import ROOT
s = open(ROOT / "report" / "references.bib", encoding="utf-8").read()
def add_doi(m):
    body = m.group(0)
    d = re.search(r"doi=\{([^}]*)\}", body)
    if not d:
        return body
    doi = d.group(1)
    note = f"doi: \\href{{https://doi.org/{doi}}}{{{doi}}}"
    if "note=" in body or "howpublished=" in body:
        body = re.sub(r"(howpublished=\{)([^}]*)\}", lambda k: f"{k.group(1)}{k.group(2)}, {note}}}", body)
        if "howpublished=" not in body:
            body = re.sub(r"(note=\{)([^}]*)\}", lambda k: f"{k.group(1)}{k.group(2)}; {note}}}", body)
        return body
    b = body.rstrip()
    assert b.endswith("}")
    return b[:-1].rstrip() + f",\n  note={{{note}}}}}"
s = re.sub(r"@\w+\{[^@]*", lambda m: add_doi(m) + "\n", s)
open(ROOT / "report" / "references_doi.bib", "w", encoding="utf-8").write(s)
print("ok")
