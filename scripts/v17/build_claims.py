"""Turn the claim skeleton into verify/claims_r2.csv.

A sentence's numbers are matched to one supplementary-table row when that row contains
them, then to a column by the words next to the number. Rows that name a module or a
cell type are used only when the sentence names the same one.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import openpyxl

REPO = Path(__file__).resolve().parents[2]
SKEL = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("claims_skeleton_r2_hints.csv")
WB_PATH = REPO / "results" / "v15" / "Supplementary_Tables_v15.xlsx"
OUT = REPO / "verify" / "claims_r2.csv"

NUM = re.compile(r"(?<![\w.])(?:[−-])?\d[\d,]*(?:\.\d+)?(?: × 10[−-]\d+)?%?")
COLORS = (
    "midnightblue greenyellow lightyellow darkorange darkgreen darkgrey darkred "
    "lightgreen royalblue darkmagenta grey60 skyblue salmon turquoise magenta "
    "yellow brown black green white ivory orange pink cyan blue red tan"
).split()
CELLS = (
    "alveolar macrophage type 2 pneumocyte airway epithelial classical monocyte "
    "non-classical monocyte vascular endothelial interstitial macrophage "
    "ciliated t cell b cell nk cell at2 macrophage monocyte fibroblast "
    "mesothelial club stromal endothelial"
).split("  ")  # placeholder, replaced below
CELLS = [
    "alveolar macrophage", "type 2 pneumocyte", "airway epithelial",
    "classical monocyte", "non-classical monocyte", "ncmonocyte", "cmonocyte",
    "vascular endothelial", "interstitial macrophage", "interstitial fibroblast",
    "ciliated", "macrophage", "monocyte", "fibroblast", "mesothelial", "club",
    "stromal", "endothelial", "t cell", "b cell", "nk", "at2", "at1", "atii",
]


def half_up(x, p):
    q = Decimal("1").scaleb(-int(p))
    return Decimal(str(x)).quantize(q, rounding=ROUND_HALF_UP)


def parse_token(tok):
    t = tok.replace("−", "-").replace(",", "").strip()
    percent = t.endswith("%")
    t = t.rstrip("%")
    m = re.fullmatch(r"([+-]?\d+(?:\.\d+)?) × 10([+-]?\d+)", t)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2)), percent, True
    try:
        return float(t), percent, False
    except ValueError:
        return None, percent, False


def same(stored, claim, precision, percent):
    if stored is None:
        return False
    try:
        stored = float(stored)
    except (TypeError, ValueError):
        return False
    if percent:
        return half_up(stored * 100, precision) == half_up(claim, precision) or half_up(stored, precision) == half_up(claim, precision)
    return half_up(stored, precision) == half_up(claim, precision)


def numbers_of(value):
    out = []
    if isinstance(value, bool) or value is None:
        return out
    if isinstance(value, (int, float)):
        out.append(float(value))
        return out
    if isinstance(value, str):
        for tok in NUM.findall(value.replace("–", " ").replace("−", "-")):
            v, _, _ = parse_token(tok)
            if v is not None:
                out.append(v)
    return out


def entities(text):
    low = text.lower()
    found = []
    for c in COLORS:
        if re.search(rf"(?<![a-z]){c}(?![a-z])", low):
            found.append(c)
    for c in CELLS:
        if c in low:
            found.append(c)
    return found


def role(context, sentence):
    pre = context.split("⟦")[0]
    if re.search(r"±\s*$", pre):
        return "sd"
    if re.search(r"(p\s*[=<>]|p-value)\s*$", pre, re.I):
        return "p"
    if re.search(r"(FDR|q-value)\s*$", pre, re.I):
        return "fdr"
    if re.search(r"(ρ|rho|Spearman|Pearson r|r\s*=)\s*$", pre, re.I):
        return "r"
    if re.search(r"(β|beta|coefficient of)\s*$", pre, re.I):
        return "beta"
    if re.search(r"\bCI\b[^.]{0,12}$", pre):
        return "ci"
    if re.search(r"(men|male)\s*$", pre, re.I):
        return "male"
    if re.search(r"(women|female)\s*$", pre, re.I):
        return "female"
    if re.search(r"(z-score|z score|z\s*=)\s*$", pre, re.I):
        return "z"
    if "range" in pre[-20:].lower() or re.search(r"[–-]\s*$", pre):
        return "range"
    if re.search(r"(donors|samples|mice|profiles|genes|features|tests|hits)\s*$", pre, re.I):
        return "n"
    return ""


ROLE_HEADER = {
    "sd": ("sd", "se", "std"),
    "p": ("p_", "p-value", "pval", "p_value", "exact_p", "p_cr1", "p_wild", "p_mixed"),
    "fdr": ("fdr", "q_"),
    "r": ("rho", "spearman", "pearson", "r_oof", "r_bar", "rbar", "recovery"),
    "beta": ("beta", "coef", "coefficient", "delta"),
    "ci": ("ci_", "lo", "hi", "lower", "upper"),
    "male": ("male", "n_male", "sex"),
    "female": ("female", "n_female"),
    "z": ("z", "coexpression"),
    "range": ("min", "max", "range", "age"),
    "n": ("n", "k", "count", "donors", "samples"),
}


def header_score(header, role_name, context):
    h = header.lower()
    score = 0
    for key in ROLE_HEADER.get(role_name, ()):
        if key in h:
            score += 3
    for word in re.findall(r"[a-z][a-z0-9_]{3,}", context.lower()):
        if word in h:
            score += 1
    return score


def load_index():
    wb = openpyxl.load_workbook(WB_PATH, read_only=True, data_only=True)
    rows = {}
    for ws in wb.worksheets:
        if ws.title == "README":
            continue
        table = list(ws.iter_rows(values_only=True))
        if not table:
            continue
        header = [str(h) if h is not None else "" for h in table[0]]
        for i, row in enumerate(table[1:], start=2):
            nums = []
            label_bits = []
            for j, value in enumerate(row):
                if isinstance(value, str):
                    label_bits.append(value)
                found = numbers_of(value)
                if not found:
                    continue
                col = openpyxl.utils.get_column_letter(j + 1)
                kind = "num" if isinstance(value, (int, float)) and not isinstance(value, bool) else "text"
                for k, num in enumerate(found):
                    nums.append({
                        "sheet": ws.title,
                        "cell": f"{col}{i}",
                        "header": header[j] if j < len(header) else "",
                        "value": num,
                        "kind": kind,
                        "k": 0 if kind == "num" else k,
                        "row": i,
                    })
            if nums:
                rows[(ws.title, i)] = {
                    "nums": nums,
                    "label": " ".join(label_bits)[:240].lower(),
                }
    return rows


def cited(para_tables):
    out = []
    for tok in (para_tables or "").split():
        m = re.fullmatch(r"S(\d+)", tok)
        if m:
            out.append(f"S{int(m.group(1)):02d}_")
    return tuple(out)


def row_ok(label, ents):
    if not ents:
        return True
    return all(e in label for e in ents)


def assign_group(claims, index):
    """Return {id: num-dict} for claims that share one row."""
    if len(claims) < 2:
        return {}
    sentence = claims[0]["sentence"]
    ents = entities(sentence)
    prefixes = cited(claims[0]["para_tables"])
    scored = []
    for key, rec in index.items():
        if ents and not row_ok(rec["label"], ents):
            continue
        used = set()
        covered = []
        for claim in claims:
            v, percent, _ = parse_token(claim["value"])
            if v is None:
                continue
            p = int(claim["precision"])
            for n, num in enumerate(rec["nums"]):
                if n in used:
                    continue
                if same(num["value"], v, p, percent):
                    used.add(n)
                    covered.append(claim["id"])
                    break
        if covered:
            bonus = 1 if prefixes and key[0].startswith(prefixes) else 0
            scored.append((len(covered) + bonus * 0.5, len(covered), bonus, key, covered))
    if not scored:
        return {}
    scored.sort(key=lambda t: (t[0], t[2]), reverse=True)
    best = scored[0]
    second = scored[1][0] if len(scored) > 1 else -1
    if best[1] < 2 or best[0] <= second:
        return {}
    rec = index[best[3]]
    used = set()
    out = {}
    for claim in claims:
        if claim["id"] not in best[4]:
            continue
        v, percent, _ = parse_token(claim["value"])
        p = int(claim["precision"])
        hits = []
        for n, num in enumerate(rec["nums"]):
            if n in used:
                continue
            if same(num["value"], v, p, percent):
                hits.append((n, num))
        if not hits:
            continue
        rname = role(claim["context"], claim["sentence"])
        hits.sort(key=lambda h: header_score(h[1]["header"], rname, claim["context"]), reverse=True)
        n, num = hits[0]
        used.add(n)
        out[claim["id"]] = num
    return out


def assign_one(claim, index):
    v, percent, _ = parse_token(claim["value"])
    if v is None:
        return None
    p = int(claim["precision"])
    ents = entities(claim["sentence"])
    prefixes = cited(claim["para_tables"])
    hits = []
    for key, rec in index.items():
        if ents and not row_ok(rec["label"], ents):
            continue
        for num in rec["nums"]:
            if same(num["value"], v, p, percent):
                score = header_score(num["header"], role(claim["context"], claim["sentence"]), claim["context"])
                if prefixes and key[0].startswith(prefixes):
                    score += 2
                if ents:
                    score += 3
                hits.append((score, num, rec["label"]))
    if not hits:
        return None
    hits.sort(key=lambda h: h[0], reverse=True)
    if len(hits) > 1 and hits[0][0] == hits[1][0]:
        return None
    if hits[0][0] < 3:
        return None
    return hits[0][1]


def formula_of(num, percent):
    if num["kind"] == "text":
        base = f"text:{num['k']}"
    else:
        base = "cell"
    if percent and num["kind"] == "num" and abs(num["value"]) <= 1.5:
        base += ":percent"
    return base


def clause_key(row):
    """Numbers in one comma-separated clause share a source row more often than a whole sentence."""
    sentence = row["sentence"]
    token = row["value"].replace("−", "-")
    at = sentence.find(token)
    if at < 0:
        at = sentence.find(token.replace(",", ""))
    if at < 0:
        return row["group"]
    left = sentence.rfind(", ", 0, at)
    # A thousands separator is a comma followed by three digits, not a clause break.
    while left >= 0 and re.match(r", \d{3}\b", sentence[left:]):
        left = sentence.rfind(", ", 0, left)
    right_candidates = [i for i in (sentence.find(", ", at + len(token)), sentence.find("; ", at)) if i != -1]
    right = min(right_candidates) if right_candidates else len(sentence)
    return row["group"] + "@" + sentence[left + 2 if left >= 0 else 0: right][:80]


def build_lookup(index):
    """Map a printed number to the rows that contain it. Very common numbers are left out."""
    buckets = defaultdict(list)
    overflow = set()
    for key, rec in index.items():
        for num in rec["nums"]:
            for p in range(0, 5):
                for scale, flag in ((1, False), (100, True)):
                    try:
                        token = f"{half_up(num['value'] * scale, p)}"
                    except Exception:
                        continue
                    slot = (p, token, flag)
                    if slot in overflow:
                        continue
                    buckets[slot].append((key, num))
                    if len(buckets[slot]) > 40:
                        overflow.add(slot)
                        del buckets[slot]
    return buckets


def hits_for(claim, buckets):
    v, percent, sci = parse_token(claim["value"])
    if v is None or sci:
        return None
    p = int(claim["precision"])
    keys = [(p, f"{half_up(v, p)}", False)]
    if percent:
        keys.append((p + 2, f"{half_up(v / 100, p + 2)}", True))
        keys.append((p, f"{half_up(v, p)}", True))
    found = []
    for key in keys:
        if key not in buckets:
            continue
        found.extend(buckets[key])
    return found


def assign_lookup(claims, buckets, index):
    """Assign a clause when one row is the unique home of its distinctive numbers."""
    per = {}
    for claim in claims:
        found = hits_for(claim, buckets)
        if not found:
            continue
        per[claim["id"]] = found
    if len(per) < 1:
        return {}
    support = defaultdict(set)
    for cid, found in per.items():
        for key, _num in found:
            support[key].add(cid)
    if not support:
        return {}
    ranked = sorted(support.items(), key=lambda kv: len(kv[1]), reverse=True)
    best_key, best_ids = ranked[0]
    best_n = len(best_ids)
    tied = [key for key, ids in support.items() if len(ids) == best_n]
    if best_n >= 2 and len(tied) == 1:
        rec = index[best_key]
    elif best_n == 1 and len(per) == 1 and len(ranked) == 1 and len(next(iter(per.values()))) <= 3:
        rec = index[best_key]
        best_ids = support[best_key]
    else:
        return {}
    out = {}
    used = set()
    for claim in claims:
        if claim["id"] not in best_ids:
            continue
        v, percent, _ = parse_token(claim["value"])
        p = int(claim["precision"])
        options = []
        for n, num in enumerate(rec["nums"]):
            if n in used:
                continue
            if same(num["value"], v, p, percent):
                options.append((n, num))
        if not options:
            continue
        rname = role(claim["context"], claim["sentence"])
        options.sort(key=lambda h: header_score(h[1]["header"], rname, claim["context"]), reverse=True)
        n, num = options[0]
        used.add(n)
        out[claim["id"]] = num
    return out


def main():
    index = load_index()
    buckets = build_lookup(index)
    rows = list(csv.DictReader(open(SKEL, encoding="utf-8")))
    groups = defaultdict(list)
    for r in rows:
        if r["kind_guess"] == "computed":
            groups[r["group"]].append(r)
    chosen = {}
    for claims in groups.values():
        chosen.update(assign_lookup(claims, buckets, index))
    # A second pass uses the comma-separated clause, which separates numbers from different tables.
    groups = defaultdict(list)
    for r in rows:
        if r["kind_guess"] == "computed" and r["id"] not in chosen:
            groups[clause_key(r)].append(r)
    for claims in groups.values():
        chosen.update(assign_lookup(claims, buckets, index))
    by_cell = {}
    for rec in index.values():
        for num in rec["nums"]:
            by_cell[(num["sheet"], num["cell"], num["k"])] = num
    for r in rows:
        if r["kind_guess"] != "computed" or r["id"] in chosen:
            continue
        cands = [c for c in r.get("candidates", "").split(" | ") if c and not c.startswith("(+")]
        if len(cands) != 1:
            continue
        m = re.match(r"([^!]+)!([A-Z]+\d+)", cands[0])
        if not m:
            continue
        # A single workbook cell at the printed precision is used when the clause match found nothing.
        for k in range(0, 6):
            num = by_cell.get((m.group(1), m.group(2), k))
            if num is None:
                continue
            v, percent, _ = parse_token(r["value"])
            if v is not None and same(num["value"], v, int(r["precision"]), percent):
                chosen[r["id"]] = num
                break
    n_comp = sum(r["kind_guess"] == "computed" for r in rows)
    print("assigned", len(chosen), "of", n_comp)
    fields = ["id", "loc", "group", "value", "precision", "kind", "rule", "sentence",
              "source_file", "source_cell", "column", "formula", "context"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            num = chosen.get(r["id"])
            v, percent, _ = parse_token(r["value"])
            w.writerow({
                "id": r["id"],
                "loc": r["loc"],
                "group": r["group"],
                "value": r["value"],
                "precision": r["precision"],
                "kind": r["kind_guess"],
                "rule": r["rule"],
                "sentence": r["sentence"],
                "source_file": "results/v15/Supplementary_Tables_v15.xlsx" if num else "",
                "source_cell": f"{num['sheet']}!{num['cell']}" if num else "",
                "column": num["header"] if num else "",
                "formula": formula_of(num, percent) if num else "",
                "context": r["context"],
            })
    open_rows = [r for r in rows if r["kind_guess"] == "computed" and r["id"] not in chosen]
    print("unassigned", len(open_rows))
    for r in open_rows[:40]:
        print(r["id"], r["loc"], r["value"], r["context"][:100].replace("\n", " "))


if __name__ == "__main__":
    main()
