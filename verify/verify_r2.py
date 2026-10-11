#!/usr/bin/env python3
"""Recompute the v0.3.0 manuscript numbers from the shipped supplementary workbook.

Reads verify/claims_r2.csv. Computed rows name a workbook cell. The stored value is
rounded half-up to the printed precision and compared with the printed value.
Design, cited and structural rows are not checked. The exit code is non-zero when a
computed row fails or has no source.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import openpyxl

REPO = Path(__file__).resolve().parents[1]
CLAIMS = REPO / "verify" / "claims_r2.csv"
REPORT = REPO / "verify" / "verify_report.csv"
WB = REPO / "results" / "v15" / "Supplementary_Tables_v15.xlsx"
NUM = re.compile(r"(?<![\w.])(?:[−-])?\d[\d,]*(?:\.\d+)?(?: × 10[−-]\d+)?%?")


def half_up(x, p):
    q = Decimal("1").scaleb(-int(p))
    return Decimal(str(x)).quantize(q, rounding=ROUND_HALF_UP)


def parse_token(tok):
    t = tok.replace("−", "-").replace(",", "").strip()
    percent = t.endswith("%")
    t = t.rstrip("%")
    m = re.fullmatch(r"([+-]?\d+(?:\.\d+)?) × 10([+-]?\d+)", t)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2)), percent
    return float(t), percent


def numbers_of(value):
    if isinstance(value, bool) or value is None:
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    out = []
    for tok in NUM.findall(str(value).replace("–", " ").replace("−", "-")):
        try:
            out.append(parse_token(tok)[0])
        except ValueError:
            continue
    return out


def stored_of(sheet, cell, formula):
    raw = sheet[cell].value
    if formula.startswith("text:"):
        k = int(formula.split(":")[1])
        nums = numbers_of(raw)
        if k >= len(nums):
            raise IndexError(f"{cell} has {len(nums)} numbers, asked for {k}")
        return nums[k]
    if raw is None or isinstance(raw, bool) or not isinstance(raw, (int, float)):
        nums = numbers_of(raw)
        if not nums:
            raise ValueError(f"{cell} is not numeric ({raw!r})")
        raw = nums[0]
    if "percent" in formula:
        return float(raw) * 100
    return float(raw)


def main() -> int:
    claims = list(csv.DictReader(CLAIMS.open(encoding="utf-8")))
    wb = openpyxl.load_workbook(WB, read_only=True, data_only=True)
    report = []
    for row in claims:
        kind = row["kind"]
        status, got = "skipped", ""
        if kind == "computed":
            if not row["source_cell"]:
                status = "unmatched"
            else:
                sheet_name, cell = row["source_cell"].split("!")
                try:
                    got_v = stored_of(wb[sheet_name], cell, row["formula"])
                    claim, percent = parse_token(row["value"])
                    p = int(row["precision"])
                    # A percent claim whose formula already scaled the stored value
                    # is compared as printed. A bare cell is compared in stored units.
                    left = half_up(got_v, p)
                    right = half_up(claim, p)
                    status = "pass" if left == right else "fail"
                    got = f"{got_v:.8g}"
                except Exception as exc:  # a missing cell is a failed check, not a crash
                    status, got = "fail", str(exc)
        report.append({
            "id": row["id"], "loc": row["loc"], "kind": kind, "value": row["value"],
            "status": status, "stored": got, "source": row["source_cell"], "column": row["column"],
        })
    with REPORT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(report[0].keys()))
        w.writeheader()
        w.writerows(report)
    by_kind = Counter(r["kind"] for r in report)
    comp = [r for r in report if r["kind"] == "computed"]
    by_status = Counter(r["status"] for r in comp)
    print("kinds", dict(by_kind))
    print("computed", dict(by_status), "of", len(comp))
    fails = [r for r in comp if r["status"] == "fail"]
    for r in fails[:30]:
        print("FAIL", r["id"], r["loc"], "printed", r["value"], "stored", r["stored"], r["source"])
    print("report", REPORT)
    return 1 if by_status.get("fail") or by_status.get("unmatched") else 0


if __name__ == "__main__":
    sys.exit(main())
