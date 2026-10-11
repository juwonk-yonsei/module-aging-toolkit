"""verify_claims_r2: recompute the numbers reported in the revised manuscript from the result tables.

Inputs
  claims_r2.csv                  one row per number printed in the manuscript and Supplementary Information
                                 (id, loc, value, precision, kind, ...); produced by extract_claims.py
  claims_expr.py                 EXPR: {claim id: expression}, NOTE: {claim id: note}
  Supplementary_Tables.xlsx      the released result tables (one sheet per table)

Every row of kind 'computed' that has an expression is recomputed and compared with the printed
value at the printed precision. Expressions are Python expressions over the helpers below, which
read the workbook only. A row whose expression returns True/False is a bound ("all p < 0.05");
a row whose expression returns a tuple or list passes only if every element prints as the value
("both effects had screen FDR 0.98").

Output: verify_report.csv (one row per claim) and a summary. The exit code is 0 only if no row
fails and every computed row has an expression or is documented as not contained in the released
tables. Rows of other kinds (design, structural, cited) are listed; where a design setting or a
cited value is stored in the tables, it is checked too.

usage: python verify_claims_r2.py claims_r2.csv Supplementary_Tables.xlsx [verify_report.csv]
"""
import csv, math, os, re, sys
from datetime import datetime
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
import numpy as np
import pandas as pd

CLAIMS = sys.argv[1] if len(sys.argv) > 1 else 'claims_r2.csv'
WB = sys.argv[2] if len(sys.argv) > 2 else 'Supplementary_Tables.xlsx'
OUT = sys.argv[3] if len(sys.argv) > 3 else 'verify_report.csv'

_SHEETS = pd.read_excel(WB, sheet_name=None)


# ----------------------------------------------------------------------------- helpers
def T(sheet):
    """The sheet as a DataFrame."""
    return _SHEETS[sheet]

def _filter(sheet, query=None, **eq):
    df = T(sheet)
    for k, val in eq.items():
        col = k.replace('__', '.')
        if isinstance(val, (list, tuple, set)):
            df = df[df[col].isin(list(val))]
        else:
            df = df[df[col] == val]
    if query:
        df = df.query(query, engine='python')
    return df

def rows(sheet, query=None, **eq):
    return _filter(sheet, query, **eq)

def v(sheet, col, query=None, **eq):
    """The single value of a column in the one row selected by the filters."""
    df = _filter(sheet, query, **eq)
    if len(df) != 1:
        raise ValueError(f'{sheet}: {len(df)} rows match {eq} {query or ""}')
    return df[col].iloc[0]

def col(sheet, column, query=None, **eq):
    return pd.to_numeric(_filter(sheet, query, **eq)[column], errors='coerce').dropna()

def med(sheet, column, query=None, **eq):
    return float(col(sheet, column, query, **eq).median())

def mean(sheet, column, query=None, **eq):
    return float(col(sheet, column, query, **eq).mean())

def mn(sheet, column, query=None, **eq):
    return float(col(sheet, column, query, **eq).min())

def mx(sheet, column, query=None, **eq):
    return float(col(sheet, column, query, **eq).max())

def n(sheet, query=None, **eq):
    """Number of rows selected by the filters."""
    return int(len(_filter(sheet, query, **eq)))

def frac(sheet, query_num, query=None, **eq):
    """Fraction of the selected rows that also satisfy query_num."""
    df = _filter(sheet, query, **eq)
    return len(df.query(query_num, engine='python')) / len(df)

def num(text, i=0):
    """The i-th number inside a text cell such as '62.4 ± 4.9 (55–73)'."""
    vals = re.findall(r'-?\d+(?:\.\d+)?', str(text).replace('−', '-').replace('–', ' '))
    return float(vals[i])

def secs(a, b):
    """Seconds from time stamp a to time stamp b."""
    pa = pd.to_datetime(str(a).replace('Z', '+00:00'))
    pb = pd.to_datetime(str(b).replace('Z', '+00:00'))
    return (pb - pa).total_seconds()

def corr(sheet, a, b, method='pearson', query=None, **eq):
    """Correlation of two columns over the selected rows; 'spearman' ranks both columns first."""
    df = _filter(sheet, query, **eq)[[a, b]].apply(pd.to_numeric, errors='coerce').dropna()
    if method == 'spearman':
        df = df.rank()
    return float(df[a].corr(df[b]))

def tq(dof, level=0.95):
    """Two-sided Student t quantile, e.g. tq(51) = 2.008."""
    from scipy import stats
    return float(stats.t.ppf(0.5 + level / 2, dof))

def tci(coef, t, dof, level=0.95):
    """(lower, upper) confidence limits from a coefficient, its t statistic and residual degrees of freedom."""
    half = tq(dof, level) * abs(coef / t)
    return (coef - half, coef + half)

def power_t(d, n1, n2, df_loss=4, alpha=0.05):
    """Power of a two-sided two-group t test for a standardized effect d (noncentrality |d| sqrt(n1 n2 / n))."""
    from scipy import stats
    n = n1 + n2
    df, nc = n - df_loss, abs(d) * math.sqrt(n1 * n2 / n)
    tc = stats.t.ppf(1 - alpha / 2, df)
    return float(stats.nct.sf(tc, df, nc) + stats.nct.cdf(-tc, df, nc))

ENV = {k: globals()[k] for k in ('T', 'rows', 'v', 'col', 'med', 'mean', 'mn', 'mx', 'n', 'frac', 'num', 'secs',
                                 'corr', 'tq', 'tci', 'power_t')}
ENV.update(abs=abs, round=round, min=min, max=max, len=len, float=float, int=int, sum=sum, all=all, any=any,
           sorted=sorted, set=set, list=list, tuple=tuple, zip=zip, str=str, math=math, pd=pd)


# ----------------------------------------------------------------------------- comparison
def printed_number(tok):
    """(value, form, significant digits of the mantissa) of a printed number such as '0.12', '−1.75', '64%',
    '19,334' or '1.6 × 10−15'."""
    t = tok.strip().rstrip(',;').replace('−', '-').replace('+', '').replace(',', '').rstrip('%')
    m = re.fullmatch(r'(-?\d+(?:\.\d+)?) × 10(-?\d+)', t)
    if m:
        mant = m.group(1).lstrip('-')
        sig = len(mant.replace('.', '').lstrip('0')) or 1
        return float(m.group(1)) * 10 ** int(m.group(2)), 'sci', sig
    return float(t), 'plain', None


def matches(value, tok, precision):
    target, form, sig = printed_number(tok)
    if form == 'sci':                           # as many significant digits as printed, e.g. 6 × 10−7, 1.6 × 10−15
        shown = f'{value:.{sig - 1}e}'
        return abs(float(shown) - target) <= 1e-9 * abs(target), shown
    q = Decimal(1).scaleb(-int(precision))
    d = Decimal(repr(float(value)))
    up = d.quantize(q, rounding=ROUND_HALF_UP)
    even = d.quantize(q, rounding=ROUND_HALF_EVEN)
    t = Decimal(repr(target)).quantize(q, rounding=ROUND_HALF_UP)
    ok = (abs(up) == abs(t) or abs(even) == abs(t)) and (target == 0 or (value >= 0) == (target >= 0) or up == 0)
    return ok, str(up)


def evaluate(expr, printed, precision):
    """Evaluate one expression; return (computed text, result)."""
    try:
        g = dict(ENV); g['__builtins__'] = {}          # helpers as globals so comprehensions see them
        val = eval(expr, g)
        if isinstance(val, (bool, np.bool_)):
            return str(bool(val)), 'pass' if bool(val) else 'FAIL'
        if isinstance(val, (list, tuple, pd.Series)):     # several cells that must all print as the value
            res = [matches(float(x), printed, precision) for x in list(val)]
            return '; '.join(sh for _, sh in res), 'pass' if res and all(ok for ok, _ in res) else 'FAIL'
        ok, shown = matches(float(val), printed, precision)
        return shown, 'pass' if ok else 'FAIL'
    except Exception as e:                                # noqa: BLE001
        return f'error: {e}'[:120], 'ERROR'


def main():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # claims_expr.py next to this script
    from claims_expr import EXPR, NOTE
    claims = list(csv.DictReader(open(CLAIMS, encoding='utf-8')))
    out, counts = [], {}
    for r in claims:
        kind = r['kind']
        expr = EXPR.get(r['id'], '')
        rec = dict(id=r['id'], document=r['document'], section=r['section'], printed=r['printed'], kind=kind,
                   expr=expr, computed='', result='', note=NOTE.get(r['id'], ''))
        if expr:                                          # computed rows, and design or cited rows with a stored value
            rec['computed'], rec['result'] = evaluate(expr, r['value'], r['precision'])
        elif kind == 'computed':
            rec['result'] = 'not_in_tables' if rec['note'].startswith('NOT IN TABLES') else 'unmapped'
        else:
            rec['result'] = 'listed'
        counts[(kind, rec['result'])] = counts.get((kind, rec['result']), 0) + 1
        out.append(rec)
    with open(OUT, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    for k in sorted(counts):
        print(f'{k[0]:>10} {k[1]:>13}: {counts[k]}')
    comp = [o for o in out if o['kind'] == 'computed']
    tally = {x: sum(o['result'] == x for o in comp) for x in ('pass', 'FAIL', 'ERROR', 'not_in_tables', 'unmapped')}
    print(f"computed claims: {len(comp)}; recomputed and matched: {tally['pass']}; failed: {tally['FAIL']}; "
          f"errors: {tally['ERROR']}; not in the released tables: {tally['not_in_tables']}; "
          f"without an expression: {tally['unmapped']}")
    bad = [o for o in out if o['result'] in ('FAIL', 'ERROR')]
    for o in bad:
        print(f"  {o['result']}: {o['id']} ({o['section']}) printed {o['printed']} computed {o['computed']} | {o['expr']}")
    sys.exit(1 if bad or tally['unmapped'] else 0)


if __name__ == '__main__':
    main()
