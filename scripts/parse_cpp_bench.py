#!/usr/bin/env python3
"""Turn cpp/qaed_bench output (paper_bench.txt) into LaTeX table rows."""
import re
import sys


def sci(x):
    """5.305 -> $5.31\\times10^{0}$"""
    if x == 0:
        return "$0$"
    import math
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10**e
    return f"$%.2f\\times10^{{%d}}$" % (m, e)


def parse(path, runs):
    """Append one dict per '=====' section found in path. A trailing 'rep=N'
    marks a repeated timing run (best-of-N); the deterministic fields are
    identical across reps and the caller takes the min over 'time'/'q'/'aed'."""
    cur = None
    for line in open(path):
        m = re.match(r"===== type=(\w+) n=(\d+) (\w+)(?: rep=\d+)? =====", line)
        if m:
            cur = {"type": m.group(1), "n": int(m.group(2)), "alg": m.group(3)}
            runs.append(cur)
            continue
        if cur is None:
            continue
        m = re.search(r"(aedq|iqrq|skew_aed|skew_iqr) time: ([\d.]+) s \(QR steps: (\d+)", line)
        if m:
            cur["time"] = float(m.group(2))
            cur["steps"] = int(m.group(3))
        m = re.search(r"time to construct Q: ([\d.]+) s(?:, time for AED: ([\d.]+) s)?", line)
        if m:
            cur["q"] = float(m.group(1))
            cur["aed"] = float(m.group(2)) if m.group(2) else None
        m = re.search(r"unitarity .* = ([\d.e+-]+)", line)
        if m:
            cur["e1"] = float(m.group(1))
        m = re.search(r"residual .* = ([\d.e+-]+)", line)
        if m:
            cur["e2"] = float(m.group(1))
        m = re.search(r"eigvec .* = ([\d.e+-]+)", line)
        if m:
            cur["e3"] = float(m.group(1))


def best(sel):
    """Collapse repeated (type, n, alg) runs to one dict with min timing."""
    out = {}
    for r in sel:
        k = (r["type"], r["n"], r["alg"])
        if k not in out:
            out[k] = dict(r)
        else:
            o = out[k]
            for f in ("time", "q", "aed"):
                if r.get(f) is not None and (o.get(f) is None or r[f] < o[f]):
                    o[f] = r[f]
    return list(out.values())


def main(paths):
    runs = []
    for p in paths:
        parse(p, runs)

    for typ in ["full", "hess", "skew"]:
        sel = best([r for r in runs if r["type"] == typ and "time" in r])
        if not sel:
            continue
        algs = [("aed", "QR+AED"), ("iqr", "QR")]
        if typ == "skew":
            algs = [("skew_aed", "QR+AED"), ("skew_iqr", "QR")]
        print(f"% ---- {typ}rand: performance ----")
        ns = sorted({r["n"] for r in sel})
        for n in ns:
            for alg, label in algs:
                r = next((r for r in sel if r["n"] == n and r["alg"] == alg), None)
                if r is None:
                    continue
                aed = sci(r["aed"]) if r.get("aed") else "N/A"
                print(f"{label} & {n} & {r['steps']} & {sci(r['time'])} & "
                      f"{sci(r['q'])} & {aed} \\\\")
        print(f"% ---- {typ}rand: stability ----")
        for n in ns:
            for alg, label in algs:
                r = next((r for r in sel if r["n"] == n and r["alg"] == alg), None)
                if r is None or "e3" not in r:
                    continue
                print(f"{label} & {n} & {sci(r['e1'])} & {sci(r['e2'])} & "
                      f"{sci(r['e3'])} \\\\")


if __name__ == "__main__":
    main(sys.argv[1:] if len(sys.argv) > 1 else ["paper_bench.txt"])
