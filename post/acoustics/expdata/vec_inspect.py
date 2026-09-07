#!/usr/bin/env python3
"""Inventory vector drawings on a PDF page: colours, path sizes, bboxes.
Usage: vec_inspect.py file.pdf page"""
import collections
import sys

import fitz

doc = fitz.open(sys.argv[1])
pg = doc[int(sys.argv[2]) - 1]
print("page rect:", pg.rect)
groups = collections.Counter()
for d in pg.get_drawings():
    col = d.get("color")
    n = sum(len(it[1:]) for it in d["items"] if it[0] == "l")
    key = (str(col), d["type"])
    groups[key] += n
for k, v in groups.most_common(25):
    print(k, "line-segs:", v)
# text blocks (axis labels) near bottom-left of page for calibration hints
for b in pg.get_text("blocks")[:40]:
    x0, y0, x1, y1, t = b[:5]
    t = t.strip().replace("\n", " | ")
    if len(t) < 40:
        print(f"text ({x0:.0f},{y0:.0f}): {t!r}")
