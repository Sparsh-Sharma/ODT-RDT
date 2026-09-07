#!/usr/bin/env python3
"""Render PDF pages to PNG: pdfrender.py file.pdf outdir p1 [p2 ...] [--dpi N]"""
import os
import sys

import fitz

args = [a for a in sys.argv[1:] if not a.startswith("--dpi")]
dpi = 200
for a in sys.argv[1:]:
    if a.startswith("--dpi"):
        dpi = int(a.split("=")[1])
doc = fitz.open(args[0])
os.makedirs(args[1], exist_ok=True)
for p in args[2:]:
    i = int(p)
    px = doc[i - 1].get_pixmap(dpi=dpi)
    out = os.path.join(args[1], f"p{i:03d}.png")
    px.save(out)
    print(out)
