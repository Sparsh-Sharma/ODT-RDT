#!/usr/bin/env python3
"""Dump text of selected PDF pages: pdftext.py file.pdf 1-6 [grep]"""
import sys

import fitz

doc = fitz.open(sys.argv[1])
rng = sys.argv[2]
pat = sys.argv[3].lower() if len(sys.argv) > 3 else None
a, b = (rng.split("-") + [rng])[:2] if "-" in rng else (rng, rng)
for i in range(int(a) - 1, min(int(b), len(doc))):
    t = doc[i].get_text()
    if pat and pat not in t.lower():
        continue
    print(f"===== page {i+1} =====")
    print(t)
