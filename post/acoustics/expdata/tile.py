#!/usr/bin/env python3
"""Crop a region of a PNG (and optionally upscale 2x):
tile.py in.png out.png x0 y0 x1 y1 [scale]"""
import sys

from PIL import Image

im = Image.open(sys.argv[1])
x0, y0, x1, y1 = map(int, sys.argv[3:7])
s = int(sys.argv[7]) if len(sys.argv) > 7 else 1
t = im.crop((x0, y0, x1, y1))
if s > 1:
    t = t.resize((t.width * s, t.height * s), Image.LANCZOS)
t.save(sys.argv[2])
print(sys.argv[2], t.size)
