#!/usr/bin/env python3
"""Overlay a calibrated axis grid on the fig-13 render (p086.png, 300 dpi)
so symbol positions can be transcribed precisely: green lines every 500 Hz
and 5 dB, red every 1000 Hz / 10 dB, labels in the margins."""
from PIL import Image, ImageDraw

X0, PXHZ = 406.0, (2000.0 - 406.0) / 5000.0
Y0, PXDB = 505.0, (2232.0 - 505.0) / 70.0


def xf(f):
    return X0 + f * PXHZ


def yd(db):
    return Y0 + (90.0 - db) * PXDB


im = Image.open("pa76_png/p086.png").convert("RGB")
dr = ImageDraw.Draw(im)
for f in range(0, 5001, 500):
    col = (255, 60, 60) if f % 1000 == 0 else (0, 170, 0)
    dr.line([(xf(f), yd(90)), (xf(f), yd(20))], fill=col, width=1)
for db in range(20, 91, 5):
    col = (255, 60, 60) if db % 10 == 0 else (0, 170, 0)
    dr.line([(xf(0), yd(db)), (xf(5000), yd(db))], fill=col, width=1)
    dr.text((xf(5000) + 8, yd(db) - 8), str(db), fill=(200, 0, 200))
im.save("pa76_png/p086_grid.png")
print("saved pa76_png/p086_grid.png", im.size)
