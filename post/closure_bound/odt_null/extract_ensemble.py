"""caro-side: read nullHIT ODT dumps -> uniform-grid line ensemble -> npz.

Runs on caro with the 'odt' conda env python (numpy available):
    python extract_ensemble.py
Output: /gpfs/caro/scratch/ws/shar_sp-AssamTea/odt_nullhit/nullhit_ensemble.npz
    lines  float32 (n_dumps, n_rlz, n_uniform, 3)   (u, v, w)
    times  float64 (n_dumps,)
Line axis is the ODT domain coordinate y in [-L/2, L/2], L = 1.
"""

import glob
import os
import re

import numpy as np

WS = "/gpfs/caro/scratch/ws/shar_sp-AssamTea/odt_nullhit"
CASE = os.environ.get("ODT_CASE", "nullHIT")
N_UNIFORM = 2048
L = 1.0


def read_dump(path):
    d = np.loadtxt(path)
    pos, u, v, w = d[:, 0], d[:, 2], d[:, 3], d[:, 4]
    t = 0.0
    with open(path) as f:
        for line in f:
            m = re.match(r"#\s*time\s*=\s*(\S+)", line)
            if m:
                t = float(m.group(1))
                break
    return t, pos, np.column_stack([u, v, w])


def main():
    rlz_dirs = sorted(glob.glob(os.path.join(WS, "data", CASE, "data",
                                             "data_*")))
    # exclude only the interactive smoke-test realization (shift 999);
    # note 999 also collides with ensemble seeds >= 999 numerically, but the
    # ensemble writes data_00999 itself, so the smoke dir was overwritten by
    # the job's rm -rf — everything present is ensemble output
    rlz_dirs = [d for d in rlz_dirs]
    n_rlz = len(rlz_dirs)
    dump_ids = sorted(int(os.path.basename(p).split("_")[1].split(".")[0])
                      for p in glob.glob(os.path.join(rlz_dirs[0],
                                                      "dmp_*.dat")))
    yg = (np.arange(N_UNIFORM) + 0.5) / N_UNIFORM * L - L / 2.0

    lines = np.empty((len(dump_ids), n_rlz, N_UNIFORM, 3), dtype=np.float32)
    times = np.zeros(len(dump_ids))
    n_bad = 0
    for r, rdir in enumerate(rlz_dirs):
        for di, d in enumerate(dump_ids):
            path = os.path.join(rdir, "dmp_%05d.dat" % d)
            try:
                t, pos, uvw = read_dump(path)
            except Exception:
                lines[di, r] = np.nan
                n_bad += 1
                continue
            times[di] = t
            for c in range(3):
                lines[di, r, :, c] = np.interp(yg, pos, uvw[:, c])
        if (r + 1) % 32 == 0:
            print("realization %d/%d" % (r + 1, n_rlz), flush=True)

    out = os.path.join(WS, "%s_ensemble.npz" % CASE.lower())
    np.savez_compressed(out, lines=lines, times=times, L=L,
                        n_bad=n_bad, case=CASE)
    print("saved", out, lines.shape, "bad reads:", n_bad)


if __name__ == "__main__":
    main()
