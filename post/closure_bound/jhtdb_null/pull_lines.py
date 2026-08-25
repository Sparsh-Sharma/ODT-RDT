"""Pull an ensemble of velocity lines from JHTDB isotropic1024coarse.

Null-test data for the LOS closure-bound estimator (SCOPING.md O5): lines
along y (the estimator's x_2 / line axis), random (x, z) positions, one
snapshot. Saves to lines.npz: array (n_lines, n_points, 3) with components
ordered (u1, u2, u3) = (u_x, u_y, u_z); u2 is longitudinal on the line.

Testing token cap: 4096 points/request -> batch 4 lines of 1024 points.
"""

import os
import time

import numpy as np
from givernylocal.turbulence_dataset import turb_dataset
from givernylocal.turbulence_toolkit import getData

TOKEN = "edu.jhu.pha.turbulence.testing-201406"
N_GRID = 1024
DX = 2.0 * np.pi / N_GRID
N_LINES = 128
LINES_PER_REQ = 4
SEED = 20260825
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.environ.get("JHTDB_SCRATCH", os.path.join(OUT_DIR, "_jhtdb_tmp"))


def main():
    os.makedirs(SCRATCH, exist_ok=True)
    dataset = turb_dataset(dataset_title="isotropic1024coarse",
                           output_path=SCRATCH, auth_token=TOKEN)

    rng = np.random.default_rng(SEED)
    # random grid-aligned (x, z) positions, distinct, decorrelated
    xz = rng.choice(N_GRID * N_GRID, size=N_LINES, replace=False)
    x_idx, z_idx = xz // N_GRID, xz % N_GRID
    y = np.arange(N_GRID) * DX

    lines = np.empty((N_LINES, N_GRID, 3), dtype=np.float32)
    n_req = N_LINES // LINES_PER_REQ
    for r in range(n_req):
        pts = np.zeros((LINES_PER_REQ * N_GRID, 3), dtype=np.float64)
        for j in range(LINES_PER_REQ):
            i = r * LINES_PER_REQ + j
            sl = slice(j * N_GRID, (j + 1) * N_GRID)
            pts[sl, 0] = x_idx[i] * DX
            pts[sl, 1] = y
            pts[sl, 2] = z_idx[i] * DX
        for attempt in range(4):
            try:
                res = np.array(getData(dataset, "velocity", 0.0, "none",
                                       "none", "field", pts, verbose=False))
                break
            except Exception as exc:
                wait = 20.0 * (attempt + 1)
                print(f"req {r}: attempt {attempt} failed ({exc}); "
                      f"retrying in {wait:.0f}s", flush=True)
                time.sleep(wait)
        else:
            raise RuntimeError(f"request {r} failed after retries")
        vals = res.reshape(-1, 3) if res.ndim == 3 else res
        for j in range(LINES_PER_REQ):
            i = r * LINES_PER_REQ + j
            lines[i] = vals[j * N_GRID:(j + 1) * N_GRID]
        print(f"req {r + 1}/{n_req} done", flush=True)

    out = os.path.join(OUT_DIR, "lines.npz")
    np.savez_compressed(out, lines=lines, x_idx=x_idx, z_idx=z_idx,
                        dx=DX, dataset="isotropic1024coarse", time=0.0,
                        seed=SEED)
    print("saved", out, lines.shape)


if __name__ == "__main__":
    main()
