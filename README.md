<div align="center">

# ODT‑RDT

### Hybrid aeroacoustics: **R**apid **D**istortion **T**heory × **O**ne‑**D**imensional **T**urbulence → Amiet leading‑edge noise

*Replacing the frozen‑isotropic upwash assumption in Amiet's theory with a physically distorted, anisotropic spectrum from a strain‑coupled ODT model — nonlinear gust distortion at a fraction of LES cost.*

<br/>

<!-- build / runtime stack -->
![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![CMake](https://img.shields.io/badge/CMake-%E2%89%A53.15-064F8C?logo=cmake&logoColor=white)
![conda](https://img.shields.io/badge/conda-Miniforge-44A833?logo=anaconda&logoColor=white)
![Cantera](https://img.shields.io/badge/Cantera-2.5.1-E5601F)
![Linux](https://img.shields.io/badge/platform-Linux%20%2F%20Ubuntu-FCC624?logo=linux&logoColor=black)
![License](https://img.shields.io/badge/license-MIT-blue)

<!-- science / topics -->
<br/>

![ODT](https://img.shields.io/badge/One--Dimensional_Turbulence-1f6feb?style=flat-square)
![RDT](https://img.shields.io/badge/Rapid_Distortion_Theory-1f6feb?style=flat-square)
![Amiet](https://img.shields.io/badge/Amiet_Theory-1f6feb?style=flat-square)
![Aeroacoustics](https://img.shields.io/badge/Aeroacoustics-8957e5?style=flat-square)
![LE Noise](https://img.shields.io/badge/Leading--Edge_Noise-8957e5?style=flat-square)
![Anisotropy](https://img.shields.io/badge/Anisotropic_Turbulence-8957e5?style=flat-square)
![SUNDIALS](https://img.shields.io/badge/SUNDIALS_CVODE-238636?style=flat-square)
![Boost](https://img.shields.io/badge/Boost-238636?style=flat-square)
![fmt](https://img.shields.io/badge/fmt-238636?style=flat-square)
![yaml-cpp](https://img.shields.io/badge/yaml--cpp-238636?style=flat-square)

</div>

---

## ⚡ Quick start — one file, one command

You do **not** need to clone anything by hand. Download the single bootstrap
script, run it, and it does the rest: installs dependencies, fetches the code,
and compiles the solver.

```bash
# 1. grab the one script you need
curl -fsSL -O https://raw.githubusercontent.com/Sparsh-Sharma/ODT-RDT/main/build_odt.sh

# 2. run it
bash build_odt.sh
```

<details>
<summary><b>Prefer a true one‑liner?</b></summary>

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/Sparsh-Sharma/ODT-RDT/main/build_odt.sh)
```

> Downloading the file first (the two‑step version above) is recommended, so you
> can re‑run or inspect it later.

</details>

That's it. When it finishes you'll have a compiled `odt.x` ready to run.

```bash
conda activate odt
cd ODT-RDT/run
./odt.x
```

> **Heads‑up:** the solver links its libraries from the `odt` conda environment,
> so always `conda activate odt` before running it.

---

## 🧩 What the script sets up for you

The bootstrap targets a **clean Ubuntu** machine and is **idempotent** — safe to
re‑run; it reuses whatever already exists.

| Stage | What happens |
| --- | --- |
| 🐍 **Conda** | Installs **Miniforge** automatically if `conda` is missing |
| 📦 **Environment** | Creates a self‑contained `odt` env with a matched compiler + linker toolchain |
| 🔬 **Core libs** | **Cantera 2.5.1**, **fmt**, **Boost**, **yaml‑cpp**, **SUNDIALS/CVODE** |
| 🛠️ **Build** | Clean out‑of‑source **CMake** configure, then compiles on **8 cores** |
| 📂 **Install** | Drops the fresh `odt.x` into `run/` and verifies the libraries resolve |
| 📊 **Post‑processing** | Adds a non‑blocking Python stack: `numpy`, `scipy`, `matplotlib`, `pyyaml` |

<details>
<summary><b>Why a pinned, self‑contained environment?</b></summary>

The solver is written against the **Cantera 2.5.x** API, which is not
source‑compatible with Cantera 3.x. Pinning the whole stack inside one conda
environment — compiler, linker, sysroot and libraries all from the same channel —
keeps everything internally consistent and avoids the toolchain/`glibc` mismatches
that arise when mixing a system compiler with conda libraries.

</details>

---

## 🔁 Developing the code

If you're editing the solver, use the fast rebuild loop instead of the full
bootstrap. It **never touches the environment or dependencies** — it just
recompiles what changed and refreshes `run/odt.x`.

```bash
# incremental: recompiles only the files you changed (seconds)
bash rebuild_odt.sh

# full recompile (after editing headers widely or CMakeLists.txt)
CLEAN=1 bash rebuild_odt.sh
```

| Script | Use it when |
| --- | --- |
| `build_odt.sh` | First setup on a new machine, or to rebuild the environment from scratch (`FRESH_ENV=1`) |
| `rebuild_odt.sh` | Your everyday edit → compile → run loop |

---

## 🧠 About the project

`ODT‑RDT` couples **Rapid Distortion Theory** physics with **One‑Dimensional
Turbulence** as a nonlinear distortion engine, feeding a physically distorted,
anisotropic upwash spectrum into **Amiet's acoustic response** for
leading‑edge noise prediction. The goal is to capture **nonlinear gust
distortion near the leading edge** — behaviour the standard frozen‑isotropic
assumption cannot represent — at a small fraction of the cost of a
scale‑resolving LES.

The ODT core is built on the [BYUignite/ODT](https://github.com/BYUignite/ODT)
implementation.

---

## 📋 Requirements

- A **Linux / Ubuntu** machine (the bootstrap installs everything else)
- `git` and `curl` — auto‑installed via `apt` if missing
- Internet access for the initial dependency download

> macOS isn't targeted by the bootstrap script (the Miniforge installer line and
> the optional `apt` calls are Linux‑specific).

---

## 📜 License

Released under the **MIT License**. The underlying ODT solver is also MIT‑licensed
(© BYU Ignite).

---

<div align="center">
<sub>Turning turbulence into a binary, one eddy at a time.</sub>
</div>
