# papers/ — reference literature

## The original PDFs

The six source PDFs live in the **claude.ai Project "ODT & LEN"**, not on disk, so they
could not be exported into this zip automatically. Bring the originals across one of
these ways:

- **You already uploaded them**, so you have local copies — copy those 6 files into
  this `papers/` folder, or
- Open the claude.ai Project → each file → download, and drop them here.

Expected filenames:

| File in the Project | Paper |
|---|---|
| `full_text.pdf` | Kerstein, Ashurst, Wunsch & Nilsen (2001), *JFM* 447 — ODT vector formulation |
| `Ashurst_Kerstein___2005___...mixing_layers.pdf` | Ashurst & Kerstein (2005), *Phys. Fluids* 17 — variable-density ODT |
| `Lignell_et_al___2013___Mesh_adaption...pdf` | Lignell et al. (2013), *TCFD* 27 — adaptive mesh |
| `full_text1.pdf` | Stephens & Lignell (2021), *SoftwareX* 13 — the ODT C++ code |
| `Proc Appl Math and Mech  2023  Medina Méndez...pdf` | Medina Méndez, Sharma, Schmidt & Klein (2023), *PAMM* 23 — far-field sound |
| `thereturntoisotropyofhomogeneousturbulence 1.pdf` | Choi & Lumley (2001), *JFM* 436 — return to isotropy |

> `.gitignore` excludes `papers/*.pdf` by default so large PDFs don't bloat the repo.
> Delete that line if you want them versioned.

## What IS in this folder

`PAPERS.md` — a full annotated bibliography with a substantial technical summary of each
paper (key equations, parameters, findings, and relevance to the project). It is written
to be read directly by Claude Code, so even before you copy the PDFs back in, a session
has the substance of the literature to work from.

**Also add** (not currently in the Project): Sharma, Klein & Schmidt (2022), *Phys.
Fluids* **34(8)** — the round-jet ODT velocity-statistics paper that the acoustics work
builds on.
