# Migration guide — folder + HPC setup for Claude Code

Step-by-step to stand up the *ODT & LEN* project in the **Claude Code desktop client**
and wire it to a local project folder and to **DLR CARO** and **DLR CARA**.

> **Mental model — important.** Claude Code runs **locally** and operates on a **local
> folder** (your git working copy of the code + notes). It does **not** "connect to an
> HPC" as a first-class feature. You reach CARO/CARA the way you already do — over
> **SSH in Claude Code's terminal** — or by **running Claude Code on the cluster login
> node**. So there are two things to connect: **one local folder** (always), and **the
> clusters via SSH** (two of them). Both are covered below.

---

## Part A — Install Claude Code

1. Install the **Claude Code desktop client** and sign in with your Anthropic /
   claude.ai account (the same one this project lives under).
2. Confirm the CLI works: open its terminal and run `claude --version`.
3. (Recommended) Have **git** installed and configured locally.

---

## Part B — Which folder to connect (the local project folder)

**Connect one folder: the git working copy of your ODT/LEN code + this kit + your
notes.** That is the folder Claude Code should open.

1. Pick a home for it, e.g. `~/research/odt-len/` (or wherever you keep code).
2. Lay it out as in `recommended-structure.md`. Minimum to start:
   ```
   odt-len/
   ├── CLAUDE.md                 ← move the kit's CLAUDE.md here
   ├── HANDOVER.md               ← from this kit
   ├── docs/                     ← from this kit
   ├── papers/                   ← from this kit (+ drop the 6 PDFs in)
   ├── code/odt/                 ← clone BYUignite/ODT (or your BTU fork) here
   ├── cases/                    ← input.yaml case definitions
   ├── notes/                    ← your research notes / drafts
   └── hpc/                      ← SLURM job scripts + ssh config (see Part C)
   ```
3. Get the code in:
   ```bash
   cd ~/research/odt-len/code
   git clone https://github.com/BYUignite/ODT.git odt
   # …or your BTU/DLR fork's remote instead
   ```
4. **Open this folder in Claude Code** — either:
   - launch the client and **"Open folder" / add `~/research/odt-len/`** as the
     workspace, or
   - from a terminal, `cd ~/research/odt-len && claude`.
5. Because `CLAUDE.md` is at the folder root, Claude Code reads the project context
   automatically at the start of every session.

That is the only local folder you need to connect. Big data does **not** come into it
(it stays on the clusters — see `.gitignore`).

---

## Part C — Connecting the HPC (CARO **and** CARA)

You are connecting **both** DLR clusters. Pick the workflow that fits how you like to work;
most people use **C1**.

### C1 — Drive the clusters from Claude Code's terminal over SSH (recommended)

Claude Code has a terminal. You run SSH/`sbatch`/`scp`/`rsync` there exactly as you would
in any terminal; Claude Code can read the output and help you write job scripts, parse
SLURM logs, and edit cases locally before you push them up.

1. **Set up SSH keys** (once per cluster), if not already:
   ```bash
   ssh-keygen -t ed25519 -C "sparsh-dlr"      # if you don't have a key
   ssh-copy-id <user>@<caro-login-host>       # CARO
   ssh-copy-id <user>@<cara-login-host>       # CARA
   ```
2. **Create `~/.ssh/config` entries** so both clusters are one word to reach. Put a copy
   in the repo at `hpc/ssh_config.sample` (no secrets) and the real one in `~/.ssh/config`:
   ```sshconfig
   Host caro
       HostName    <caro-login-hostname>     # from your DLR/GWDG onboarding
       User        <your-username>
       IdentityFile ~/.ssh/id_ed25519
       # ProxyJump <dlr-gateway>             # if CARO is only reachable via a DLR jump host

   Host cara
       HostName    <cara-login-hostname>     # from your DLR/ZIH onboarding
       User        <your-username>
       IdentityFile ~/.ssh/id_ed25519
       # ProxyJump <dlr-gateway>             # if a jump host is required
   ```
   Then connecting is just `ssh caro` or `ssh cara`.
3. **Confirm access and scheduler** on each:
   ```bash
   ssh caro 'sinfo -s; squeue --me; module avail 2>&1 | head'
   ssh cara 'sinfo -s; squeue --me; module avail 2>&1 | head'
   ```
   Both are SLURM clusters, so `sbatch`, `squeue`, `sacct`, `scancel` all apply.
4. **Typical loop** (from the Claude Code terminal):
   - edit `cases/<case>/input.yaml` and the SLURM script locally,
   - push to the cluster: `rsync -av cases/<case> caro:~/odt/run/<case>/`,
   - submit: `ssh caro 'cd ~/odt/run/<case> && sbatch slrmJob_array.sh'`,
   - watch: `ssh caro 'squeue --me'`,
   - pull post-processed results back (small files only):
     `rsync -av caro:~/odt/data/<case>/post/ results/<case>/`.

> **Which cluster for what?** Both are open to all DLR institutes and are functionally
> interchangeable for ODT (embarrassingly parallel realizations). Practical split: use
> whichever has queue availability; prefer **CARA** if you ever need the **A100 GPU**
> nodes; **CARO** (Göttingen, 64-core nodes) is a fine default for CPU realization
> sweeps. Keep one cluster as primary to avoid data sprawl.

### C2 — Run Claude Code **on** the cluster login node

If you'd rather have Claude Code sitting next to the data (large post-processing, editing
files in place on the cluster filesystem):

1. `ssh caro` (or `cara`), then install/run the Claude Code CLI on the login node and
   `cd` into your project/scratch directory there.
2. Everything is local to the cluster then — no rsync round-trips — but you're subject to
   login-node etiquette (no heavy compute on the login node; always `sbatch`).
3. Keep the git repo as the source of truth and `git pull`/`push` between your laptop and
   the cluster copy.

> Use **C1** for day-to-day development (edit locally, submit remotely) and **C2**
> occasionally for big on-cluster post-processing.

### C3 — Editing cluster files live from the desktop bridge (optional)

The Claude (Cowork) desktop app can mount a **local** folder for Claude to edit directly.
That is for local files, not the HPC. If you keep a synced/mounted copy of results
locally (e.g. via `rsync` or a mount), you can point Claude at that folder — but the
canonical way to reach CARO/CARA remains SSH (C1/C2).

---

## Part D — First-run checklist

- [ ] Claude Code installed, signed in, `claude --version` works.
- [ ] `~/research/odt-len/` created and laid out (Part B).
- [ ] `CLAUDE.md` at the repo root; `HANDOVER.md`, `docs/`, `papers/` present.
- [ ] 6 PDFs copied into `papers/` (see `papers/README.md`).
- [ ] ODT code cloned into `code/odt/`; note its git remote in `CLAUDE.md`.
- [ ] `~/.ssh/config` has `caro` and `cara`; `ssh caro` and `ssh cara` both work.
- [ ] `sinfo`/`squeue` succeed on both clusters.
- [ ] Fill in the placeholders in `docs/hpc-cheatsheet.md` (hostnames, account, modules).
- [ ] Open a session and paste `prompts/first-session-prompt.md`.

> **Placeholders to fill from your DLR onboarding** (not publicly documented, so left
> blank on purpose): CARO/CARA login hostnames, any DLR jump host, your SLURM
> account/project ID, partition/QOS names, and the module names for your compiler +
> MPI + Cantera + CMake.
