# Running fist-sampler on UH CARYA

Practical notes for building and running the SAM-3.0 (and other BES-Fluctuations) tasks
on the University of Houston [CARYA HPC cluster](https://www.uh.edu/research/rcdc/support-and-services/user-guide/getting-started-clusters.php).

## Cluster specs (confirmed via `sinfo`)

| Partition | Max walltime | Nodes |
|---|---|---|
| `batch` (default) | 14 days | 305 |
| `gpu`             | 7 days  | 40  |

Single-threaded FIST-Sampler jobs run happily on `batch` with `-n 1`.

---

## Prerequisites

- Active CougarNet ID on a project with a CARYA allocation.
- **UH VPN** if connecting from off-campus (mandatory).
- Windows users: use PuTTY / MobaXterm / XShell — WSL2 + VPN is flaky on CARYA.
- **VSCode Remote is not allowed on the HPC cluster.** Edit locally, sync with `rsync`/`scp`.

---

## 1. One-time setup

### 1.1 Connect
```bash
# Connect UH VPN first if off-campus.
ssh <cougarnet_id>@carya.rcdc.uh.edu
```

### 1.2 Work under `/project`, not `$HOME`

`$HOME` is a hard **10 GB** quota. This repo plus the MUSIC hypersurfaces and the
`T-μ` remap outputs will blow past that. Put everything under your group's project space:

```bash
cd /project/<PI_lastname>        # e.g., /project/vovchenko
mkdir -p $USER && cd $USER
```

### 1.3 Clone the repo (with submodules)

```bash
git clone --recursive https://github.com/vlvovch/fist-sampler.git
cd fist-sampler
# If you forgot --recursive:
#   git submodule update --init --recursive
```

### 1.4 Build — on a **compute node**, not the login node

Login nodes forbid compute-heavy work. Grab a short interactive allocation for the build:

```bash
salloc -t 1:00:00 -n 8 -N 1
ml avail                          # inspect what's offered
ml add gcc cmake                  # names may differ — e.g. gcc/11.2.0
# (intel-oneapi is also available if you prefer)

cd /project/.../fist-sampler
mkdir -p build && cd build
cmake ..
make -j8 BES-SAM3-cumulants       # or `make -j8` to build everything
exit                              # release the build allocation
```

Binary lands at `build/tasks/BESFluctuations/SAM3-cumulants/BES-SAM3-cumulants`.

### 1.5 Stage the large MUSIC hypersurfaces

`.gitignore` excludes `input/hydro/*/*/*.dat`, so surfaces aren't in git. Copy from your
workstation:

```bash
# From your Mac/laptop (substitute energies you need):
rsync -avP input/hydro/AuAu.27/ \
  <cougarnet_id>@carya.rcdc.uh.edu:/project/.../fist-sampler/input/hydro/AuAu.27/
```

Each hypersurface is ~250–420 MB.

---

## 2. Per-run workflow (one SLURM allocation per ensemble)

**Always use `tmux`.** From the RCDC guide:

> "If you get an interactive allocation … then disconnect from the cluster, for example
> by putting your laptop to sleep, your allocation will be terminated and your job killed."

### 2.1 Start tmux on the login node

```bash
tmux new -s sam3_27_GCE       # one session per (energy, ensemble)
```

### 2.2 Request an interactive allocation and launch

Inside tmux. The three `=1` flags guarantee **exactly 1 CPU core** is reserved:

```bash
salloc -p batch -t 14-00:00:00 --ntasks=1 --cpus-per-task=1 --nodes=1 --mem=4G

# Sanity-check the allocation:
echo "ntasks=$SLURM_NTASKS  cpus_per_task=$SLURM_CPUS_PER_TASK"
scontrol show job $SLURM_JOB_ID | grep -E 'NumCPUs|NumNodes'
# Expect: NumCPUs=1  NumNodes=1
ml add gcc cmake              # same modules you built with
cd /project/.../fist-sampler

# Convenience wrapper — takes <energy> [ensemble: GCE|B|BQS] [centrality]:
tasks/BESFluctuations/SAM3-cumulants/run-SAM3.sh 27 GCE
```

The script resolves paths from its own location, checks that the binary/input/
hypersurface exist, and redirects stdout+stderr to
`results/SAM3/AuAu.<energy>.<centrality>.EVHRG.<ensemble>.log`.

Detach with `Ctrl-b d`. Reconnect later:
```bash
ssh <cougarnet_id>@carya.rcdc.uh.edu
tmux attach -t sam3_27_GCE
```

### 2.3 Ensemble variants

The wrapper takes the ensemble as its second arg; output filenames are auto-suffixed
with `.GCE`, `.B`, or `.BQS` so they never collide. Run each in its own tmux session:

```bash
tmux new -s sam3_27_GCE   # -> salloc + run-SAM3.sh 27 GCE
tmux new -s sam3_27_B     # -> salloc + run-SAM3.sh 27 B
tmux new -s sam3_27_BQS   # -> salloc + run-SAM3.sh 27 BQS
```

All three run concurrently on separate nodes.

---

## 3. Monitoring & control

```bash
squeue -u $USER                    # your running jobs
scancel <jobid>                    # kill one
sbalance balance statement project <projectname>   # hours remaining
tmux ls                            # list your tmux sessions

# Peek at progress without reattaching:
tail -n 3 results/SAM3/AuAu.27.C0-5.EVHRG.GCE.log
grep '# Events:' results/SAM3/AuAu.27.C0-5.EVHRG.GCE.SAM3-cumulants.dat | head -1
```

The cumulant output files are refreshed every 1,000 events, so intermediate analysis is
always safe.

### Pulling results back to your workstation

```bash
# From your Mac:
rsync -avP <cougarnet_id>@carya.rcdc.uh.edu:/project/.../fist-sampler/results/SAM3/ \
  ./results/SAM3.carya/
```

---

## 4. Resource planning

| | 1 run | 3 ensembles × 1 energy | 3 ensembles × 7 BES energies |
|---|---|---|---|
| Max wall time | 14 days | 14 days (parallel) | 14 days (parallel) |
| CPU-hours | 336 | ~1,008 | ~7,000 |

In practice you'll `scancel` well before the wall — ~500k–1M events per ensemble is
plenty for publication statistics, which is more like **~50–150 CPU-h per run** based
on the ~800 events/min rate observed for 27 GeV on a typical x86_64 core.

Check your allocation balance before committing long runs:
```bash
sbalance balance statement project <projectname>
```

---

## 5. Gotchas / common pitfalls

- **Laptop sleep kills the job** unless you're inside `tmux`.
- **10 GB home quota** — output to `/project/...`, never `$HOME`.
- **VSCode Remote is disallowed** on CARYA.
- **Login node = no computation.** Always `salloc` before running anything non-trivial,
  including the build.
- **Re-use the same compiler modules** for build and run. Mismatched libstdc++
  versions can surface as runtime `GLIBCXX_*` errors.
- **Allocation award ID** may be required as `-A <award_id>` if your project isn't
  your default. If `salloc` complains, add `-A <projectname>`.
- **Hypersurface path is relative in input files.** Either `cd` to the project root
  before running, or override `--hypersurface_file=<absolute path>`.

---

## 6. Quick reference card

```bash
# one-time: connect + build
ssh <cougarnet_id>@carya.rcdc.uh.edu
cd /project/<PI_lastname>/$USER/fist-sampler
salloc -t 1:00:00 -n 8 -N 1; ml add gcc cmake
mkdir -p build && cd build && cmake .. && make -j8 BES-SAM3-cumulants
exit

# each run:
tmux new -s sam3_<energy>_<ens>
salloc -p batch -t 14-00:00:00 --ntasks=1 --cpus-per-task=1 --nodes=1 --mem=4G
ml add gcc cmake
tasks/BESFluctuations/SAM3-cumulants/run-SAM3.sh <energy> <GCE|B|BQS>
# Ctrl-b d to detach
```
