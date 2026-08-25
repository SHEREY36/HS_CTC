# HS_CTC: standalone spherocylinder closure coefficients

HS_CTC performs binary classical-trajectory calculations for smooth, dissipative
spherocylinders. The closure pipeline estimates a Maxwellian routing baseline and
16 non-Maxwellian correction coefficients directly from collision data. It does
not fit against homogeneous-cooling, shear-flow, Fourier-flow, or other DEM
ensemble simulations.

## Closure

The tabulated model is

```text
f_tr = [3 theta/(3 theta + 2)] C_M(theta, alpha, AR)
       [1 + sum_a beta_a(theta, alpha, AR) X_a].
```

The 12 core features are `a2_tr`, `a2_rot`, `a11`, the six independent Gram
entries of `Pi`, `R`, and `Q`, and the three Gram entries of `q_tr` and `q_rot`.
The four separately flagged higher-order candidates are `a3_tr`, `a3_rot`,
`a21`, and `a12`. See `scripts/closure_common.py` for the exact normalization
and projection scores.

For each collision,

```text
delta_tr    = Etr_final_elastic - Etr_final_inelastic
delta_rot   = Erot_final_elastic - Erot_final_inelastic
delta_total = delta_tr + delta_rot
f_M         = sum(delta_tr) / sum(delta_total)
C_M         = f_M / [3 theta/(3 theta + 2)]
```

Per-event ratios are diagnostic only and are never averaged to obtain `f_M`.

## Build and run

Requirements are gfortran with OpenMP and Python 3 with NumPy.

```bash
make -C build clean all
python3 -m unittest discover -s tests -v
```

The extended executable interface is:

```text
build/SphCyl alpha kTm kTI AR output_dir [seed] [nsamples] [output_mode]
```

`output_mode` is `closure`, `legacy`, or `both`. The default is `both` for
backward compatibility. Closure production uses:

```bash
OMP_NUM_THREADS=20 build/SphCyl 0.8 1.0 1.0 2.0 \
  results/example/shard_000 12345 20000 closure
python3 scripts/validate_closure_run.py results/example/shard_000 --mark-success
python3 scripts/estimate_closure.py \
  --run-dir results/example/shard_000 --output-dir coefficients/example
```

`system_input.dat` supplies diameter, mass, and material properties. The command
line overrides alpha, temperatures, aspect ratio, seed, and sample count.

## Closure event schema

Every closure shard contains:

- `closure_events.bin`: little-endian float64 records with 38 columns;
- `metadata.txt`: schema version, exact column order, scales, parameters, seed,
  and expected record count;
- `validation.json` and `_SUCCESS` after validation.

The record contains both pre-collision particle velocities, both laboratory-frame
angular velocities, both particle axes, positive translational/rotational/total
dissipation, elastic and inelastic energies, elastic replay error, contact count,
impact geometry, and a unique event ID. The binary format is approximately 30.4
MB per 100,000 events. Event-ID-specific random streams make pre-collision states
identical across alpha for a fixed `(theta, AR, shard)` even under OpenMP.

## Full sweep on Slurm

The simulated grid has 12 inelastic alpha values (`0.50` through `0.95`, plus
`0.975` and `0.99`), 20 temperature ratios, and five aspect ratios: 1,200 nodes
per shard. Alpha 1 is not simulated because its routing ratio is `0/0`; the
aggregation script computes its smooth limit from the four near-elastic anchors.

After cloning on the HPC:

```bash
git clone git@github.com:SHEREY36/HS_CTC.git
cd HS_CTC
make -C build clean all
python3 -m venv .venv --system-site-packages
source .venv/bin/activate
python -c "import numpy"
mkdir -p logs results manifests coefficients closure_output
```

Pilot stage, 20,000 events per node:

```bash
python3 scripts/make_sweep_manifest.py \
  --stage pilot --samples 20000 --shard-id 0 --output manifests/pilot.csv
bash hpc/submit_manifest.sh manifests/pilot.csv hpc/sweep_array.slurm 50
```

When all pilot jobs pass, estimate each node. The submission record printed by
the wrapper contains all Slurm job IDs.

```bash
python3 scripts/check_manifest_success.py manifests/pilot.csv
bash hpc/submit_manifest.sh manifests/pilot.nodes.csv hpc/estimate_array.slurm 50
# After the estimation arrays finish:
python3 scripts/check_manifest_success.py manifests/pilot.nodes.csv --field coefficient_dir
```

Add the 80,000-event production shard and re-run the estimators:

```bash
python3 scripts/make_sweep_manifest.py \
  --stage production --samples 80000 --shard-id 1 --output manifests/production.csv
bash hpc/submit_manifest.sh manifests/production.csv hpc/sweep_array.slurm 50
python3 scripts/check_manifest_success.py manifests/production.csv
bash hpc/submit_manifest.sh manifests/production.nodes.csv hpc/estimate_array.slurm 50
# After the estimation arrays finish:
python3 scripts/check_manifest_success.py manifests/production.nodes.csv --field coefficient_dir
```

Aggregate the current estimates:

```bash
sbatch hpc/aggregate.slurm
```

If `closure_output/continuation_manifest.csv` contains tasks, submit it, rerun
the affected node estimators with `closure_output/continuation_nodes.csv`, and
aggregate again. Continuation shards contain
100,000 events, stop at one million events per node, and are requested when the
baseline or coefficient confidence-interval criterion fails. Tail contribution
remains a reported reliability diagnostic but does not automatically add jobs.

```bash
bash hpc/submit_manifest.sh closure_output/continuation_manifest.csv hpc/sweep_array.slurm 50
python3 scripts/check_manifest_success.py closure_output/continuation_manifest.csv
bash hpc/submit_manifest.sh closure_output/continuation_nodes.csv hpc/estimate_array.slurm 50
python3 scripts/check_manifest_success.py closure_output/continuation_nodes.csv \
  --field coefficient_dir
sbatch hpc/aggregate.slurm
```

For the final table, require all simulated nodes:

```bash
sbatch --export=ALL,REQUIRE_COMPLETE=1 hpc/aggregate.slurm
```

Final products are:

- `closure_coefficients_long.csv`: raw and alpha-limit coefficients with
  statistical and extrapolation-systematic errors;
- `closure_grid.npz`: dense `(alpha, theta, AR, coefficient)` values, confidence
  bounds, and separate uncertainty arrays;
- `qa_summary.csv`: precision and shard-stability state;
- `continuation_manifest.csv`: only nodes that still need events.

Use `squeue -u "$USER"` and `sacct -j JOB_ID` to monitor. Large event shards
belong in project or scratch storage, not Git.

## Git handoff

Existing local simulation directories are ignored. Stage source explicitly:

```bash
git status --short
git add model build scripts hpc tests README.md .gitignore requirements.txt \
  system_input.dat run_parameter_sweep.py submit_array.sh
git commit -m "Add standalone 16-term closure estimation pipeline"
git push origin master
```
