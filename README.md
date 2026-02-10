# HS_CTC — Hard-Sphere Spherocylinder Collision Trajectory Code

Fortran DEM (Discrete Element Method) code for simulating binary collisions between two spherocylinder particles, used to study the **coefficient of tangential curvature (CTC)** in granular gas models.

---

## Physics Overview

Each simulation event fires two spherocylinders at each other with initial velocities sampled from Maxwell–Boltzmann distributions at prescribed translational and rotational temperatures. A Hertzian contact spring–dashpot force model with a given coefficient of restitution **α** is integrated until the particles separate. Post-collision scattering angles and energy partitioning are recorded.

**Particle geometry** — a spherocylinder of diameter *D* and cylinder length *L*:
```
AR = (L + D) / D    (aspect ratio; AR = 1 is a sphere)
```

**Contact force** (normal only):
```
F_N = K_N δ_n  +  C_N v_n          (Hertz spring + viscous damper)
K_N = (4/3) √(D/4) E*              C_N = 2 β √(m/2 · K_N)
```
where *E\** is the reduced Young's modulus and *β* is derived from the coefficient of restitution *α*.

**Integration** — velocity-Verlet with quaternion-based rotation update.

---

## Repository Structure

```
HS_CTC/
├── model/                  Fortran source
│   ├── main.f90            Main collision loop (OpenMP parallelised)
│   ├── read_input.f90      Input parsing (file + command-line)
│   ├── initialize.f90      Geometry, material props, file handles
│   ├── init_part.f90       Particle state initialisation, velocity sampling
│   ├── calc_force_dem.f90  Contact force, OVERLAP_PP, PROJECTED_AREA
│   ├── integrate_eom.f90   Velocity-Verlet time-stepper
│   ├── measure_dem.f90     Per-collision output (buffered)
│   ├── measure_final.f90   End-of-run statistics
│   └── mod/
│       ├── constants_mod.f90   π, SMALL_NUM, cross-product
│       ├── particles_mod.f90   Particle state arrays
│       ├── run_param_mod.f90   NTRY, NSAMPLES, dt
│       └── output_mod.f90      Output buffers, FLUSH_BUFFERS
├── build/
│   └── spherocylinder.mak  Makefile (gfortran + OpenMP)
├── run_parameter_sweep.py  Sequential parameter sweep (Python)
├── submit_array.sh         SLURM job-array script (HPC cluster)
└── system_input.dat        Runtime input file (geometry + material)
```

---

## Building

Requires **gfortran ≥ 9** with OpenMP support.

```bash
cd build
make -f spherocylinder.mak          # builds ../build/SphCyl
cp SphCyl ..                        # copy to project root
```

To rebuild from scratch:

```bash
make -f spherocylinder.mak clean && make -f spherocylinder.mak
```

**Compiler flags** (in `build/spherocylinder.mak`):
```
FCFLAGS = -g -fopenmp -O3
```

---

## Input File — `system_input.dat`

```
NSAMPLES                   ! number of collision events to record
DIA  LCYL  MASS            ! geometry (LCYL overridden by AR in CLI mode)
EYoung  GPoisson           ! material properties
ALPHA_PP                   ! coefficient of restitution
ToE  kTm  kTI              ! T(emperature)/E(nergy) mode, values
```

Example (used in legacy / test mode):
```
80000
1.D0    0.0D0    1.D0
8.9D9   0.3D0
0.80D0
T  1.0D0  1.0D0
```

---

## Running

### Command-line mode (recommended for sweeps)

```bash
export OMP_NUM_THREADS=20
./SphCyl  <alpha>  <kTm>  <kTI>  <AR>  <output_dir>
```

| Argument | Meaning |
|----------|---------|
| `alpha`  | Coefficient of restitution (0 < α ≤ 1) |
| `kTm`    | Translational temperature (kT/m units) |
| `kTI`    | Rotational temperature (kT/I units) |
| `AR`     | Aspect ratio; `LCYL = (AR−1)·DIA` |
| `output_dir` | Directory for output files (created if absent) |

`NSAMPLES`, `DIA`, `MASS`, `EYoung`, `GPoisson` are still read from `system_input.dat`.

**Example:**
```bash
./SphCyl 0.80 1.0 1.0 2.0 results/alpha0.80_r1.0_AR2.0
```

### Legacy mode (no CLI args)

```bash
cd run_directory          # directory that contains system_input.dat
/path/to/SphCyl
```

---

## Parameter Sweep

### Temperature ratio parameterisation

`kTI` is fixed at **1.0** to avoid singularities. The ratio
```
r = kTm / kTI = kTm
```
is swept in [0.1, 0.2, …, 2.0] (step 0.1, 20 values).

### Default sweep space

| Parameter | Range | Count |
|-----------|-------|-------|
| α (alpha) | 0.50 – 1.00, step 0.05 | 11 |
| r = kTm/kTI | 0.1 – 2.0, step 0.1 | 20 |
| AR | 1.0, 1.5, 2.0, 3.0, 4.0 | 5 |
| **Total** | | **1 100** |

### Sequential sweep (local)

```bash
python run_parameter_sweep.py --threads 20 --output-dir results

# Override parameter ranges at the command line:
python run_parameter_sweep.py --alpha 0.7 0.8 0.9 --ar 1.0 2.0
```

### SLURM job array (HPC cluster)

```bash
sbatch submit_array.sh          # launches 1100 independent jobs
squeue -u $USER                 # monitor progress
```

Each job uses 20 OpenMP threads on 1 node (20 CPU cores).
Adjust `--cpus-per-task` and `--time` in `submit_array.sh` as needed.

---

## Output Files

All files are written to the specified `output_dir`.

| File | Columns | Description |
|------|---------|-------------|
| `chi.txt` | b, χ, ψ | Impact parameter, scattering angle (translational), scattering angle (rotational) |
| `Ef.txt` | Et₀, Er1₀, Er2₀, Et_f, Er1_f, Er2_f, b_c | Pre/post-collision energies + contact impact parameter |
| `EnergyCons.txt` | E_f/E_0 | Energy conservation ratio (should be ≤ 1 for inelastic) |
| `NPhit.txt` | N | Number of contact points per collision |
| `PreRotEnergy.txt` | \|v_rel\|₀, \|v_rel\|_f | Initial and final relative speed |
| `projArea.txt` | A_proj | Projected area at first contact |
| `csx.txt` | σ | Collision cross-section estimate |

---

## Parallelisation

The outer collision loop (each NTRY trial is independent) is parallelised with **OpenMP**:

- Thread count set via `OMP_NUM_THREADS` environment variable.
- Each thread uses a unique random seed: `seed = 12345 + thread_id·100000 + NTRY`.
- Output is buffered per-thread (buffer size 1000) and flushed atomically.
- Expected parallel efficiency ~85 % on 20 cores (~17× speedup).

---

## Validation

Quick sanity checks after a run:

```bash
# Energy conservation: all values should be ≤ 1
awk '{if ($1 > 1.0001) print "FAIL: line " NR, $1}' EnergyCons.txt | head

# Line counts should equal NSAMPLES
wc -l chi.txt Ef.txt EnergyCons.txt
```

---

## License

This code was developed as part of a doctoral research project. All rights reserved.
