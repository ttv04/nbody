# nbody

A C++20 gravitational N-body sandbox. The core library swaps force solvers and time integrators. Simulator apps write JSON trajectories; a raylib visualizer plays them back; a notebook plots energy stability.

## Features

- Direct pairwise gravity with softening (`NaiveSolver`, \(O(N^2)\), OpenMP over targets)
- Barnes–Hut orthtree gravity (`BarnesHutSolver`: quadtree in 2D, octree in 3D)
- Integrators: Euler, Euler–Cromer, leapfrog (KDK), velocity Verlet
- Built-in two-body and three-body circular setups
- JSON result files with time, energy, positions, and velocities
- 2D playback: camera, trajectories, velocity arrows, Barnes–Hut bounds
- Energy-drift plots in `analysis/analysis.ipynb`

Initial-condition generators exist in the library (uniform boxes, galaxies, point-mass rings) but do not yet have their own simulator binaries.

## Layout

```text
include/nbody/         public headers
src/core/              Simulation loop
src/state/             bodies, State, energy, acceleration buffer
src/solvers/           naive + Barnes–Hut
src/integrators/       time steppers
src/generators/        initial conditions
src/data/              JSON load / save
apps/simulators/       twobodies, threebodies, shared CLI
apps/visualizer/       raylib playback
apps/bh_tree_print/    dump a small Barnes–Hut tree
scripts/               compare integrators and write result JSON
analysis/              notebook + results/
```

The physics library does not depend on raylib. Only the visualizer links it.

## Requirements

- CMake 3.28+
- C++20 compiler (GCC 13 works)
- Git (used only if CMake has to download dependencies)
- On Linux / WSL: X11 and OpenGL development packages

```bash
sudo apt install build-essential cmake libx11-dev libxrandr-dev \
    libxinerama-dev libxcursor-dev libxi-dev libgl1-mesa-dev
```

CMake will use a system install if present, otherwise `FetchContent`:

| Library | Fallback |
|---|---|
| [raylib 5.5](https://github.com/raysan5/raylib) | visualizer only |
| [argparse v3.2](https://github.com/p-ranav/argparse) | CLI apps |
| [nlohmann/json v3.12.0](https://github.com/nlohmann/json) | result files |

OpenMP is used when CMake finds it.

## Build

```bash
cmake -B build
cmake --build build
```

Binaries land under `build/apps/`:

```text
build/apps/simulators/twobodies
build/apps/simulators/threebodies
build/apps/visualizer/visualizer
```

## How to use this repo

The usual loop is: **build → simulate → play back → plot energy**.

### 1. Run a simulation

Both simulators share the same flags. `t-end` is required. Defaults: `--dt 0.1`, `--g 1`, `--soft 0.5`, `--integrator velocity-verlet`, `--solver naive`.

```bash
# Two equal masses on a circular orbit
./build/apps/simulators/twobodies 2000 \
    --integrator leapfrog \
    --output analysis/results/twobodies_leapfrog.json

# Three equal masses on an equilateral orbit
./build/apps/simulators/threebodies 2000 \
    --integrator velocity-verlet \
    --solver naive \
    --output analysis/results/threebodies_velocity-verlet.json
```

Useful flags:

| Flag | Meaning |
|---|---|
| `--integrator` | `euler`, `euler-cromer`, `leapfrog`, `velocity-verlet` |
| `--solver` | `naive` or `barnes-hut` |
| `--dt` | timestep |
| `--t-start` | skip-and-record start time |
| `--g`, `--soft` | gravity constant and softening |
| `--theta`, `--max-leaf` | Barnes–Hut only |
| `--output` | JSON path (default: `<app>.json` in the cwd) |

Start with a modest `t-end`. Step count is `(t_end - t_start) / dt`, and every step is written to disk. `t-end=100000` with `dt=0.1` is a million frames per integrator.

For 2–3 bodies, OpenMP thread launch can dominate the force eval. Prefer:

```bash
OMP_NUM_THREADS=1 ./build/apps/simulators/twobodies 2000 --output analysis/results/twobodies_leapfrog.json
```

### 2. Compare integrators

```bash
# from the repo root
OMP_NUM_THREADS=1 ./scripts/compare_twobodies.sh 2000
OMP_NUM_THREADS=1 ./scripts/compare_threebodies.sh 2000
```

Each script runs all four integrators and writes:

```text
analysis/results/<prefix>_{euler,euler-cromer,leapfrog,velocity-verlet}.json
```

Pass a smaller first argument if you only want a smoke test. The scripts’ default `t-end` is currently `100000`.

### 3. Play a result back

The visualizer is **playback**, not a live sandbox. It needs a JSON file and is 2D-only.

```bash
./build/apps/visualizer/visualizer analysis/results/twobodies_leapfrog.json
./build/apps/visualizer/visualizer analysis/results/threebodies_leapfrog.json --fps 60 --trajectory-path-length 2000
```

| Control | Action |
|---|---|
| Space | pause / play |
| Left / Right | step one frame |
| R | restart |
| V | velocity arrows |
| P | trajectory paths |
| T | Barnes–Hut cell bounds (Barnes–Hut results only) |
| Drag | pan |
| Scroll | zoom |

`--trajectory-path-length -1` keeps the full trail.

### 4. Plot energy drift

From `analysis/`, so `results/...` resolves:

```bash
cd analysis
jupyter notebook analysis.ipynb
```

The notebook loads `results/<prefix>_<integrator>.json` and plots total energy plus relative drift \(|(E(t)-E_0)/E_0|\). Two-body runs are the clean integrator comparison; the equal-mass three-body orbit is dynamically unstable, so close encounters can dominate the curves.

`analysis/results/` is gitignored. Regenerate files after changing the solver or integrators.

## Force solvers

**Naive** sums Newton’s law between every pair. Softening is added to the distance so close encounters stay finite. The outer loop is one target body per thread, so each thread writes only its own acceleration row.

**Barnes–Hut** rebuilds an indexed orthtree each force evaluation. The root is an axis-aligned cube covering every body. Internal nodes have \(2^d\) children. A leaf stores at most `--max-leaf` particles (default 12). Distant groups are approximated when \(s/r < \theta\).

## Integrators

| Name | Order | Notes |
|---|---|---|
| `euler` | 1 | explicit, not symplectic; energy typically grows |
| `euler-cromer` | 1 | symplectic Euler; better long-term energy |
| `leapfrog` | 2 | KDK, symplectic |
| `velocity-verlet` | 2 | same family as leapfrog; first-step acceleration is computed once and reused |

RK4 is declared in the headers but not implemented yet.

## Writing another experiment

1. Add an app under `apps/simulators/` that builds a `State`, wraps it in `Simulation`, and records a `Result`.
2. Call `add_sim_args` / `sim_options_from` for the shared CLI.
3. Save with `save_result`.
4. Point the visualizer at the JSON, or add a `load_runs('<prefix>')` cell in the notebook.

Generators already in `include/nbody/generators.hpp`:

- `UniformMass2D` / `UniformMass3D`
- `UniformMassGalaxy2D` / `UniformMassGalaxy3D`
- `PointMassRing2D`

## Status

Working today: two-body and three-body runs, four integrators, naive and Barnes–Hut solvers, 2D JSON playback, energy plots.

Still open:

- RK4
- Simulator apps for the ring / galaxy generators
- 3D playback (`viz3d` exists; the visualizer app is 2D-only)
- Live (non-file) visualization
