# nbody

A C++ N-body gravity sandbox. The core library can swap force solvers and time integrators. A raylib visualizer draws the live simulation.

## Features

- Direct (naive) pairwise gravity with softening
- Euler, Euler–Cromer, and leapfrog (KDK) integrators
- Uniform 2D disk / ring initial conditions
- Interactive window with timestep and simulated time

## Layout

```text
include/nbody/   public headers
src/core/        simulation loop
src/state/       bodies, state, acceleration buffer
src/solvers/     force evaluators
src/integrators/ time steppers
src/generators/  initial conditions
apps/visualizer/ raylib app
```

## Requirements

- CMake 3.28+
- C++20 compiler (GCC 13 works)
- Git (used only if CMake has to download raylib)
- On Linux / WSL: X11 and OpenGL development packages

```bash
sudo apt install build-essential cmake libx11-dev libxrandr-dev \
    libxinerama-dev libxcursor-dev libxi-dev libgl1-mesa-dev
```

Raylib is **not** committed in this repo. CMake will:

1. Use a system / `CMAKE_PREFIX_PATH` install via `find_package(raylib)`
2. Otherwise download [raylib 5.5](https://github.com/raysan5/raylib) with `FetchContent` into the build directory

The physics library does not depend on raylib. Only `apps/visualizer` links it.

## Build

```bash
cmake -B build
cmake --build build
```

## Run

```bash
./build/apps/visualizer/visualizer
```

Close the window to quit. The overlay shows FPS, `dt`, and simulated time `t`.
