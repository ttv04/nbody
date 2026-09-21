#pragma once

#include "nbody/solvers.hpp"
#include <string>

// Constants rather than macros: this header is included ahead of the nbody
// headers, and a macro named G collides with their `double G` parameters.
inline constexpr double G = 1.0;
inline constexpr double SOFT = 0.5;
inline constexpr double DT = 0.1;
inline constexpr double T_END = 2000.0;
// Steps batched per frame, so a 20000 step run does not take 20000 frames.
inline constexpr int STEPS_PER_FRAME = 10;
inline constexpr unsigned int SEED = 42u;

// One sweep of the ring setup over euler, euler-cromer and leapfrog, all
// replaying the same initial state. Each integrator writes
// <results>/<file_prefix>_<integrator>.json as soon as it finishes.
struct RingExperiment {
    const Solver& solver;
    // Recorded in the JSON, and shown in the window title.
    std::string solver_label;
    std::string file_prefix;
    // Only Barnes-Hut can report tree bounds, so [T] does nothing without it.
    const BarnesHutSolver* bounds_source = nullptr;
};

void run_ring_experiment(const RingExperiment& experiment);
