#include "nbody/solvers.hpp"
#include <omp.h>

Solver::Solver(double softening, double gravitational_constant)
    : softening_(softening),
      gravitational_constant_(gravitational_constant) {}

NaiveSolver::NaiveSolver(double softening, double gravitational_constant)
    : Solver(softening, gravitational_constant) {}

void NaiveSolver::solve(
    const State& state,
    AccelerationField& acceleration
) const {
    if (
        acceleration.body_count() != state.size()
        || acceleration.dimensions() != state.dim
    ) {
        acceleration.resize(state.size(), state.dim);
    } else {
        acceleration.reset();
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < state.size(); ++i) {
        const Body& body_a = state.bodies[i];
        for (std::size_t j = 0; j < state.size(); ++j) {
            if (i == j) {
                continue;
            }
            const Body& body_b = state.bodies[j];

            accumulate_acceleration(
                acceleration,
                i,
                body_a.r,
                body_b.r,
                body_b.m,
                state.dim,
                softening_,
                gravitational_constant_
            );
        }
    }
}
