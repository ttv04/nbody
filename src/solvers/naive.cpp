#include "nbody/solvers.hpp"

#include <cmath>

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

    for (std::size_t i = 0; i < state.size(); ++i) {
        const Body& body_a = state.bodies[i];

        for (std::size_t j = i + 1; j < state.size(); ++j) {
            const Body& body_b = state.bodies[j];

            double r2 = 0;
            for (std::size_t d = 0; d < state.dim; ++d) {
                const double displacement = body_b.r[d] - body_a.r[d];
                r2 += displacement * displacement;
            }

            const double softened_r2 =
                r2 + softening_ * softening_;

            if (softened_r2 == 0.0) {
                continue;
            }

            const double inverse_distance = 1.0 / std::sqrt(softened_r2);
            const double coefficient =
                gravitational_constant_
                * inverse_distance
                * inverse_distance
                * inverse_distance;

            for (std::size_t d = 0; d < state.dim; ++d) {
                const double displacement = body_b.r[d] - body_a.r[d];

                acceleration(i, d) +=
                    coefficient * body_b.m * displacement;
                acceleration(j, d) -=
                    coefficient * body_a.m * displacement;
            }
        }
    }
}