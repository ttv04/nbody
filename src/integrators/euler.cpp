#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"

void EulerIntegrator::step(
    State& state,
    const Solver& solver,
    AccelerationField& acceleration,
    double dt
) const {
    solver.solve(state, acceleration);

    for (std::size_t i = 0; i < state.size(); ++i) {
        Body& body = state.bodies[i];

        for (std::size_t d = 0; d < state.dim; ++d) {
            body.r[d] += body.v[d] * dt;
            body.v[d] += acceleration(i, d) * dt;
        }
    }
}