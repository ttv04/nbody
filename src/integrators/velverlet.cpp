#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"

VelocityVerletIntegrator::VelocityVerletIntegrator() : first(true) {}

void VelocityVerletIntegrator::step(
    State& state,
    const Solver& solver,
    AccelerationField& acceleration,
    double dt
) const {
    if (!first) {
        solver.solve(state, acceleration);
    }
    else {
        first = true;
    }

    for (size_t i = 0; i < state.size(); i++) {
        Body& body = state.bodies[i];
        for (size_t d = 0; d < state.dim; d++) {
            body.r[d] = body.r[d] + body.v[d] * dt + 0.5 * acceleration(i, d) * dt * dt;
            body.v[d] = body.v[d] + 0.5 * dt * acceleration(i, d);
        }
    }

    solver.solve(state, acceleration);

    for (size_t i = 0; i < state.size(); i++) {
        Body& body = state.bodies[i];
        for (size_t d = 0; d < state.dim; d++) {
            body.v[d] += 0.5 * dt * acceleration(i, d);
        }
    }

}