#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"

/*
class LeapfrogIntegrator: public Integrator {
    public:
        void step(
            State& state,
            const Solver& solver,
            AccelerationField& acceleration,
            double dt
        ) const;
    };
*/

void LeapfrogIntegrator::step(
    State& state,
    const Solver& solver,
    AccelerationField& acceleration,
    double dt
) const {
    solver.solve(state, acceleration);

    // KDK form
    
    // v_{n + 1/2} and x_{n + 1}
    for (size_t i = 0; i < state.size(); i++) {
        Body& body = state.bodies[i];
        for (size_t d = 0; d < body.dimensions(); d++) {
            body.v[d] = body.v[d] + 0.5 * dt * acceleration(i, d);
            body.r[d] = body.r[d] + dt * body.v[d];
        }
    }

    // v_{n + 1}
    solver.solve(state, acceleration);
    for (size_t i = 0; i < state.size(); i++) {
        Body& body = state.bodies[i];
        for (size_t d = 0; d < body.dimensions(); d++) {
            body.v[d] = body.v[d] + 0.5 * dt * acceleration(i, d);
        }
    }
}