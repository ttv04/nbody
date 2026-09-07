#pragma once

#include "integrators.hpp"
#include "solvers.hpp"


class Simulation {
public:
    Simulation(
        Integrator& integrator,
        const Solver& solver,
        const State& initial_state
    );

    void step(double dt);

    State& state();
    const State& state() const;

private:
    State current_state_;
    AccelerationField acceleration_;
    Integrator& integrator_;
    const Solver& solver_;
};