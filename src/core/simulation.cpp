#include "nbody/core.hpp"

Simulation::Simulation(
    Integrator& integrator,
    const Solver& solver,
    const State& initial_state
)
    : current_state_(initial_state),
      acceleration_(current_state_.size(), current_state_.dim),
      integrator_(integrator),
      solver_(solver) {}

State& Simulation::state() {
    return current_state_;
}

const State& Simulation::state() const {
    return current_state_;
}

void Simulation::step(double dt) {
    integrator_.step(current_state_, solver_, acceleration_, dt);
}