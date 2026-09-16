#pragma once

#include "state.hpp"

class Solver;

class Integrator {
public:
    virtual ~Integrator() = default;

    virtual void step(
        State& state,
        const Solver& solver,
        AccelerationField& acceleration,
        double dt
    ) const = 0;
};

class EulerIntegrator : public Integrator {
public:
    void step(
        State& state,
        const Solver& solver,
        AccelerationField& acceleration,
        double dt
    ) const;
};

class EulerCromerIntegrator : public Integrator {
public:
    void step(
        State& state,
        const Solver& solver,
        AccelerationField& acceleration,
        double dt
    ) const;
};

class LeapfrogIntegrator: public Integrator {
public:
    void step(
        State& state,
        const Solver& solver,
        AccelerationField& acceleration,
        double dt
    ) const;
};

class RK4Integrator: public Integrator {
public:
    void step(
        State& state,
        const Solver& solver,
        AccelerationField& acceleration,
        double dt
    ) const;
};