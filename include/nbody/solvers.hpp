#pragma once

#include "state.hpp"

class Solver {
public:
    Solver(double softening = 0.0, double gravitational_constant = 6.67430e-11);
    virtual ~Solver() = default;

    virtual void solve(
        const State& state,
        AccelerationField& acceleration
    ) const = 0;

protected:
    double softening_;
    double gravitational_constant_;
};

class NaiveSolver : public Solver {
public:
    NaiveSolver(
        double softening,
        double gravitational_constant = 6.67430e-11
    );

    void solve(
        const State& state,
        AccelerationField& acceleration
    ) const;
};

// Barnes-Hut and FMM solvers can implement the same Solver interface.