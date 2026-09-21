#include "nbody/solvers.hpp"
#include "ring_experiment.hpp"

constexpr double THETA = 0.5;

int main() {
    BarnesHutSolver solver(SOFT, G, THETA, 1);

    run_ring_experiment({
        solver,
        "barnes-hut",
        "ring",
        &solver
    });

    return 0;
}
