#include "nbody/solvers.hpp"
#include "ring_experiment.hpp"

// Same ring setup as ring.cpp, but with exact pairwise forces, so the energy
// drift it reports comes from the integrator alone.
int main() {
    NaiveSolver solver(SOFT, G);

    run_ring_experiment({
        solver,
        "naive",
        "ring_naive"
    });

    return 0;
}
