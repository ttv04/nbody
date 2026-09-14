#include "nbody/solvers.hpp"
#include "nbody/state.hpp"

#include <iostream>

int main() {
    State state({
        Body({0.0, 0.0}, {0.0, 0.0}, 1.0),
        Body({1.0, 0.0}, {0.0, 0.0}, 1.0),
        Body({0.0, 1.0}, {0.0, 0.0}, 1.0),
        Body({1.0, 1.0}, {0.0, 0.0}, 1.0),
        Body({0.1, 0.1}, {0.0, 0.0}, 1.0),
    });

    std::cout << "max_particles_per_leaf = 1\n";
    BHIndexedOrthoTree tree1(state, 1);
    tree1.print();

    std::cout << "\nmax_particles_per_leaf = 2\n";
    BHIndexedOrthoTree tree2(state, 2);
    tree2.print();

    State state2({
        Body({0}, {1}, 2),
        Body({50}, {1}, 2),
        Body({100}, {1}, 2)
    });
    BHIndexedOrthoTree tree3(state2, 1);
    tree3.print();
    return 0;
}
