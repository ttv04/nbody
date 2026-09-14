#include "nbody/solvers.hpp"

#include <cmath>

void accumulate_acceleration(
    AccelerationField& acceleration,
    size_t target,
    const std::vector<double>& target_pos,
    const std::vector<double>& source_pos,
    double source_mass,
    size_t dim,
    double softening,
    double grav
) {
    double r2 = 0.0;
    for (size_t d = 0; d < dim; d++) {
        const double displacement = source_pos[d] - target_pos[d];
        r2 += displacement * displacement;
    }

    const double softened_r2 = r2 + softening * softening;
    if (softened_r2 == 0.0) {
        return;
    }

    const double inverse_distance = 1.0 / std::sqrt(softened_r2);
    const double coefficient =
        grav * inverse_distance * inverse_distance * inverse_distance;

    for (size_t d = 0; d < dim; d++) {
        const double displacement = source_pos[d] - target_pos[d];
        acceleration(target, d) += coefficient * source_mass * displacement;
    }
}
