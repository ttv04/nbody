#pragma once

#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"

void update_camera_3d(Camera3D& cam, float min_distance, float max_distance);
void draw_particles_3d(const Camera3D& cam, const State& state);
void draw_velocity_arrows_3d(const State& state, double scale);
void draw_trajectory_paths_3d(
    const std::vector<std::vector<std::vector<double>>>& displacement,
    const std::vector<Config>& styles,
    std::size_t frame,
    std::size_t max_segments = 2000
);
void draw_tree_bounds_3d(const Camera3D& cam, const std::vector<BHBound>& bounds);
