#pragma once

#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"

void update_camera(Camera2D& cam, float min_zoom, float max_zoom);
void draw_particles(const Camera2D& cam, const State& state);
void draw_trajectory_paths(
    const Camera2D& cam,
    const std::vector<std::vector<std::vector<double>>>& displacement,
    const std::vector<Config>& styles,
    std::size_t frame,
    std::size_t max_segments = 2000
);
void draw_velocity_arrows(const Camera2D& cam, const State& state, double scale);
void draw_tree_bounds(const Camera2D& cam, const std::vector<BHBound>& bounds);
void draw_hud(double dt, double t);
