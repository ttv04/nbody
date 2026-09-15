#pragma once

#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"

void update_camera(Camera2D& cam);
void draw_particles(const Camera2D& cam, const State& state);
void draw_tree_bounds(const Camera2D& cam, const std::vector<BHBound>& bounds);
void draw_hud(double dt, double t);
