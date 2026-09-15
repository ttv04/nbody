#pragma once

#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"

void update_camera_3d(Camera3D& cam);
void draw_particles_3d(const Camera3D& cam, const State& state);
void draw_tree_bounds_3d(const Camera3D& cam, const std::vector<BHBound>& bounds);
