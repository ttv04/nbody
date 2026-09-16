#include <cmath>
#include <cstddef>
#include <random>

#include "nbody/core.hpp"
#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"
#include "viz2d.hpp"

#define G 1.0
#define THETA 0.75
#define MPPL 1
#define SOFT 2.0
#define DT 0.08

int main() {
    const std::size_t particle_n = 2500;
    const double saturn_mass = 120000.0;
    const double particle_mass = 0.4;
    const float saturn_radius = 48.0f;
    const double cluster_x = 85.0;
    const double cluster_y = 0.0;
    const double cluster_radius = 14.0;

    State state(particle_n + 1, 2, particle_mass);

    state.bodies[0] = Body(
        {0.0, 0.0},
        {0.0, 0.0},
        saturn_mass,
        Config{255, 196, 92, 255, saturn_radius}
    );

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    std::uniform_real_distribution<double> tint(0.0, 1.0);

    std::size_t filled = 1;
    while (filled < state.size()) {
        const double nx = unit(rng);
        const double ny = unit(rng);
        if (nx * nx + ny * ny > 1.0) {
            continue;
        }

        const double x = cluster_x + nx * cluster_radius;
        const double y = cluster_y + ny * cluster_radius;
        const double r = std::sqrt(x * x + y * y);
        if (r < 1.0) {
            continue;
        }

        const double v = std::sqrt(G * saturn_mass / r);

        Body& body = state.bodies[filled];
        body.r[0] = x;
        body.r[1] = y;
        body.v[0] = -v * y / r;
        body.v[1] = v * x / r;
        body.m = particle_mass;

        const unsigned char shade = static_cast<unsigned char>(190 + 50 * tint(rng));
        body.config = Config{shade, shade, 255, 255, 0.35f};
        ++filled;
    }

    LeapfrogIntegrator integrator;
    BarnesHutSolver solver(SOFT, G, THETA, MPPL);
    Simulation sim{integrator, solver, state};

    InitWindow(1200, 800, "Saturn formation");
    SetTargetFPS(60);

    Camera2D cam{};
    cam.offset = {
        GetScreenWidth() * 0.5f,
        GetScreenHeight() * 0.5f
    };
    cam.target = {0.0f, 0.0f};
    cam.rotation = 0.0f;
    cam.zoom = 1.2f;

    bool paused = false;
    bool show_bounds = false;
    double time = 0.0;
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (IsKeyPressed(KEY_T)) {
            show_bounds = !show_bounds;
        }
        if (!paused) {
            sim.step(DT);
            Body& saturn = sim.state().bodies[0];
            saturn.r[0] = 0.0;
            saturn.r[1] = 0.0;
            saturn.v[0] = 0.0;
            saturn.v[1] = 0.0;
            time += DT;
        }

        update_camera(cam);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(cam);
        if (show_bounds) {
            draw_tree_bounds(cam, solver.last_bounds());
        }
        draw_particles(cam, sim.state());
        EndMode2D();

        draw_hud(DT, time);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
