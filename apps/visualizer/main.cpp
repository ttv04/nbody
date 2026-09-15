#include <cstdlib>
#include <iostream>
#include "nbody/core.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "nbody/generators.hpp"
#include "raylib.h"
#include "viz2d.hpp"
#include "viz3d.hpp"

#define DT 1
#define FPS 100

int main() {
    std::cout << "N-body simulation\n";

    LeapfrogIntegrator integrator;
    BarnesHutSolver solver(1.0, 1.0);
/*
    UniformMass2D generator(1000, 1000, 10, 0, 0);
    State initial_state = generator.getRandomState();
*/
    
    UniformMassGalaxy3D generator(1000, 500, 10000, 500, 500, 500, 1, 250);
    State initial_state = generator.getRandomState();
    Body sun{
        {500, 500, 500},
        {0.0, 0.0, 0.0},
        10000,
        Config{255, 165, 0, 255, 50}
    };

    // initial_state.bodies.push_back(sun);
    /*
    UniformMassGalaxy2D generator(1000, 500, 10000, 500, 500, 1, 250);
    State initial_state = generator.getRandomState();
    */
    /*
    Body sun{
        {500, 500},
        {0.0, 0.0},
        100000,
        Config{255, 165, 0, 255, 20}
    };

    initial_state.bodies.push_back(sun);
    */
    /*
    {
        Body(
            {780.0, 400.0},
            {0.0, 11.327},
            40000.0
        ),
        Body(
            {510.0, 555.885},
            {-9.809, -5.664},
            40000.0
        ),
        Body(
            {510.0, 244.115},
            {9.809, -5.664},
            40000.0
        )
    };
    */

    Simulation simulation(integrator, solver, initial_state);

    if (simulation.state().dim > 3) {
        std::cout << "ONLY SUPPORT DIMENSION <= 3\n";
        return 0;
    }

    const bool use_3d = simulation.state().dim == 3;

    InitWindow(1200, 800, use_3d ? "N-body simulation 3D" : "N-body simulation");
    SetTargetFPS(FPS);

    Camera2D cam2d{};
    cam2d.target = {0.0f, 0.0f};
    cam2d.offset = {GetScreenWidth() / 2.0f, GetScreenHeight() / 2.0f};
    cam2d.rotation = 0.0f;
    cam2d.zoom = 0.7f;

    Camera3D cam3d{};
    cam3d.position = {0.0f, -1400.0f, 900.0f};
    cam3d.target = {0.0f, 0.0f, 0.0f};
    cam3d.up = {0.0f, 1.0f, 0.0f};
    cam3d.fovy = 45.0f;
    cam3d.projection = CAMERA_PERSPECTIVE;

    double sim_time = 0.0;

    bool pause = false;
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            pause = !pause;
        }

        if (!pause) {
            simulation.step(DT);
            sim_time += DT;
        }

        if (use_3d) {
            update_camera_3d(cam3d);
        } else {
            update_camera(cam2d);
        }

        BeginDrawing();
        ClearBackground(BLACK);

        if (use_3d) {
            BeginMode3D(cam3d);
            draw_tree_bounds_3d(cam3d, solver.last_bounds());
            draw_particles_3d(cam3d, simulation.state());
            EndMode3D();
        } else {
            BeginMode2D(cam2d);
            draw_tree_bounds(cam2d, solver.last_bounds());
            draw_particles(cam2d, simulation.state());
            EndMode2D();
        }

        draw_hud(DT, sim_time);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
