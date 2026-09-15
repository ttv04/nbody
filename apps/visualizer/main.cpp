#include <cstdlib>
#include <iostream>
#include "nbody/core.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "nbody/generators.hpp"
#include "raylib.h"
#include "viz2d.hpp"

#define DT 0.1
#define FPS 100

int main() {
    std::cout << "N-body simulation\n";

    LeapfrogIntegrator integrator;
    BarnesHutSolver solver(1.0, 1.0);
    UniformMassGalaxy2D generator(1000, 500, 10000, 500, 500, 1, 250);
    State initial_state = generator.getRandomState();
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

    InitWindow(1200, 800, "N-body simulation");
    SetTargetFPS(FPS);

    Camera2D cam{};
    cam.target = {500.0f, 500.0f};
    cam.offset = {GetScreenWidth() / 2.0f, GetScreenHeight() / 2.0f};
    cam.rotation = 0.0f;
    cam.zoom = 1.0f;

    double sim_time = 0.0;

    while (!WindowShouldClose()) {
        simulation.step(DT);
        sim_time += DT;
        update_camera(cam);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(cam);
        // draw_tree_bounds(cam, solver.last_bounds());
        draw_particles(cam, simulation.state());
        EndMode2D();

        draw_hud(DT, sim_time);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
