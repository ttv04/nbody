#include <cstdlib>
#include <iostream>
#include "nbody/core.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "nbody/generators.hpp"
#include "raylib.h"

#define DT 0.01
#define FPS 100

int main() {
    std::cout << "N-body simulation\n";

    LeapfrogIntegrator integrator;
    BarnesHutSolver solver(1.0, 1.0);
    UniformMassGalaxy2D generator(1000, 500, 10000, 500, 500, 1, 250);
    State initial_state = generator.getRandomState();
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

    double sim_time = 0.0;

    while (!WindowShouldClose()) {
        simulation.step(DT);
        sim_time += DT;

        BeginDrawing();
        ClearBackground(BLACK);
        DrawFPS(10, 10);
        DrawText(TextFormat("dt: %.4f", DT), 10, 32, 20, RAYWHITE);
        DrawText(TextFormat("t: %.2f", sim_time), 10, 54, 20, RAYWHITE);

        for (const Body& body : simulation.state().bodies) {
            DrawCircleV(
                Vector2{(float)body.r[0], (float)body.r[1]},
                2.0f,
                RED
            );
        }

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
