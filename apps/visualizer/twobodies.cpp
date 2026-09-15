#include "nbody/core.hpp"
#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"
#include "viz2d.hpp"

#define G 1.0
#define THETA 0.75
#define MPPL 12
#define SOFT 0.0
#define DT 0.1

int main() {
    Body body1{ {500, 500}, {0,  0.707}, 100, Config{230, 41, 55, 255, 8} };
    Body body2{ {600, 500}, {0, -0.707}, 100, Config{80, 180, 255, 255, 8} };

    State initital_state{
        {body1, body2}
    };

    LeapfrogIntegrator integrator;
    NaiveSolver solver(SOFT, G);
    Simulation sim{
        integrator,
        solver,
        initital_state
    };

    int h = GetScreenHeight();
    int w = GetScreenWidth();

    InitWindow(w, h, "2-body simulation");
    SetTargetFPS(60);

    Camera2D cam{
        Vector2{
            float(h) / 2.0f,
            float(w) / 2.0f
        },
        Vector2{
            0.0,
            0.0
        },
        0.0,
        1.0f
    };

    bool paused = false;
    double time = 0.0;
    while(!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (!paused) {
            sim.step(DT);
            time += DT;
        }

        update_camera(cam);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(cam);
        draw_particles(cam, sim.state());
        EndMode2D();

        draw_hud(DT, time);

        EndDrawing();
    }

    CloseWindow();
    
    return 0;
}
