#include <cmath>
#include <cstddef>
#include <optional>

#include "nbody/core.hpp"
#include "nbody/generators.hpp"
#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"
#include "viz2d.hpp"

#define G 1.0
#define SOFT 0.5
#define DT 0.4
#define SEED 42u

double total_energy(const State& state) {
    double kinetic = 0.0;
    for (const Body& body : state.bodies) {
        double v2 = 0.0;
        for (double v : body.v) {
            v2 += v * v;
        }
        kinetic += 0.5 * body.m * v2;
    }

    // Softened the same way the solver softens the force.
    double potential = 0.0;
    for (std::size_t i = 0; i < state.size(); i++) {
        for (std::size_t j = i + 1; j < state.size(); j++) {
            double r2 = SOFT * SOFT;
            for (std::size_t d = 0; d < state.dim; d++) {
                const double dr = state.bodies[j].r[d] - state.bodies[i].r[d];
                r2 += dr * dr;
            }
            potential -= G * state.bodies[i].m * state.bodies[j].m / std::sqrt(r2);
        }
    }

    return kinetic + potential;
}

int main() {
    const std::size_t particle_n = 600;
    const double central_mass = 1000.0;
    const double particle_mass = 0.05;
    const double inner_radius = 40.0;
    const double outer_radius = 60.0;

    PointMassRing2D generator(
        particle_n,
        central_mass,
        particle_mass,
        inner_radius,
        outer_radius,
        G
    );
    const State initial_state = generator.getState(SEED);
    const double initial_energy = total_energy(initial_state);

    EulerIntegrator euler;
    EulerCromerIntegrator euler_cromer;
    LeapfrogIntegrator leapfrog;
    Integrator* integrators[] = {&euler, &euler_cromer, &leapfrog};
    const char* names[] = {"euler", "euler-cromer", "leapfrog"};
    int current = 2;

    NaiveSolver solver(SOFT, G);

    // Simulation holds references, so switching integrators rebuilds it.
    std::optional<Simulation> sim;
    sim.emplace(*integrators[current], solver, initial_state);

    InitWindow(1200, 800, "Ring orbiting a point mass");
    SetTargetFPS(60);

    Camera2D cam{};
    cam.offset = {
        GetScreenWidth() * 0.5f,
        GetScreenHeight() * 0.5f
    };
    cam.target = {0.0f, 0.0f};
    cam.rotation = 0.0f;
    cam.zoom = 5.0f;

    bool paused = false;
    double time = 0.0;
    double energy = initial_energy;
    int frame = 0;

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }

        int requested = current;
        if (IsKeyPressed(KEY_ONE)) {
            requested = 0;
        }
        if (IsKeyPressed(KEY_TWO)) {
            requested = 1;
        }
        if (IsKeyPressed(KEY_THREE)) {
            requested = 2;
        }

        if (requested != current || IsKeyPressed(KEY_R)) {
            current = requested;
            sim.emplace(*integrators[current], solver, initial_state);
            time = 0.0;
            energy = initial_energy;
        }

        if (!paused) {
            sim->step(DT);
            time += DT;
        }

        // O(N^2), so it does not run every frame.
        if (++frame % 10 == 0) {
            energy = total_energy(sim->state());
        }

        update_camera(cam);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(cam);
        draw_particles(cam, sim->state());
        EndMode2D();

        draw_hud(DT, time);
        DrawText(
            TextFormat("integrator: %s", names[current]),
            10, 76, 20, RAYWHITE
        );
        DrawText(
            TextFormat(
                "E: %.3f   drift: %+.2f%%",
                energy,
                100.0 * (energy - initial_energy) / std::fabs(initial_energy)
            ),
            10, 98, 20, RAYWHITE
        );
        DrawText(
            "[1] euler  [2] euler-cromer  [3] leapfrog  [R] restart  [space] pause",
            10, 120, 18, GRAY
        );

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
