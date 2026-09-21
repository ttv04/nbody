#include "ring_experiment.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "nbody/core.hpp"
#include "nbody/data.hpp"
#include "nbody/generators.hpp"
#include "nbody/integrators.hpp"
#include "nbody/state.hpp"
#include "raylib.h"
#include "viz2d.hpp"

void run_ring_experiment(const RingExperiment& experiment) {
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
    const double initial_energy = initial_state.total_energy(G, SOFT);

    EulerIntegrator euler;
    EulerCromerIntegrator euler_cromer;
    LeapfrogIntegrator leapfrog;
    Integrator* integrators[] = {&euler, &euler_cromer, &leapfrog};
    const char* names[] = {"euler", "euler-cromer", "leapfrog"};
    // Filename stems the analysis notebook reads.
    const char* keys[] = {"euler", "euler_cromer", "leapfrog"};
    const int run_n = 3;

    const std::string title =
        "Ring orbiting a point mass (" + experiment.solver_label + ")";
    InitWindow(1200, 800, title.c_str());
    SetTargetFPS(60);

    Camera2D cam{};
    cam.offset = {
        GetScreenWidth() * 0.5f,
        GetScreenHeight() * 0.5f
    };
    cam.target = {0.0f, 0.0f};
    cam.rotation = 0.0f;
    cam.zoom = 5.0f;

    const std::filesystem::path out_dir(RESULTS_DIR);
    std::filesystem::create_directories(out_dir);

    // Stepping on a counter instead of accumulating DT keeps t exact.
    const long long total_steps = static_cast<long long>(T_END / DT);

    bool paused = false;
    bool show_bounds = false;
    bool aborted = false;

    // Every integrator replays the same initial state, back to back.
    for (int run = 0; run < run_n && !aborted; ++run) {
        Simulation sim(*integrators[run], experiment.solver, initial_state);

        long long step = 0;
        double time = 0.0;
        double energy = initial_energy;
        std::vector<double> sample_t{time};
        std::vector<double> sample_e{energy};

        while (step < total_steps) {
            if (WindowShouldClose()) {
                aborted = true;
                break;
            }
            if (IsKeyPressed(KEY_SPACE)) {
                paused = !paused;
            }
            if (IsKeyPressed(KEY_T)) {
                show_bounds = !show_bounds;
            }

            if (!paused) {
                for (int i = 0; i < STEPS_PER_FRAME && step < total_steps; ++i) {
                    sim.step(DT);
                    ++step;
                }
                time = step * DT;

                // O(N^2), so it stays out of the inner step loop.
                energy = sim.state().total_energy(G, SOFT);
                sample_t.push_back(time);
                sample_e.push_back(energy);
            }

            update_camera(cam);

            BeginDrawing();
            ClearBackground(BLACK);

            BeginMode2D(cam);
            draw_particles(cam, sim.state());
            if (show_bounds && experiment.bounds_source != nullptr) {
                draw_tree_bounds(cam, experiment.bounds_source->last_bounds());
            }
            EndMode2D();

            draw_hud(DT, time);
            DrawText(
                TextFormat(
                    "run %d/%d: %s   t: %.0f/%.0f",
                    run + 1, run_n, names[run], time, T_END
                ),
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
                experiment.bounds_source != nullptr
                    ? "[T] tree bounds  [space] pause"
                    : "[space] pause",
                10, 120, 18, GRAY
            );

            EndDrawing();
        }

        save_data(
            SimData(
                names[run],
                experiment.solver_label,
                initial_state.size(),
                sample_t,
                sample_e
            ),
            (out_dir
             / (experiment.file_prefix + "_" + keys[run] + ".json")).string()
        );
    }

    CloseWindow();
}
