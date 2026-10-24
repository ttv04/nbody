#include "argparse/argparse.hpp"
#include "nbody/data.hpp"
#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include "raylib.h"
#include "raymath.h"
#include "viz2d.hpp"
#include "viz3d.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void validate_result(const Result& result) {
    if (result.dim != 2 && result.dim != 3) {
        throw std::runtime_error(
            "This visualizer supports 2D and 3D results; file dimension is "
            + std::to_string(result.dim)
        );
    }
    if (result.N == 0 || result.displacement.empty()) {
        throw std::runtime_error("The result contains no bodies or frames");
    }
    if (result.masses.size() != result.N || result.style.size() != result.N) {
        throw std::runtime_error("Mass and style counts must match N");
    }
    if (result.t.size() != result.displacement.size()
        || result.velocity.size() != result.displacement.size()) {
        throw std::runtime_error("Time, displacement, and velocity frame counts must match");
    }

    for (std::size_t frame = 0; frame < result.displacement.size(); frame++) {
        if (result.displacement[frame].size() != result.N
            || result.velocity[frame].size() != result.N) {
            throw std::runtime_error("Each frame must contain N bodies");
        }
        for (std::size_t body = 0; body < result.N; body++) {
            if (result.displacement[frame][body].size() != result.dim
                || result.velocity[frame][body].size() != result.dim) {
                throw std::runtime_error("Each position and velocity must match dim");
            }
        }
    }
}

State make_state(const Result& result) {
    State state(result.N, result.dim);
    for (std::size_t body = 0; body < result.N; body++) {
        state.bodies[body].m = result.masses[body];
        state.bodies[body].config = result.style[body];
    }
    return state;
}

void apply_frame(const Result& result, std::size_t frame, State& state) {
    for (std::size_t body = 0; body < result.N; body++) {
        state.bodies[body].r = result.displacement[frame][body];
        state.bodies[body].v = result.velocity[frame][body];
    }
}

Camera2D make_camera(const Result& result) {
    double min_x = result.displacement.front().front()[0];
    double max_x = min_x;
    double min_y = result.displacement.front().front()[1];
    double max_y = min_y;

    for (const auto& frame : result.displacement) {
        for (const auto& position : frame) {
            min_x = std::min(min_x, position[0]);
            max_x = std::max(max_x, position[0]);
            min_y = std::min(min_y, position[1]);
            max_y = std::max(max_y, position[1]);
        }
    }

    double width = max_x - min_x;
    double height = max_y - min_y;
    const double span = std::max(width, height);
    const double extent = span > 0.0 ? span : 1.0;
    if (width <= 0.0) {
        width = extent;
    }
    if (height <= 0.0) {
        height = extent;
    }

    Camera2D camera{};
    camera.target = {
        static_cast<float>((min_x + max_x) * 0.5),
        static_cast<float>((min_y + max_y) * 0.5)
    };
    camera.offset = {
        GetScreenWidth() * 0.5f,
        GetScreenHeight() * 0.5f
    };
    camera.rotation = 0.0f;
    const double fitted = 0.8 * std::min(
        GetScreenWidth() / width,
        GetScreenHeight() / height
    );
    camera.zoom = std::isfinite(fitted) && fitted > 0.0
        ? static_cast<float>(fitted)
        : 1.0f;
    return camera;
}

Camera3D make_camera_3d(const Result& result, float& fitted_distance) {
    double min_x = result.displacement.front().front()[0];
    double max_x = min_x;
    double min_y = result.displacement.front().front()[1];
    double max_y = min_y;
    double min_z = result.displacement.front().front()[2];
    double max_z = min_z;

    for (const auto& frame : result.displacement) {
        for (const auto& position : frame) {
            min_x = std::min(min_x, position[0]);
            max_x = std::max(max_x, position[0]);
            min_y = std::min(min_y, position[1]);
            max_y = std::max(max_y, position[1]);
            min_z = std::min(min_z, position[2]);
            max_z = std::max(max_z, position[2]);
        }
    }

    const double dx = max_x - min_x;
    const double dy = max_y - min_y;
    const double dz = max_z - min_z;
    double radius = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(radius > 0.0)) {
        radius = 1.0;
    }

    const Vector3 target{
        static_cast<float>((min_x + max_x) * 0.5),
        static_cast<float>((min_y + max_y) * 0.5),
        static_cast<float>((min_z + max_z) * 0.5)
    };
    const float fovy = 45.0f;
    fitted_distance = static_cast<float>(radius / std::tan(fovy * DEG2RAD * 0.5));
    if (!(fitted_distance > 0.0f)) {
        fitted_distance = 10.0f;
    }

    const Vector3 direction = Vector3Normalize({0.35f, 0.5f, 0.79f});
    Camera3D camera{};
    camera.target = target;
    camera.position = Vector3Add(target, Vector3Scale(direction, fitted_distance));
    camera.up = {0.0f, 1.0f, 0.0f};
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;
    return camera;
}

void play_3d(
    const Result& result,
    unsigned int fps,
    std::size_t max_segments
) {
    const bool supports_tree = result.solver == "barnes-hut";
    InitWindow(1200, 800, result.title.c_str());
    SetTargetFPS(static_cast<int>(fps));

    float fitted_distance = 1.0f;
    Camera3D camera = make_camera_3d(result, fitted_distance);
    const float min_distance = fitted_distance * 0.02f;
    const float max_distance = fitted_distance * 40.0f;
    State state = make_state(result);
    BarnesHutSolver tree_solver(
        result.soft,
        result.g,
        result.solver_theta,
        result.solver_max_leaf
    );
    AccelerationField tree_acceleration;
    std::size_t frame = 0;
    std::size_t tree_frame = std::numeric_limits<std::size_t>::max();
    bool paused = false;
    bool show_velocity = false;
    bool show_tree = false;
    bool show_paths = false;

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (IsKeyPressed(KEY_R)) {
            frame = 0;
            paused = false;
        }
        if (IsKeyPressed(KEY_RIGHT) && frame + 1 < result.t.size()) {
            frame++;
            paused = true;
        }
        if (IsKeyPressed(KEY_LEFT) && frame > 0) {
            frame--;
            paused = true;
        }
        if (IsKeyPressed(KEY_V)) {
            show_velocity = !show_velocity;
        }
        if (supports_tree && IsKeyPressed(KEY_T)) {
            show_tree = !show_tree;
        }
        if (IsKeyPressed(KEY_P)) {
            show_paths = !show_paths;
        }

        apply_frame(result, frame, state);
        if (show_tree && tree_frame != frame) {
            tree_solver.solve(state, tree_acceleration);
            tree_frame = frame;
        }
        update_camera_3d(camera, min_distance, max_distance);

        BeginDrawing();
        ClearBackground(BLACK);
        BeginMode3D(camera);
        if (show_tree) {
            draw_tree_bounds_3d(camera, tree_solver.last_bounds());
        }
        if (show_paths) {
            draw_trajectory_paths_3d(
                result.displacement,
                result.style,
                frame,
                max_segments
            );
        }
        draw_particles_3d(camera, state);
        if (show_velocity) {
            draw_velocity_arrows_3d(state, 0.4);
        }
        EndMode3D();

        draw_hud(result.dt, result.t[frame]);
        DrawText(
            TextFormat("frame: %zu / %zu", frame + 1, result.t.size()),
            10,
            76,
            20,
            RAYWHITE
        );
        DrawText(
            TextFormat("integrator: %s", result.integrator.c_str()),
            10,
            98,
            20,
            RAYWHITE
        );
        DrawText(
            paused ? "paused | SPACE play | arrows step | R restart"
                   : "SPACE pause | arrows step | R restart",
            10,
            120,
            18,
            LIGHTGRAY
        );
        DrawText("left drag pans | right drag orbits | wheel zooms", 10, 142, 18, LIGHTGRAY);
        DrawText(
            TextFormat("V velocity arrows: %s", show_velocity ? "ON" : "OFF"),
            10,
            164,
            18,
            show_velocity ? SKYBLUE : LIGHTGRAY
        );
        DrawText(
            TextFormat("P trajectory paths: %s", show_paths ? "ON" : "OFF"),
            10,
            186,
            18,
            show_paths ? YELLOW : LIGHTGRAY
        );
        if (supports_tree) {
            DrawText(
                TextFormat("T tree bounds: %s", show_tree ? "ON" : "OFF"),
                10,
                208,
                18,
                show_tree ? GREEN : LIGHTGRAY
            );
        }
        EndDrawing();

        if (!paused && frame + 1 < result.t.size()) {
            frame++;
        } else if (frame + 1 == result.t.size()) {
            paused = true;
        }
    }

    CloseWindow();
}

}  // namespace

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("visualizer");
    program.add_argument("file").help("Path to a simulation result JSON file");
    program.add_argument("--fps")
        .scan<'u', unsigned int>()
        .default_value(60u)
        .help("Playback frames per second");
    program.add_argument("--trajectory-path-length")
        .scan<'i', int>()
        .default_value(2000)
        .help("Maximum trajectory segments per body; use -1 for complete history");
    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n' << program;
        return 1;
    }

    try {
        const unsigned int fps = program.get<unsigned int>("--fps");
        if (fps == 0) {
            throw std::runtime_error("--fps must be greater than zero");
        }

        const int requested_max_segments =
            program.get<int>("--trajectory-path-length");
        if (requested_max_segments == 0 || requested_max_segments < -1) {
            throw std::runtime_error(
                "--trajectory-path-length must be positive or -1"
            );
        }
        const std::size_t max_segments = requested_max_segments == -1
            ? std::numeric_limits<std::size_t>::max()
            : static_cast<std::size_t>(requested_max_segments);

        const Result result = load_result(program.get<std::string>("file"));
        validate_result(result);
        const bool supports_tree = result.solver == "barnes-hut";
        if (supports_tree && result.solver_max_leaf == 0) {
            throw std::runtime_error("Barnes-Hut max_leaf must be greater than zero");
        }
        if (result.dim == 3) {
            play_3d(result, fps, max_segments);
            return 0;
        }

        InitWindow(1200, 800, result.title.c_str());
        SetTargetFPS(static_cast<int>(fps));

        Camera2D camera = make_camera(result);
        const float min_zoom = camera.zoom * 0.05f;
        const float max_zoom = camera.zoom * 10000.0f;
        State state = make_state(result);
        BarnesHutSolver tree_solver(
            result.soft,
            result.g,
            result.solver_theta,
            result.solver_max_leaf
        );
        AccelerationField tree_acceleration;
        std::size_t frame = 0;
        std::size_t tree_frame = std::numeric_limits<std::size_t>::max();
        bool paused = false;
        bool show_velocity = false;
        bool show_tree = false;
        bool show_paths = false;

        while (!WindowShouldClose()) {
            if (IsKeyPressed(KEY_SPACE)) {
                paused = !paused;
            }
            if (IsKeyPressed(KEY_R)) {
                frame = 0;
                paused = false;
            }
            if (IsKeyPressed(KEY_RIGHT) && frame + 1 < result.t.size()) {
                frame++;
                paused = true;
            }
            if (IsKeyPressed(KEY_LEFT) && frame > 0) {
                frame--;
                paused = true;
            }
            if (IsKeyPressed(KEY_V)) {
                show_velocity = !show_velocity;
            }
            if (supports_tree && IsKeyPressed(KEY_T)) {
                show_tree = !show_tree;
            }
            if (IsKeyPressed(KEY_P)) {
                show_paths = !show_paths;
            }

            apply_frame(result, frame, state);
            if (show_tree && tree_frame != frame) {
                tree_solver.solve(state, tree_acceleration);
                tree_frame = frame;
            }
            update_camera(camera, min_zoom, max_zoom);

            BeginDrawing();
            ClearBackground(BLACK);
            BeginMode2D(camera);
            if (show_tree) {
                draw_tree_bounds(camera, tree_solver.last_bounds());
            }
            if (show_paths) {
                draw_trajectory_paths(
                    camera,
                    result.displacement,
                    result.style,
                    frame,
                    max_segments
                );
            }
            draw_particles(camera, state);
            if (show_velocity) {
                draw_velocity_arrows(camera, state, 5.0);
            }
            EndMode2D();

            draw_hud(result.dt, result.t[frame]);
            DrawText(
                TextFormat("frame: %zu / %zu", frame + 1, result.t.size()),
                10,
                76,
                20,
                RAYWHITE
            );
            DrawText(
                TextFormat("integrator: %s", result.integrator.c_str()),
                10,
                98,
                20,
                RAYWHITE
            );
            DrawText(
                paused ? "paused | SPACE play | arrows step | R restart"
                       : "SPACE pause | arrows step | R restart",
                10,
                120,
                18,
                LIGHTGRAY
            );
            DrawText(
                TextFormat("V velocity arrows: %s", show_velocity ? "ON" : "OFF"),
                10,
                142,
                18,
                show_velocity ? SKYBLUE : LIGHTGRAY
            );
            DrawText(
                TextFormat("P trajectory paths: %s", show_paths ? "ON" : "OFF"),
                10,
                164,
                18,
                show_paths ? YELLOW : LIGHTGRAY
            );
            if (supports_tree) {
                DrawText(
                    TextFormat("T tree bounds: %s", show_tree ? "ON" : "OFF"),
                    10,
                    186,
                    18,
                    show_tree ? GREEN : LIGHTGRAY
                );
            }
            EndDrawing();

            if (!paused && frame + 1 < result.t.size()) {
                frame++;
            } else if (frame + 1 == result.t.size()) {
                paused = true;
            }
        }

        CloseWindow();
    } catch (const std::exception& error) {
        std::cerr << "visualizer: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
