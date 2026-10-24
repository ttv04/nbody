#include "viz2d.hpp"
#include "raymath.h"

void update_camera(Camera2D& cam, float min_zoom, float max_zoom) {
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) || IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
        const Vector2 delta = GetMouseDelta();
        cam.target.x -= delta.x / cam.zoom;
        cam.target.y -= delta.y / cam.zoom;
    }

    const float wheel = GetMouseWheelMove();
    if (wheel == 0.0f) {
        return;
    }

    const Vector2 before = GetScreenToWorld2D(GetMousePosition(), cam);
    cam.zoom *= (wheel > 0.0f) ? 1.1f : (1.0f / 1.1f);
    cam.zoom = Clamp(cam.zoom, min_zoom, max_zoom);
    const Vector2 after = GetScreenToWorld2D(GetMousePosition(), cam);
    cam.target.x += before.x - after.x;
    cam.target.y += before.y - after.y;
}

void draw_particles(const Camera2D& cam, const State& state) {
    for (const Body& body : state.bodies) {
        if (body.r.size() != 2) {
            continue;
        }

        const float min_world_radius = 0.8f / cam.zoom;
        DrawCircleV(
            Vector2{(float)body.r[0], (float)body.r[1]},
            fmaxf(body.config.radius, min_world_radius),
            Color{body.config.r, body.config.g, body.config.b, body.config.a}
        );
    }
}

void draw_trajectory_paths(
    const Camera2D& cam,
    const std::vector<std::vector<std::vector<double>>>& displacement,
    const std::vector<Config>& styles,
    std::size_t frame,
    std::size_t max_segments
) {
    if (frame == 0 || displacement.empty()) {
        return;
    }

    const std::size_t first_segment =
        frame > max_segments ? frame - max_segments + 1 : 1;
    const float thickness = 1.5f / cam.zoom;

    for (std::size_t body = 0; body < styles.size(); body++) {
        const Config& style = styles[body];
        const Color color{style.r, style.g, style.b, 130};

        for (std::size_t segment = first_segment; segment <= frame; segment++) {
            const auto& previous = displacement[segment - 1][body];
            const auto& current = displacement[segment][body];
            DrawLineEx(
                Vector2{
                    static_cast<float>(previous[0]),
                    static_cast<float>(previous[1])
                },
                Vector2{
                    static_cast<float>(current[0]),
                    static_cast<float>(current[1])
                },
                thickness,
                color
            );
        }
    }
}

void draw_velocity_arrows(const Camera2D& cam, const State& state, double scale) {
    const float thickness = 2.0f / cam.zoom;
    const float max_head_length = 10.0f / cam.zoom;
    const float max_head_width = 7.0f / cam.zoom;
    const Color color{120, 210, 255, 220};

    for (const Body& body : state.bodies) {
        if (body.r.size() != 2 || body.v.size() != 2) {
            continue;
        }

        const Vector2 position{
            static_cast<float>(body.r[0]),
            static_cast<float>(body.r[1])
        };
        const Vector2 velocity{
            static_cast<float>(body.v[0]),
            static_cast<float>(body.v[1])
        };
        const float speed = Vector2Length(velocity);
        if (speed <= 1e-6f) {
            continue;
        }

        const Vector2 direction = Vector2Scale(velocity, 1.0f / speed);
        const float displayed_radius = fmaxf(body.config.radius, 0.8f / cam.zoom);
        const Vector2 start = Vector2Add(
            position,
            Vector2Scale(direction, displayed_radius)
        );
        const Vector2 tip = Vector2Add(
            start,
            Vector2Scale(velocity, static_cast<float>(scale))
        );

        const float arrow_length = speed * static_cast<float>(scale);
        const float head_length = fminf(max_head_length, arrow_length * 0.4f);
        const float head_width = fminf(max_head_width, arrow_length * 0.3f);
        const Vector2 head_base = Vector2Subtract(
            tip,
            Vector2Scale(direction, head_length)
        );
        const Vector2 perpendicular{-direction.y, direction.x};
        const Vector2 left = Vector2Add(
            head_base,
            Vector2Scale(perpendicular, head_width * 0.5f)
        );
        const Vector2 right = Vector2Subtract(
            head_base,
            Vector2Scale(perpendicular, head_width * 0.5f)
        );

        DrawLineEx(start, head_base, thickness, color);
        DrawTriangle(tip, right, left, color);
    }
}

void draw_tree_bounds(const Camera2D& cam, const std::vector<BHBound>& bounds) {
    const float thickness = 1.0f / cam.zoom;
    const Color color{80, 200, 120, 160};

    for (const BHBound& bound : bounds) {
        if (bound.lo.size() < 2 || bound.hi.size() < 2) {
            continue;
        }

        const Rectangle rect{
            (float)bound.lo[0],
            (float)bound.lo[1],
            (float)(bound.hi[0] - bound.lo[0]),
            (float)(bound.hi[1] - bound.lo[1])
        };
        DrawRectangleLinesEx(rect, thickness, color);
    }
}

void draw_hud(double dt, double t) {
    DrawFPS(10, 10);
    DrawText(TextFormat("dt: %.4f", dt), 10, 32, 20, RAYWHITE);
    DrawText(TextFormat("t: %.2f", t), 10, 54, 20, RAYWHITE);
}
