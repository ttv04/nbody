#include "viz2d.hpp"
#include "raymath.h"

void update_camera(Camera2D& cam) {
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
    cam.zoom = Clamp(cam.zoom, 0.05f, 1000.0f);
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
