#include "viz3d.hpp"
#include "raymath.h"

// Credit to Cursor Grok 4.6

namespace {

float world_radius_from_pixels(const Camera3D& cam, const Vector3& position, float pixel_radius) {
    const float dist = Vector3Distance(cam.position, position);
    const float half_height = GetScreenHeight() * 0.5f;
    if (half_height <= 0.0f) {
        return pixel_radius;
    }
    return pixel_radius * dist * tanf(cam.fovy * DEG2RAD * 0.5f) / half_height;
}

}  // namespace

void update_camera_3d(Camera3D& cam) {
    Vector3 offset = Vector3Subtract(cam.position, cam.target);

    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
        const Vector2 delta = GetMouseDelta();
        const float yaw = -delta.x * 0.005f;
        const float pitch = -delta.y * 0.005f;

        offset = Vector3RotateByAxisAngle(offset, cam.up, yaw);

        Vector3 forward = Vector3Subtract(cam.target, cam.position);
        if (Vector3LengthSqr(forward) > 0.0f) {
            const Vector3 right = Vector3Normalize(Vector3CrossProduct(forward, cam.up));
            const Vector3 pitched = Vector3RotateByAxisAngle(offset, right, pitch);
            const float alignment = fabsf(
                Vector3DotProduct(Vector3Normalize(pitched), Vector3Normalize(cam.up))
            );
            if (alignment < 0.995f) {
                offset = pitched;
            }
        }

        cam.position = Vector3Add(cam.target, offset);
    } else if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        const Vector2 delta = GetMouseDelta();
        const float dist = Vector3Length(offset);
        Vector3 forward = Vector3Subtract(cam.target, cam.position);
        if (Vector3LengthSqr(forward) > 0.0f && dist > 0.0f) {
            forward = Vector3Normalize(forward);
            const Vector3 right = Vector3Normalize(Vector3CrossProduct(forward, cam.up));
            const Vector3 up = Vector3Normalize(Vector3CrossProduct(right, forward));
            const float pan_scale = dist * 0.001f;
            const Vector3 move = Vector3Add(
                Vector3Scale(right, -delta.x * pan_scale),
                Vector3Scale(up, delta.y * pan_scale)
            );
            cam.target = Vector3Add(cam.target, move);
            cam.position = Vector3Add(cam.position, move);
        }
    }

    const float wheel = GetMouseWheelMove();
    if (wheel == 0.0f) {
        return;
    }

    float dist = Vector3Length(Vector3Subtract(cam.position, cam.target));
    dist *= (wheel > 0.0f) ? 0.9f : (1.0f / 0.9f);
    dist = Clamp(dist, 5.0f, 20000.0f);
    const Vector3 dir = Vector3Normalize(Vector3Subtract(cam.position, cam.target));
    cam.position = Vector3Add(cam.target, Vector3Scale(dir, dist));
}

void draw_particles_3d(const Camera3D& cam, const State& state) {
    for (const Body& body : state.bodies) {
        if (body.r.size() < 3) {
            continue;
        }

        const Vector3 position{
            (float)body.r[0],
            (float)body.r[1],
            (float)body.r[2]
        };
        DrawSphereEx(
            position,
            world_radius_from_pixels(cam, position, body.config.radius),
            6,
            8,
            Color{body.config.r, body.config.g, body.config.b, body.config.a}
        );
    }
}

void draw_tree_bounds_3d(const Camera3D& cam, const std::vector<BHBound>& bounds) {
    (void)cam;
    const Color color{80, 200, 120, 160};

    for (const BHBound& bound : bounds) {
        if (bound.lo.size() < 3 || bound.hi.size() < 3) {
            continue;
        }

        const Vector3 center{
            (float)(0.5 * (bound.lo[0] + bound.hi[0])),
            (float)(0.5 * (bound.lo[1] + bound.hi[1])),
            (float)(0.5 * (bound.lo[2] + bound.hi[2]))
        };
        DrawCubeWires(
            center,
            (float)(bound.hi[0] - bound.lo[0]),
            (float)(bound.hi[1] - bound.lo[1]),
            (float)(bound.hi[2] - bound.lo[2]),
            color
        );
    }
}
