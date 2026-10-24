#include "nbody/system.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

System::System(Body primary, double gravitational_constant)
    : primary_(std::move(primary)),
      gravitational_constant_(gravitational_constant) {
    require_gravitational_constant(gravitational_constant_);
    Body& body = std::get<Body>(primary_);
    require_body(body);
    body.r = {0.0, 0.0};
    body.v = {0.0, 0.0};
}

System::System(System primary, double gravitational_constant)
    : primary_(std::make_unique<System>(std::move(primary))),
      gravitational_constant_(gravitational_constant) {
    require_gravitational_constant(gravitational_constant_);
}

System::System(const System& other)
    : primary_(copy_node(other.primary_)),
      gravitational_constant_(other.gravitational_constant_) {
    satellites_.reserve(other.satellites_.size());
    for (const Satellite& satellite : other.satellites_) {
        satellites_.push_back(Satellite{
            copy_node(satellite.object),
            satellite.radius,
            satellite.theta,
            satellite.retrograde
        });
    }
}

System& System::operator=(const System& other) {
    if (this != &other) {
        System copy(other);
        *this = std::move(copy);
    }
    return *this;
}

System::System(System&& other) noexcept = default;
System& System::operator=(System&& other) noexcept = default;
System::~System() = default;

System& System::add(Body body, double radius, double theta, bool retrograde) {
    require_body(body);
    require_orbit(radius, theta);
    // (radius, theta) is this body's position. Its velocity is kept and the
    // orbit velocity is added in realize().
    body.r = {0.0, 0.0};
    satellites_.push_back(Satellite{std::move(body), radius, theta, retrograde});
    return *this;
}

System& System::add(System system, double radius, double theta, bool retrograde) {
    require_orbit(radius, theta);
    satellites_.push_back(Satellite{
        std::make_unique<System>(std::move(system)),
        radius,
        theta,
        retrograde
    });
    return *this;
}

double System::mass() const {
    double total = mass_of(primary_);
    for (const Satellite& satellite : satellites_) {
        total += mass_of(satellite.object);
    }
    return total;
}

State System::realize() const {
    std::vector<Body> bodies = local_bodies(primary_);
    const double primary_mass = mass_of(primary_);

    for (const Satellite& satellite : satellites_) {
        const double cos_theta = std::cos(satellite.theta);
        const double sin_theta = std::sin(satellite.theta);
        // Circular speed about the primary center of mass. The primary's
        // mass is already there: a body primary sits at the origin, and a
        // system primary was centered by its own realize(). Satellites do
        // not pull on the primary or on each other in this initial condition.
        const double speed = std::sqrt(
            gravitational_constant_ * primary_mass / satellite.radius
        );
        const double direction = satellite.retrograde ? -1.0 : 1.0;

        std::vector<Body> satellite_bodies = local_bodies(satellite.object);
        for (Body& body : satellite_bodies) {
            body.r[0] += satellite.radius * cos_theta;
            body.r[1] += satellite.radius * sin_theta;
            body.v[0] += direction * -speed * sin_theta;
            body.v[1] += direction * speed * cos_theta;
            bodies.push_back(std::move(body));
        }
    }

    // Place this system's own center of mass at the origin at rest, so a
    // parent orbit is measured from the barycenter rather than the primary.
    center_on_mass(bodies);
    return State(std::move(bodies));
}

void System::center_on_mass(std::vector<Body>& bodies) {
    double total_mass = 0.0;
    double com_x = 0.0;
    double com_y = 0.0;
    double com_vx = 0.0;
    double com_vy = 0.0;
    for (const Body& body : bodies) {
        total_mass += body.m;
        com_x += body.m * body.r[0];
        com_y += body.m * body.r[1];
        com_vx += body.m * body.v[0];
        com_vy += body.m * body.v[1];
    }
    if (!(total_mass > 0.0)) {
        return;
    }

    com_x /= total_mass;
    com_y /= total_mass;
    com_vx /= total_mass;
    com_vy /= total_mass;
    for (Body& body : bodies) {
        body.r[0] -= com_x;
        body.r[1] -= com_y;
        body.v[0] -= com_vx;
        body.v[1] -= com_vy;
    }
}

void System::require_body(const Body& body) {
    if (body.r.size() != 2 || body.v.size() != 2) {
        throw std::invalid_argument(
            "System is 2D only and requires bodies with 2 position and 2 velocity components"
        );
    }
    if (!std::isfinite(body.m) || body.m < 0.0) {
        throw std::invalid_argument("body mass must be finite and non-negative");
    }
}

void System::require_gravitational_constant(double gravitational_constant) {
    if (!std::isfinite(gravitational_constant) || gravitational_constant < 0.0) {
        throw std::invalid_argument(
            "gravitational constant must be finite and non-negative"
        );
    }
}

void System::require_orbit(double radius, double theta) {
    if (!std::isfinite(radius) || !(radius > 0.0)) {
        throw std::invalid_argument("orbital radius must be finite and positive");
    }
    if (!std::isfinite(theta)) {
        throw std::invalid_argument("orbital angle must be finite");
    }
}

System::Node System::copy_node(const Node& node) {
    if (const Body* body = std::get_if<Body>(&node)) {
        return *body;
    }
    const std::unique_ptr<System>& system = std::get<std::unique_ptr<System>>(node);
    return std::make_unique<System>(*system);
}

double System::mass_of(const Node& node) {
    if (const Body* body = std::get_if<Body>(&node)) {
        return body->m;
    }
    return std::get<std::unique_ptr<System>>(node)->mass();
}

std::vector<Body> System::local_bodies(const Node& node) {
    if (const Body* body = std::get_if<Body>(&node)) {
        return {*body};
    }
    State state = std::get<std::unique_ptr<System>>(node)->realize();
    return std::move(state.bodies);
}
