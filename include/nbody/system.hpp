#pragma once

#include "nbody/state.hpp"

#include <memory>
#include <variant>
#include <vector>

// 2D initial condition on circular orbits around a primary.
//
// realize() returns every system with its center of mass at the origin and
// its center-of-mass velocity equal to zero. A satellite, body or nested
// system, is placed so that its center of mass sits at (radius, theta) from
// the primary's center of mass. Every body in that satellite then receives
// the same circular-orbit velocity, sqrt(G * M_primary / radius), on top of
// the velocity it already had. Pass retrograde to orbit clockwise.
//
// realize() flattens the tree into a State. Every body must be 2D.
class System {
public:
    System(Body primary, double gravitational_constant);
    System(System primary, double gravitational_constant);

    System(const System& other);
    System& operator=(const System& other);
    System(System&& other) noexcept;
    System& operator=(System&& other) noexcept;
    ~System();

    System& add(Body body, double radius, double theta, bool retrograde = false);
    System& add(System system, double radius, double theta, bool retrograde = false);

    double mass() const;
    State realize() const;

private:
    using Node = std::variant<Body, std::unique_ptr<System>>;

    struct Satellite {
        Node object;
        double radius;
        double theta;
        bool retrograde;
    };

    static void require_body(const Body& body);
    static void require_gravitational_constant(double gravitational_constant);
    static void require_orbit(double radius, double theta);
    static Node copy_node(const Node& node);
    static double mass_of(const Node& node);
    static std::vector<Body> local_bodies(const Node& node);
    static void center_on_mass(std::vector<Body>& bodies);

    Node primary_;
    std::vector<Satellite> satellites_;
    double gravitational_constant_;
};
