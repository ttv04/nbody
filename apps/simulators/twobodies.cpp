#include "nbody/core.hpp"
#include "nbody/data.hpp"
#include "nbody/state.hpp"
#include "argparse/argparse.hpp"
#include "sim_options.hpp"
#include <exception>
#include <cmath>
#include <numbers>
const double PI = std::numbers::pi;

int main(int argc, const char* argv[]) {
    argparse::ArgumentParser program("twobodies");

    add_sim_args(program);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    SimOptions o = sim_options_from(program);

    const double orbit_radius = 50.0;
    const double iden_mass = 50.0;
    const double iden_radius = 5.0;

    // For three equal masses on an equilateral triangle, gravity points toward
    // the center. Match it to v^2 / R, including the solver's softening:
    //
    //   v = R * sqrt(2 G m / (2 R^2 + soft^2)^(3/2))
    //
    // Gravity itself is the centripetal force; no additional force is needed.
    const double softened_side_squared =
        4.0 * orbit_radius * orbit_radius + o.soft * o.soft;
    const double angular_velocity = std::sqrt(
        2.0 * o.g * iden_mass
        / std::pow(softened_side_squared, 1.5)
    );
    const double orbit_speed = orbit_radius * angular_velocity;
    

    Body body1(2, iden_mass, {
        Config(255, 0, 0, 255, iden_radius)
    });

    body1.r[0] = orbit_radius * cos(0);
    body1.r[1] = orbit_radius * sin(0);
    body1.v[0] = orbit_speed * cos(PI / 2);
    body1.v[1] = orbit_speed * sin(PI / 2);

    Body body2(2, iden_mass, {
        Config(0, 255, 0, 255, iden_radius)
    });

    body2.r[0] = orbit_radius * cos(PI);
    body2.r[1] = orbit_radius * sin(PI);
    body2.v[0] = orbit_speed * cos(PI + PI / 2);
    body2.v[1] = orbit_speed * sin(PI + PI / 2);

    const State init_state({
        body1,
        body2
    });

    Simulation sim(*o.integrator, *o.solver, init_state);

    // Run the prestep before start_time
    int steps = (int)(o.t_start / o.dt); 
    
    for (int i = 0; i < steps; i++) {
        sim.step(o.dt);
    }

    // Now record
    Result result{
        "Two Bodies",
        o.integrator_name,
        o.solver_name,
        o.theta,
        o.max_leaf,
        2,
        2,
        o.g,
        o.soft,
        o.t_start,
        o.t_end,
        o.dt,
        static_cast<double>(o.seed),
        std::vector<double>({body1.m, body2.m}),
        std::vector<Config>({body1.config, body2.config}),
        std::vector<double>(),
        std::vector<double>(),
        std::vector<std::vector<std::vector<double>>>(),
        std::vector<std::vector<std::vector<double>>>(),
    };

    steps = (int)((o.t_end - o.t_start) / o.dt);
    for (int i = 0; i < steps; i++) {
        sim.step(o.dt);
        const State& state = sim.state();

        std::vector<std::vector<double>> displacement;
        std::vector<std::vector<double>> velocity;

        displacement.reserve(state.size());
        velocity.reserve(state.size());

        for (std::size_t body_index = 0; body_index < state.size(); body_index++) {
            const Body& body = state.bodies[body_index];
            displacement.push_back(body.r);
            velocity.push_back(body.v);
        }

        result.t.push_back(o.t_start + o.dt * (i + 1));
        result.energy.push_back(state.total_energy(o.g, o.soft));
        result.displacement.push_back(displacement);
        result.velocity.push_back(velocity);
    }
    
    save_result(result, o.output);
}