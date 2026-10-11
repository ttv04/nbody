#include "sim_options.hpp"
#include "argparse/argparse.hpp"
#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"
#include <memory>
#include <string>
#include <string_view>

namespace {

std::string default_output_name(const argparse::ArgumentParser& program) {
    const std::string usage = program.usage();
    constexpr std::string_view header = "Usage: ";
    const auto name_begin = header.size();
    const auto name_end = usage.find_first_of(" \n", name_begin);
    const auto name = usage.substr(
        name_begin,
        name_end == std::string::npos ? std::string::npos : name_end - name_begin
    );
    return name + ".json";
}

}

void add_sim_args(argparse::ArgumentParser& program) {
    program.add_argument("t-end")
        .scan<'g', double>()
        .help("The end time (or duration) of the simulation\n");
    
    program.add_argument("--integrator")
        .default_value(std::string("velocity-verlet"))
        .choices("velocity-verlet", "euler", "euler-cromer", "leapfrog")
        .help("Integrator selector. There are 4 integrators:\n"
            "- Euler: euler\n"
            "- Euler-Cromer: euler-cromer\n"
            "- Leapfrog: leapfrog\n"
            "- Velocity Verlet: velocity-verlet\n"
            "By default, Velocity Verlet will be selected"
        );
    program.add_argument("--solver")
        .default_value(std::string("naive"))
        .choices("naive", "barnes-hut")
        .help("Solver selector. There are 2 solvers:\n"
            "- Naive: naive\n"
            "\t A very simple force solver, with the calculation of a pairwise bodies's force. Time complexity O(N^2)"
            "- Barnes-Hut: barnes-hut\n"
            "\t An efficient force solver but the accuracy is slightly off. The accuracy depend on the choosen variable theta. Time complexity O(NlgN)"
            "By default, the naive method is selected"
        );
    
    program.add_argument("--output")
        .help("Specify the output name and directory")
        .default_value(default_output_name(program));

    program.add_argument("--theta")
        .scan<'g', double>()
        .default_value(0.75)
        .help("Barnes-Hut opening angle. Used only with --solver barnes-hut");
    
    program.add_argument("--max-leaf")
        .scan<'u', size_t>()
        .default_value(size_t{12})
        .help("Maximum bodies per leaf. Used only with --solver barnes-hut");

    program.add_argument("--g")
        .scan<'g', double>()
        .default_value(1.0)
        .help("The gravitational constant for the simulation");
    program.add_argument("--soft")
        .scan<'g', double>()
        .default_value(0.5)
        .help("Modify the softening factor (so that the interaction force between 2 very close bodies won't jump)\n");
    program.add_argument("--dt")
        .scan<'g', double>()
        .default_value(0.1)
        .help("Time step for each state transition\n");
    program.add_argument("--t-start")
        .scan<'g', double>()
        .default_value(0.0)
        .help("The start time of the simulation. By default it will start at t = 0\n");
    program.add_argument("--seed")
        .scan<'u', unsigned>()
        .default_value(42u)
        .help("Set the seed of the simulation");
}

SimOptions sim_options_from(const argparse::ArgumentParser& program) {
    std::unique_ptr<Integrator> integrator;
    std::unique_ptr<Solver> solver;

    auto lit_integrator = program.get<std::string>("--integrator");
    auto lit_solver = program.get<std::string>("--solver");
    auto lit_start_time = program.get<double>("--t-start");
    auto lit_end_time = program.get<double>("t-end");
    auto lit_g = program.get<double>("--g");
    auto lit_soft = program.get<double>("--soft");
    auto lit_seed = program.get<unsigned int>("--seed");
    auto lit_theta = program.get<double>("--theta");
    auto lit_max_leaf = program.get<size_t>("--max-leaf");
    auto lit_dt = program.get<double>("--dt");
    auto lit_output = program.get<std::string>("--output");
    const double g = lit_g;
    const double t_start = lit_start_time;
    const double t_end = lit_end_time;
    const double soft = lit_soft;
    const double dt = lit_dt;
    const double theta = lit_theta;
    const size_t max_leaf = lit_max_leaf;
    const unsigned int seed = lit_seed;
    const std::string output = lit_output;

    if (lit_integrator == "euler") {
        integrator = std::make_unique<EulerIntegrator>();
    }
    else if (lit_integrator == "euler-cromer") {
        integrator = std::make_unique<EulerCromerIntegrator>();
    }
    else if (lit_integrator == "leapfrog") {
        integrator = std::make_unique<LeapfrogIntegrator>();
    }
    else if (lit_integrator == "velocity-verlet") {
        integrator = std::make_unique<VelocityVerletIntegrator>();
    }

    if (lit_solver == "naive") {
        solver = std::make_unique<NaiveSolver>(soft, g);
    }
    else if (lit_solver == "barnes-hut") {
        solver = std::make_unique<BarnesHutSolver>(soft, g, theta, max_leaf);
    }

    SimOptions options = {
        t_start,
        t_end,
        dt,
        g,
        soft,
        theta,
        max_leaf,
        seed,
        std::move(integrator),
        lit_integrator,
        std::move(solver),
        lit_solver,
        output
    };

    return options;
}