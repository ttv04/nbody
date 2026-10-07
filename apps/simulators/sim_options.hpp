#include <cstddef>
#include <string>
#include <memory>
#include "nbody/integrators.hpp"
#include "nbody/solvers.hpp"
#include "argparse/argparse.hpp"

struct SimOptions {
    double t_start, t_end, dt, g, soft, theta;
    size_t max_leaf;
    unsigned seed;
    std::unique_ptr<Integrator> integrator;
    std::string integrator_name;
    std::unique_ptr<Solver> solver; 
    std::string solver_name;
    std::string output;
};

void add_sim_args(argparse::ArgumentParser& program);
SimOptions sim_options_from(const argparse::ArgumentParser& program);
