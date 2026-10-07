#pragma once

#include <cstddef>
#include <string>
#include <vector>
#include "nbody/state.hpp"

struct Result {
    std::string title;
    std::string integrator;
    std::string solver;
    double solver_theta = 0.75;
    std::size_t solver_max_leaf = 12;
    std::size_t N;
    std::size_t dim;
    double g;
    double soft;
    double t_start;
    double t_end;
    double dt;
    double seed;
    std::vector<double> masses; // Should be of size N
    std::vector<Config> style; // Should be of size N
    std::vector<double> t; // Should be of size (t_end - t_start) / dt
    std::vector<double> energy; // Should be of size (t_end - t_start) / dt
    std::vector<std::vector<std::vector<double>>> displacement; // Should be of size N * dim * (t_end - t_start) / dt
    std::vector<std::vector<std::vector<double>>> velocity; // Should be of size N * dim * (t_end - t_start) / dt
};

struct SimData {
    std::string integrator;
    std::string solver;
    std::size_t N;
    std::vector<double> t;
    std::vector<double> e;
    SimData(
        const std::string& integrator, 
        const std::string& solver,
        std::size_t N,
        const std::vector<double>& t,
        const std::vector<double>& e
    ) : integrator(integrator),
        solver(solver),
        N(N),
        t(t),
        e(e) {}
};

void save_data(
    const SimData& data,
    const std::string& filename
);

void save_result(
    const Result& result,
    const std::string& filename
);

Result load_result(const std::string& filename);