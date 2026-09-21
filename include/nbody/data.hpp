#pragma once

#include <cstddef>
#include <string>
#include <vector>

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