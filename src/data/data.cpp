#include "nlohmann/json.hpp"
#include "nbody/data.hpp"
#include <fstream>
using json = nlohmann::json;

void save_data(
    const SimData& data,
    const std::string& filename) {

    json j;
    j["integrator"] = data.integrator;
    j["solver"] = data.solver;
    j["t"] = data.t;
    j["e"] = data.e;
    j["N"] = data.N;

    std::ofstream fo(filename);
    
    fo << j << std::endl;
}