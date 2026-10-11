#include "nlohmann/json.hpp"
#include "nbody/data.hpp"
#include <fstream>
#include <stdexcept>
using json = nlohmann::json;

void save_data(
    const SimData& data,
    const std::string& filename) {

    json j;
    j["integrator"] = data.integrator;
    j["solver"] = data.solver;
    j["e"] = data.e;
    j["N"] = data.N;

    std::ofstream fo(filename);
    
    fo << j << std::endl;
}

void save_result(
    const Result& result,
    const std::string& filename) {

    json styles = json::array();
    for (const Config& config : result.style) {
        styles.push_back({
            {"color", {
                config.r,
                config.g,
                config.b,
                config.a
            }},
            {"radius", config.radius}
        });
    }

    json solver = {{"name", result.solver}};
    if (result.solver == "barnes-hut") {
        solver["theta"] = result.solver_theta;
        solver["max_leaf"] = result.solver_max_leaf;
    }

    json j = {
        {"name", result.title},
        {"g", result.g},
        {"soft", result.soft},
        {"dt", result.dt},
        {"t_start", result.t_start},
        {"t_end", result.t_end},
        {"seed", result.seed},
        {"integrator", result.integrator},
        {"solver", solver},
        {"dim", result.dim},
        {"N", result.N},
        {"masses", result.masses},
        {"style", styles},
        {"t", result.t},
        {"energy", result.energy},
        {"displacement", result.displacement},
        {"velocity", result.velocity},
    };
    
    std::ofstream fo(filename);

    fo << j << std::endl;

}

Result load_result(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Could not open result file: " + filename);
    }

    json j;
    input >> j;

    Result result{};
    result.title = j.at("name").get<std::string>();
    result.integrator = j.at("integrator").get<std::string>();

    const json& solver = j.at("solver");
    if (solver.is_object()) {
        result.solver = solver.at("name").get<std::string>();
        result.solver_theta = solver.value("theta", 0.75);
        result.solver_max_leaf =
            solver.value("max_leaf", static_cast<std::size_t>(12));
    } else {
        result.solver = solver.get<std::string>();
        result.solver_theta = 0.75;
        result.solver_max_leaf = 12;
    }

    result.N = j.at("N").get<std::size_t>();
    result.dim = j.at("dim").get<std::size_t>();
    result.g = j.at("g").get<double>();
    result.soft = j.at("soft").get<double>();
    result.t_start = j.at("t_start").get<double>();
    result.t_end = j.at("t_end").get<double>();
    result.dt = j.at("dt").get<double>();
    result.seed = j.at("seed").get<double>();
    result.masses = j.at("masses").get<std::vector<double>>();

    for (const json& style : j.at("style")) {
        const json& color = style.at("color");
        if (!color.is_array() || color.size() != 4) {
            throw std::runtime_error("Each style color must contain RGBA values");
        }

        const auto r = color.at(0).get<unsigned int>();
        const auto g = color.at(1).get<unsigned int>();
        const auto b = color.at(2).get<unsigned int>();
        const auto a = color.at(3).get<unsigned int>();
        if (r > 255 || g > 255 || b > 255 || a > 255) {
            throw std::runtime_error("Style color values must be between 0 and 255");
        }

        result.style.emplace_back(
            static_cast<unsigned char>(r),
            static_cast<unsigned char>(g),
            static_cast<unsigned char>(b),
            static_cast<unsigned char>(a),
            style.at("radius").get<float>()
        );
    }

    result.t = j.at("t").get<std::vector<double>>();
    result.energy = j.at("energy").get<std::vector<double>>();
    result.displacement =
        j.at("displacement").get<std::vector<std::vector<std::vector<double>>>>();
    result.velocity =
        j.at("velocity").get<std::vector<std::vector<std::vector<double>>>>();

    return result;
}