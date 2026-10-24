#include "nbody/generators.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int radial_samples = 2048;
constexpr double two_pi = 2.0 * std::numbers::pi;

double nfw_mass_factor(double x) {
    return std::log1p(x) - x / (1.0 + x);
}

struct Equilibrium {
    std::vector<double> radius;
    std::vector<double> inward_acceleration;
    std::vector<double> potential;
    std::vector<double> halo_sigma;
    std::vector<double> bulge_sigma;
    std::vector<double> circular_speed;
    std::vector<double> radial_sigma;
    std::vector<double> azimuthal_sigma;
    std::vector<double> vertical_sigma;
    std::vector<double> azimuthal_speed;
};

double lookup(
    const std::vector<double>& radius,
    const std::vector<double>& values,
    double r
) {
    if (r <= radius.front()) {
        return values.front();
    }
    if (r >= radius.back()) {
        return values.back();
    }
    const auto upper = std::lower_bound(radius.begin(), radius.end(), r);
    const std::size_t hi = static_cast<std::size_t>(upper - radius.begin());
    const std::size_t lo = hi - 1;
    const double span = radius[hi] - radius[lo];
    const double fraction = span > 0.0 ? (r - radius[lo]) / span : 0.0;
    return values[lo] + fraction * (values[hi] - values[lo]);
}

void rotate_x(std::vector<double>& value, double angle) {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double y = value[1] * cosine - value[2] * sine;
    const double z = value[1] * sine + value[2] * cosine;
    value[1] = y;
    value[2] = z;
}

void rotate_z(std::vector<double>& value, double angle) {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double x = value[0] * cosine - value[1] * sine;
    const double y = value[0] * sine + value[1] * cosine;
    value[0] = x;
    value[1] = y;
}

void require_finite_nonnegative(double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and non-negative");
    }
}

void require_positive(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and positive");
    }
}

}  // namespace

CompositeGalaxy3D::CompositeGalaxy3D(
    GalaxyProfile profile,
    double gravitational_constant,
    double softening,
    double center_x,
    double center_y,
    double center_z
)
    : profile_(profile),
      gravitational_constant_(gravitational_constant),
      softening_(softening),
      center_x_(center_x),
      center_y_(center_y),
      center_z_(center_z) {
    if (!std::isfinite(gravitational_constant_) || gravitational_constant_ < 0.0) {
        throw std::invalid_argument("gravitational constant must be finite and non-negative");
    }
    require_finite_nonnegative(softening_, "softening");
    if (!std::isfinite(center_x_) || !std::isfinite(center_y_) || !std::isfinite(center_z_)) {
        throw std::invalid_argument("galaxy center must be finite");
    }
    if (!std::isfinite(profile_.inclination) || !std::isfinite(profile_.position_angle)) {
        throw std::invalid_argument("galaxy orientation must be finite");
    }

    const bool any_particles =
        profile_.halo_count > 0 || profile_.disk_count > 0 || profile_.bulge_count > 0;
    if (!any_particles) {
        throw std::invalid_argument("a galaxy needs at least one particle");
    }

    if (profile_.halo_count > 0) {
        require_positive(profile_.halo_mass, "halo mass");
        require_positive(profile_.halo_scale_radius, "halo scale radius");
        require_positive(profile_.halo_virial_radius, "halo virial radius");
    }
    if (profile_.bulge_count > 0) {
        require_positive(profile_.bulge_mass, "bulge mass");
        require_positive(profile_.bulge_scale_radius, "bulge scale radius");
        require_positive(profile_.bulge_outer_radius, "bulge outer radius");
    }
    if (profile_.disk_count > 0) {
        require_positive(profile_.disk_mass, "disk mass");
        require_positive(profile_.disk_scale_radius, "disk scale radius");
        require_positive(profile_.disk_scale_height, "disk scale height");
        require_positive(profile_.disk_outer_radius, "disk outer radius");
        require_positive(profile_.toomre_q, "Toomre Q");
    }
}

State CompositeGalaxy3D::getRandomState() {
    std::random_device random;
    return getState(random());
}

State CompositeGalaxy3D::getState(unsigned int seed) {
    const GalaxyProfile& profile = profile_;
    const double grav = gravitational_constant_;
    const double softening2 = softening_ * softening_;

    const bool has_halo = profile.halo_count > 0;
    const bool has_bulge = profile.bulge_count > 0;
    const bool has_disk = profile.disk_count > 0;

    const double concentration = has_halo
        ? profile.halo_virial_radius / profile.halo_scale_radius
        : 1.0;
    const double nfw_virial_factor = has_halo ? nfw_mass_factor(concentration) : 1.0;
    const double halo_density_scale = has_halo
        ? profile.halo_mass / (
            4.0 * std::numbers::pi
            * profile.halo_scale_radius * profile.halo_scale_radius * profile.halo_scale_radius
            * nfw_virial_factor
        )
        : 0.0;

    const double bulge_enclosed_fraction = has_bulge
        ? std::pow(
            profile.bulge_outer_radius
                / (profile.bulge_outer_radius + profile.bulge_scale_radius),
            2.0
        )
        : 1.0;
    const double bulge_infinite_mass = has_bulge
        ? profile.bulge_mass / bulge_enclosed_fraction
        : 0.0;

    const double disk_x_outer = has_disk
        ? profile.disk_outer_radius / profile.disk_scale_radius
        : 1.0;
    const double disk_enclosed_fraction = has_disk
        ? 1.0 - std::exp(-disk_x_outer) * (1.0 + disk_x_outer)
        : 1.0;
    const double disk_infinite_mass = has_disk
        ? profile.disk_mass / disk_enclosed_fraction
        : 0.0;
    const double disk_central_density = has_disk
        ? disk_infinite_mass / (
            2.0 * std::numbers::pi
            * profile.disk_scale_radius * profile.disk_scale_radius
        )
        : 0.0;

    auto halo_mass_within = [&](double r) {
        if (!has_halo || r <= 0.0) {
            return 0.0;
        }
        if (r >= profile.halo_virial_radius) {
            return profile.halo_mass;
        }
        return profile.halo_mass
            * nfw_mass_factor(r / profile.halo_scale_radius)
            / nfw_virial_factor;
    };
    auto bulge_mass_within = [&](double r) {
        if (!has_bulge || r <= 0.0) {
            return 0.0;
        }
        if (r >= profile.bulge_outer_radius) {
            return profile.bulge_mass;
        }
        const double fraction = r / (r + profile.bulge_scale_radius);
        return bulge_infinite_mass * fraction * fraction;
    };
    auto halo_density = [&](double r) {
        if (!has_halo || r <= 0.0 || r >= profile.halo_virial_radius) {
            return 0.0;
        }
        const double x = r / profile.halo_scale_radius;
        return halo_density_scale / (x * (1.0 + x) * (1.0 + x));
    };
    auto bulge_density = [&](double r) {
        if (!has_bulge || r <= 0.0 || r >= profile.bulge_outer_radius) {
            return 0.0;
        }
        const double radius_plus_scale = r + profile.bulge_scale_radius;
        return bulge_infinite_mass / (2.0 * std::numbers::pi)
            * profile.bulge_scale_radius
            / (r * radius_plus_scale * radius_plus_scale * radius_plus_scale);
    };
    auto disk_surface_density = [&](double cylindrical_radius) {
        if (!has_disk || cylindrical_radius >= profile.disk_outer_radius) {
            return 0.0;
        }
        return disk_central_density * std::exp(-cylindrical_radius / profile.disk_scale_radius);
    };
    auto disk_speed_squared = [&](double cylindrical_radius) {
        if (!has_disk || cylindrical_radius <= 0.0) {
            return 0.0;
        }
        const double y = cylindrical_radius / (2.0 * profile.disk_scale_radius);
        if (y < 1.0e-8) {
            return 0.0;
        }
        const double bessel_term =
            std::cyl_bessel_i(0.0, y) * std::cyl_bessel_k(0.0, y)
            - std::cyl_bessel_i(1.0, y) * std::cyl_bessel_k(1.0, y);
        const double newtonian =
            4.0 * std::numbers::pi * grav * disk_central_density
            * profile.disk_scale_radius * y * y * bessel_term;
        const double radius2 = cylindrical_radius * cylindrical_radius;
        const double softened = radius2 / (radius2 + softening2);
        return std::max(0.0, newtonian) * softened * std::sqrt(softened);
    };

    double inner_radius = 0.0;
    double outer_radius = 0.0;
    if (has_halo) {
        inner_radius = profile.halo_scale_radius;
        outer_radius = std::max(outer_radius, profile.halo_virial_radius);
    }
    if (has_bulge) {
        inner_radius = inner_radius == 0.0
            ? profile.bulge_scale_radius
            : std::min(inner_radius, profile.bulge_scale_radius);
        outer_radius = std::max(outer_radius, profile.bulge_outer_radius);
    }
    if (has_disk) {
        inner_radius = inner_radius == 0.0
            ? profile.disk_scale_radius
            : std::min(inner_radius, profile.disk_scale_radius);
        outer_radius = std::max(outer_radius, profile.disk_outer_radius);
    }
    inner_radius *= 1.0e-4;
    outer_radius *= 1.05;

    Equilibrium equilibrium;
    equilibrium.radius.resize(radial_samples);
    equilibrium.inward_acceleration.resize(radial_samples);
    equilibrium.potential.resize(radial_samples);
    equilibrium.halo_sigma.resize(radial_samples);
    equilibrium.bulge_sigma.resize(radial_samples);
    equilibrium.circular_speed.resize(radial_samples);
    equilibrium.radial_sigma.resize(radial_samples);
    equilibrium.azimuthal_sigma.resize(radial_samples);
    equilibrium.vertical_sigma.resize(radial_samples);
    equilibrium.azimuthal_speed.resize(radial_samples);

    const double log_inner = std::log(inner_radius);
    const double log_outer = std::log(outer_radius);
    for (int i = 0; i < radial_samples; i++) {
        const double fraction = static_cast<double>(i) / (radial_samples - 1);
        const double radius = std::exp(log_inner + (log_outer - log_inner) * fraction);
        equilibrium.radius[i] = radius;

        const double radius2 = radius * radius;
        const double softened_radius = std::pow(radius2 + softening2, 1.5);
        const double spherical_mass = halo_mass_within(radius) + bulge_mass_within(radius);
        const double spherical_acceleration =
            grav * spherical_mass * radius / softened_radius;
        const double disk_acceleration = disk_speed_squared(radius) / radius;
        equilibrium.inward_acceleration[i] = spherical_acceleration + disk_acceleration;
        equilibrium.circular_speed[i] = std::sqrt(
            std::max(0.0, equilibrium.inward_acceleration[i] * radius)
        );
    }

    equilibrium.potential.back() =
        -grav * (profile.halo_mass * has_halo + profile.bulge_mass * has_bulge
            + profile.disk_mass * has_disk)
        / equilibrium.radius.back();
    for (int i = radial_samples - 2; i >= 0; i--) {
        const double dr = equilibrium.radius[i + 1] - equilibrium.radius[i];
        const double mean_acceleration = 0.5 * (
            equilibrium.inward_acceleration[i] + equilibrium.inward_acceleration[i + 1]
        );
        equilibrium.potential[i] =
            equilibrium.potential[i + 1] - mean_acceleration * dr;
    }

    auto dispersion_from_jeans = [&](auto density, std::vector<double>& sigma) {
        double accumulated = 0.0;
        sigma.back() = 0.0;
        for (int i = radial_samples - 2; i >= 0; i--) {
            const double dr = equilibrium.radius[i + 1] - equilibrium.radius[i];
            const double density_here = density(equilibrium.radius[i]);
            const double density_out = density(equilibrium.radius[i + 1]);
            accumulated += 0.5 * dr * (
                density_here * equilibrium.inward_acceleration[i]
                + density_out * equilibrium.inward_acceleration[i + 1]
            );
            sigma[i] = density_here > 0.0
                ? std::sqrt(std::max(0.0, accumulated / density_here))
                : 0.0;
        }
    };
    if (has_halo) {
        dispersion_from_jeans(halo_density, equilibrium.halo_sigma);
    }
    if (has_bulge) {
        dispersion_from_jeans(bulge_density, equilibrium.bulge_sigma);
    }

    if (has_disk) {
        std::vector<double> omega(radial_samples, 0.0);
        std::vector<double> kappa(radial_samples, 0.0);
        for (int i = 0; i < radial_samples; i++) {
            omega[i] = equilibrium.circular_speed[i] / equilibrium.radius[i];
        }
        for (int i = 1; i < radial_samples - 1; i++) {
            const double dr = equilibrium.radius[i + 1] - equilibrium.radius[i - 1];
            const double omega_gradient = (omega[i + 1] - omega[i - 1]) / dr;
            const double kappa_squared =
                4.0 * omega[i] * omega[i]
                + 2.0 * omega[i] * equilibrium.radius[i] * omega_gradient;
            const double floor = 0.05 * omega[i] * omega[i];
            kappa[i] = std::sqrt(std::max(kappa_squared, floor));
        }
        kappa.front() = kappa[1];
        kappa.back() = kappa[radial_samples - 2];

        std::vector<double> surface_sigma_squared(radial_samples, 0.0);
        for (int i = 0; i < radial_samples; i++) {
            const double surface = disk_surface_density(equilibrium.radius[i]);
            const double kappa_here = std::max(kappa[i], 1.0e-8);
            equilibrium.radial_sigma[i] =
                profile.toomre_q * 3.36 * grav * surface / kappa_here;
            const double sheet_sigma = std::sqrt(
                std::numbers::pi * grav * surface * profile.disk_scale_height
            );
            const double in_plane_sigma = equilibrium.radial_sigma[i];
            equilibrium.vertical_sigma[i] = std::max(sheet_sigma, 0.6 * in_plane_sigma);
            const double omega_here = std::max(omega[i], 1.0e-8);
            equilibrium.azimuthal_sigma[i] =
                equilibrium.radial_sigma[i] * kappa_here / (2.0 * omega_here);
            surface_sigma_squared[i] =
                surface * equilibrium.radial_sigma[i] * equilibrium.radial_sigma[i];
        }

        for (int i = 1; i < radial_samples - 1; i++) {
            const double radial_sigma = equilibrium.radial_sigma[i];
            const double circular = equilibrium.circular_speed[i];
            if (radial_sigma <= 0.0 || circular <= 0.0) {
                equilibrium.azimuthal_speed[i] = circular;
                continue;
            }
            const double dr = equilibrium.radius[i + 1] - equilibrium.radius[i - 1];
            const double derivative =
                (surface_sigma_squared[i + 1] - surface_sigma_squared[i - 1]) / dr;
            const double logarithmic_slope = surface_sigma_squared[i] > 0.0
                ? equilibrium.radius[i] * derivative / surface_sigma_squared[i]
                : 0.0;
            const double anisotropy =
                equilibrium.azimuthal_sigma[i] / radial_sigma;
            const double lag_coefficient =
                -logarithmic_slope - (1.0 - anisotropy * anisotropy);
            const double streaming_squared =
                circular * circular - radial_sigma * radial_sigma * lag_coefficient;
            const double clamped = std::clamp(streaming_squared, 0.0, circular * circular);
            equilibrium.azimuthal_speed[i] = std::sqrt(clamped);
        }
        equilibrium.azimuthal_speed.front() = equilibrium.azimuthal_speed[1];
        equilibrium.azimuthal_speed.back() = equilibrium.azimuthal_speed[radial_samples - 2];
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> azimuth(0.0, two_pi);
    std::normal_distribution<double> normal(0.0, 1.0);

    auto escape_speed = [&](double radius) {
        const double potential = lookup(equilibrium.radius, equilibrium.potential, radius);
        return std::sqrt(std::max(0.0, -2.0 * potential));
    };
    auto draw_bound = [&](double sigma, double limit, double& vx, double& vy, double& vz) {
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;
        if (sigma <= 0.0 || limit <= 0.0) {
            return;
        }
        for (int attempt = 0; attempt < 32; attempt++) {
            const double x = sigma * normal(rng);
            const double y = sigma * normal(rng);
            const double z = sigma * normal(rng);
            if (x * x + y * y + z * z < limit * limit) {
                vx = x;
                vy = y;
                vz = z;
                return;
            }
        }
    };

    const std::size_t total =
        profile.halo_count + profile.bulge_count + profile.disk_count;
    std::vector<Body> bodies;
    bodies.reserve(total);

    auto place_sphere = [&](std::size_t count, double mass, auto sample_radius, auto sigma_of, Config style) {
        const double particle_mass = mass / static_cast<double>(count);
        for (std::size_t n = 0; n < count; n++) {
            Body body(3, particle_mass, style);
            const double radius = sample_radius();
            const double mu = 2.0 * unit(rng) - 1.0;
            const double phi = azimuth(rng);
            const double planar = std::sqrt(std::max(0.0, 1.0 - mu * mu));
            body.r[0] = radius * planar * std::cos(phi);
            body.r[1] = radius * planar * std::sin(phi);
            body.r[2] = radius * mu;

            const double sigma = lookup(equilibrium.radius, sigma_of, radius);
            draw_bound(
                sigma,
                escape_speed(radius),
                body.v[0],
                body.v[1],
                body.v[2]
            );
            bodies.push_back(std::move(body));
        }
    };

    if (has_halo) {
        place_sphere(
            profile.halo_count,
            profile.halo_mass,
            [&]() {
                const double target = unit(rng);
                double lo = 0.0;
                double hi = concentration;
                for (int i = 0; i < 60; i++) {
                    const double mid = 0.5 * (lo + hi);
                    if (nfw_mass_factor(mid) / nfw_virial_factor < target) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return 0.5 * (lo + hi) * profile.halo_scale_radius;
            },
            equilibrium.halo_sigma,
            profile.halo_style
        );
    }

    if (has_bulge) {
        place_sphere(
            profile.bulge_count,
            profile.bulge_mass,
            [&]() {
                const double scaled = std::sqrt(unit(rng) * bulge_enclosed_fraction);
                return profile.bulge_scale_radius * scaled / (1.0 - scaled);
            },
            equilibrium.bulge_sigma,
            profile.bulge_style
        );
    }

    if (has_disk) {
        const double particle_mass =
            profile.disk_mass / static_cast<double>(profile.disk_count);
        const double z_limit = std::tanh(6.0);
        for (std::size_t n = 0; n < profile.disk_count; n++) {
            Body body(3, particle_mass, profile.disk_style);
            const double target = unit(rng);
            double lo = 0.0;
            double hi = disk_x_outer;
            for (int i = 0; i < 60; i++) {
                const double mid = 0.5 * (lo + hi);
                const double enclosed =
                    (1.0 - std::exp(-mid) * (1.0 + mid)) / disk_enclosed_fraction;
                if (enclosed < target) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            const double cylindrical = 0.5 * (lo + hi) * profile.disk_scale_radius;
            const double phi = azimuth(rng);
            const double vertical_unit = (2.0 * unit(rng) - 1.0) * z_limit;
            const double z = profile.disk_scale_height * std::atanh(vertical_unit);

            body.r[0] = cylindrical * std::cos(phi);
            body.r[1] = cylindrical * std::sin(phi);
            body.r[2] = z;

            const double radial_sigma =
                lookup(equilibrium.radius, equilibrium.radial_sigma, cylindrical);
            const double azimuthal_sigma =
                lookup(equilibrium.radius, equilibrium.azimuthal_sigma, cylindrical);
            const double vertical_sigma =
                lookup(equilibrium.radius, equilibrium.vertical_sigma, cylindrical);
            const double streaming =
                lookup(equilibrium.radius, equilibrium.azimuthal_speed, cylindrical);
            const double limit = escape_speed(cylindrical);

            double radial_velocity = 0.0;
            double azimuthal_velocity = streaming;
            double vertical_velocity = 0.0;
            for (int attempt = 0; attempt < 32; attempt++) {
                radial_velocity = radial_sigma * normal(rng);
                azimuthal_velocity = streaming + azimuthal_sigma * normal(rng);
                vertical_velocity = vertical_sigma * normal(rng);
                const double speed2 =
                    radial_velocity * radial_velocity
                    + azimuthal_velocity * azimuthal_velocity
                    + vertical_velocity * vertical_velocity;
                if (limit <= 0.0 || speed2 < limit * limit) {
                    break;
                }
                radial_velocity = 0.0;
                azimuthal_velocity = streaming;
                vertical_velocity = 0.0;
            }

            body.v[0] =
                radial_velocity * std::cos(phi) - azimuthal_velocity * std::sin(phi);
            body.v[1] =
                radial_velocity * std::sin(phi) + azimuthal_velocity * std::cos(phi);
            body.v[2] = vertical_velocity;
            bodies.push_back(std::move(body));
        }
    }

    double mass = 0.0;
    double com_x = 0.0;
    double com_y = 0.0;
    double com_z = 0.0;
    double com_vx = 0.0;
    double com_vy = 0.0;
    double com_vz = 0.0;
    for (const Body& body : bodies) {
        mass += body.m;
        com_x += body.m * body.r[0];
        com_y += body.m * body.r[1];
        com_z += body.m * body.r[2];
        com_vx += body.m * body.v[0];
        com_vy += body.m * body.v[1];
        com_vz += body.m * body.v[2];
    }
    com_x /= mass;
    com_y /= mass;
    com_z /= mass;
    com_vx /= mass;
    com_vy /= mass;
    com_vz /= mass;

    for (Body& body : bodies) {
        body.r[0] -= com_x;
        body.r[1] -= com_y;
        body.r[2] -= com_z;
        body.v[0] -= com_vx;
        body.v[1] -= com_vy;
        body.v[2] -= com_vz;
        rotate_z(body.r, profile.position_angle);
        rotate_z(body.v, profile.position_angle);
        rotate_x(body.r, profile.inclination);
        rotate_x(body.v, profile.inclination);
        body.r[0] += center_x_;
        body.r[1] += center_y_;
        body.r[2] += center_z_;
    }

    return State(std::move(bodies));
}
