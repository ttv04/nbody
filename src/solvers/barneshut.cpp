#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include <cmath>
#include <iostream>
#include <string>

BHNode::BHNode(const State& state, std::vector<double> lo, std::vector<double> hi)
    : state_(&state),
      is_leaf_(true),
      mass_(0.0),
      width_(0.0),
      lo_(std::move(lo)),
      hi_(std::move(hi)),
      center_(state.dim, 0.0) {
    for (size_t d = 0; d < lo_.size(); d++) {
        const double side = hi_[d] - lo_[d];
        if (side > width_) {
            width_ = side;
        }
    }
}

BHNode::~BHNode() {
    for (BHNode* child : children_) {
        delete child;
    }
}

void BHNode::accelerationSum(
    size_t particle_index,
    double theta,
    double softening,
    double grav,
    AccelerationField& acceleration
) {
    const Body& body = state_->bodies[particle_index];

    if (is_leaf_) {
        for (size_t j : particles_) {
            if (j == particle_index) {
                continue;
            }
            const Body& other = state_->bodies[j];
            accumulate_acceleration(
                acceleration,
                particle_index,
                body.r,
                other.r,
                other.m,
                state_->dim,
                softening,
                grav
            );
        }
        return;
    }

    double r2 = 0.0;
    for (size_t d = 0; d < state_->dim; d++) {
        const double displacement = center_[d] - body.r[d];
        r2 += displacement * displacement;
    }

    if (r2 > 0.0 && width_ * width_ < theta * theta * r2) {
        accumulate_acceleration(
            acceleration,
            particle_index,
            body.r,
            center_,
            mass_,
            state_->dim,
            softening,
            grav
        );
        return;
    }

    for (BHNode* child : children_) {
        if (child == nullptr) {
            continue;
        }
        child->accelerationSum(
            particle_index,
            theta,
            softening,
            grav,
            acceleration
        );
    }
}

void BHIndexedOrthoTree::insertHelper(size_t index, BHNode* node) {
    if (node->is_leaf_) {
        if (node->particles_.size() < max_particles_per_leaf_) {
            node->particles_.push_back(index);
            return;
        }

        std::vector<size_t> old = node->particles_;
        node->particles_.clear();
        node->is_leaf_ = false;
        node->children_.assign(size_t{1} << state_.dim, nullptr);

        for (size_t ind : old) {
            insertHelper(ind, node);
        }
        insertHelper(index, node);
        return;
    }

    size_t mask = 0;
    for (size_t d = 0; d < state_.dim; d++) {
        const double mid = 0.5 * (node->lo_[d] + node->hi_[d]);
        if (state_.bodies[index].r[d] >= mid) {
            mask |= (size_t{1} << d);
        }
    }

    if (node->children_[mask] == nullptr) {
        std::vector<double> child_lo(state_.dim);
        std::vector<double> child_hi(state_.dim);
        for (size_t d = 0; d < state_.dim; d++) {
            const double lo = node->lo_[d];
            const double hi = node->hi_[d];
            const double mid = 0.5 * (lo + hi);
            if (mask & (size_t{1} << d)) {
                child_lo[d] = mid;
                child_hi[d] = hi;
            } else {
                child_lo[d] = lo;
                child_hi[d] = mid;
            }
        }

        BHNode* child = new BHNode(state_, std::move(child_lo), std::move(child_hi));
        node->children_[mask] = child;
    }

    insertHelper(index, node->children_[mask]);
}

void BHIndexedOrthoTree::insert(size_t index) {
    insertHelper(index, root_);
}

void BHIndexedOrthoTree::collect_bounds(std::vector<BHBound>& out) const {
    out.clear();
    collectBoundsHelper(root_, out);
}

void BHIndexedOrthoTree::collectBoundsHelper(
    const BHNode* node,
    std::vector<BHBound>& out
) const {
    if (node == nullptr) {
        return;
    }

    out.push_back(BHBound{node->lo_, node->hi_});
    if (node->is_leaf_) {
        return;
    }

    for (const BHNode* child : node->children_) {
        collectBoundsHelper(child, out);
    }
}

void BHIndexedOrthoTree::print() const {
    printHelper(root_, 0);
}

void BHIndexedOrthoTree::printHelper(const BHNode* node, int depth) const {
    const std::string indent(static_cast<size_t>(depth) * 2, ' ');

    if (node == nullptr) {
        std::cout << indent << "(empty)\n";
        return;
    }

    std::cout << indent;
    if (node->is_leaf_) {
        std::cout << "Leaf";
    } else {
        std::cout << "Internal";
    }

    std::cout << " bound=";
    for (size_t d = 0; d < state_.dim; d++) {
        if (d > 0) {
            std::cout << " x ";
        }
        std::cout << "[" << node->lo_[d] << ", " << node->hi_[d] << "]";
    }

    if (node->is_leaf_) {
        std::cout << " particles=[";
        for (size_t i = 0; i < node->particles_.size(); i++) {
            if (i > 0) {
                std::cout << ", ";
            }
            const size_t ind = node->particles_[i];
            std::cout << ind << " (";
            for (size_t d = 0; d < state_.dim; d++) {
                if (d > 0) {
                    std::cout << ", ";
                }
                std::cout << state_.bodies[ind].r[d];
            }
            std::cout << ")";
        }
        std::cout << "]\n";
        return;
    }

    std::cout << "\n";
    for (size_t i = 0; i < node->children_.size(); i++) {
        std::cout << indent << "  child " << i << ":\n";
        printHelper(node->children_[i], depth + 2);
    }
}

void BHIndexedOrthoTree::calculateCenterOfMass(BHNode* node) {
    if (node == nullptr) {
        return;
    }

    node->center_.assign(state_.dim, 0.0);

    if (node->is_leaf_) {
        double total_mass = 0.0;
        for (size_t index : node->particles_) {
            total_mass += state_.bodies[index].m;
        }
        node->mass_ = total_mass;
        if (total_mass == 0.0) {
            return;
        }
        for (size_t d = 0; d < state_.dim; d++) {
            double total_mass_distance = 0.0;
            for (size_t index : node->particles_) {
                total_mass_distance +=
                    state_.bodies[index].m * state_.bodies[index].r[d];
            }
            node->center_[d] = total_mass_distance / total_mass;
        }
        return;
    }

    double total_mass = 0.0;
    for (BHNode* child : node->children_) {
        if (child == nullptr) {
            continue;
        }
        calculateCenterOfMass(child);
        total_mass += child->mass_;
    }
    node->mass_ = total_mass;
    if (total_mass == 0.0) {
        return;
    }
    for (size_t d = 0; d < state_.dim; d++) {
        double total_mass_distance = 0.0;
        for (BHNode* child : node->children_) {
            if (child == nullptr) {
                continue;
            }
            total_mass_distance += child->center_[d] * child->mass_;
        }
        node->center_[d] = total_mass_distance / total_mass;
    }
}

BHIndexedOrthoTree::BHIndexedOrthoTree(
    const State& state,
    size_t max_particles_per_leaf
)
    : state_(state),
      max_particles_per_leaf_(max_particles_per_leaf),
      root_(nullptr) {

    std::vector<double> lo(state_.dim);
    std::vector<double> hi(state_.dim);
    for (size_t d = 0; d < state_.bodies[0].dimensions(); d++) {
        lo[d] = state_.bodies[0].r[d];
        hi[d] = lo[d];
    }
    for (const Body& body : state_.bodies) {
        for (size_t d = 0; d < body.dimensions(); d++) {
            lo[d] = fmin(lo[d], body.r[d]);
            hi[d] = fmax(hi[d], body.r[d]);
        }
    }

    double box_min = lo[0];
    double box_max = hi[0];
    for (size_t d = 1; d < state_.dim; d++) {
        box_min = fmin(box_min, lo[d]);
        box_max = fmax(box_max, hi[d]);
    }
    for (size_t d = 0; d < state_.dim; d++) {
        lo[d] = box_min;
        hi[d] = box_max;
    }

    root_ = new BHNode(state_, std::move(lo), std::move(hi));

    for (size_t i = 0; i < state.size(); i++) {
        insert(i);
    }

    calculateCenterOfMass(root_);
}

BHIndexedOrthoTree::~BHIndexedOrthoTree() {
    delete root_;
}

void BHIndexedOrthoTree::calculateAcceleration(
    size_t particle_index,
    double theta,
    double softening,
    double grav,
    AccelerationField& acceleration
) {
    root_->accelerationSum(particle_index, theta, softening, grav, acceleration);
}

BarnesHutSolver::BarnesHutSolver(
    double softening,
    double gravitational_constant,
    double theta,
    size_t max_particles_per_leaf
)
    : Solver(softening, gravitational_constant),
      theta_(theta),
      max_particles_per_leaf_(max_particles_per_leaf < 1 ? 1 : max_particles_per_leaf) {}

void BarnesHutSolver::solve(
    const State& state,
    AccelerationField& acceleration
) const {
    if (
        acceleration.body_count() != state.size()
        || acceleration.dimensions() != state.dim
    ) {
        acceleration.resize(state.size(), state.dim);
    } else {
        acceleration.reset();
    }

    BHIndexedOrthoTree tree(state, max_particles_per_leaf_);
    tree.collect_bounds(last_bounds_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < static_cast<int>(state.size()); i++) {
        tree.calculateAcceleration(
            static_cast<size_t>(i),
            theta_,
            softening_,
            gravitational_constant_,
            acceleration
        );
    }
}

const std::vector<BHBound>& BarnesHutSolver::last_bounds() const {
    return last_bounds_;
}
