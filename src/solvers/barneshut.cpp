#include "nbody/solvers.hpp"
#include "nbody/state.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#define THETA 0.75
#define MAX_PARTICLES_PER_LEAF 12

BHBound::BHBound(std::vector<std::vector<double>> bounds)
    : bounds(bounds),
      width_(0.0) {
    for (size_t d = 0; d < this->bounds.size(); d++) {
        const double side = this->bounds[d][1] - this->bounds[d][0];
        if (side > width_) {
            width_ = side;
        }
    }
}

double BHBound::l(size_t dimension) const {
    return bounds[dimension][0];
}

double BHBound::r(size_t dimension) const {
    return bounds[dimension][1];
}

double BHBound::width() const {
    return width_;
}

BHLeafNode::BHLeafNode(
    const State& state,
    const BHBound& bound
)
    : BHNode(bound),
      state_(state) {}

bool BHLeafNode::is_leaf() const {
    return true;
}

void BHLeafNode::accelerationSum(
    size_t particle_index,
    double theta,
    double softening,
    double grav,
    AccelerationField& acceleration
) {
    (void)theta;

    const Body& body = state_.bodies[particle_index];

    for (size_t j : particle_indexes_) {
        if (j == particle_index) {
            continue;
        }

        const Body& other = state_.bodies[j];
        accumulate_acceleration(
            acceleration,
            particle_index,
            body.r,
            other.r,
            other.m,
            state_.dim,
            softening,
            grav
        );
    }
}

BHInternalNode::BHInternalNode(
    const State& state,
    const BHBound& bound
)
    : BHNode(bound),
      state_(state) {
    children_.resize(size_t{1} << state_.dim);
}

bool BHInternalNode::is_leaf() const {
    return false;
}

void BHInternalNode::accelerationSum(
    size_t particle_index,
    double theta,
    double softening,
    double grav,
    AccelerationField& acceleration
) {
    const Body& body = state_.bodies[particle_index];

    double r2 = 0.0;
    if (center_.size() == state_.dim) {
        for (size_t d = 0; d < state_.dim; d++) {
            const double displacement = center_[d] - body.r[d];
            r2 += displacement * displacement;
        }
    }

    const double width = bound_.width();

    if (r2 > 0.0 && width * width < theta * theta * r2) {
        accumulate_acceleration(
            acceleration,
            particle_index,
            body.r,
            center_,
            mass_,
            state_.dim,
            softening,
            grav
        );
        return;
    }

    for (size_t i = 0; i < children_.size(); i++) {
        if (children_[i] == nullptr) {
            continue;
        }
        children_[i]->accelerationSum(
            particle_index,
            theta,
            softening,
            grav,
            acceleration
        );
    }
}

void BHIndexedOrthoTree::insertHelper(
    size_t index, 
    std::unique_ptr<BHNode>& node) {

    if (node->is_leaf()) {
        BHLeafNode* leaf = static_cast<BHLeafNode*>(node.get());
        if (leaf->particle_indexes_.size() < max_particles_per_leaf_) {
            leaf->particle_indexes_.push_back(index);
        }
        else {
            std::unique_ptr<BHNode> newNode =
                std::make_unique<BHInternalNode>(state_, leaf->bound_);

            for (size_t ind : leaf->particle_indexes_) {
                insertHelper(ind, newNode);
            }
            insertHelper(index, newNode);
            node = std::move(newNode);
        }
    }
    else {
        BHInternalNode* internal = static_cast<BHInternalNode*>(node.get());

        std::vector<bool> bits(state_.dim);
        for (size_t d = 0; d < state_.dim; d++) {
            const double mid =
                0.5 * (internal->bound_.l(d) + internal->bound_.r(d));
            bits[d] = state_.bodies[index].r[d] >= mid;
        }

        size_t tree_index = 0;
        for (size_t i = 0; i < state_.dim; i++) {
            tree_index += static_cast<size_t>(bits[i]) * (size_t{1} << i);
        }

        if (internal->children_[tree_index] == nullptr) {
            std::vector<std::vector<double>> child_b(
                state_.dim,
                std::vector<double>(2)
            );
            for (size_t d = 0; d < state_.dim; d++) {
                const double lo = internal->bound_.l(d);
                const double hi = internal->bound_.r(d);
                const double mid = 0.5 * (lo + hi);
                if (bits[d]) {
                    child_b[d][0] = mid;
                    child_b[d][1] = hi;
                } else {
                    child_b[d][0] = lo;
                    child_b[d][1] = mid;
                }
            }

            internal->children_[tree_index] =
                std::make_unique<BHLeafNode>(state_, BHBound(child_b));
        }

        insertHelper(index, internal->children_[tree_index]);
    }
};

void BHIndexedOrthoTree::insert(size_t index) {
    insertHelper(index, root_);
};

void BHIndexedOrthoTree::print() const {
    printHelper(root_.get(), 0);
}

void BHIndexedOrthoTree::printHelper(const BHNode* node, int depth) const {
    const std::string indent(static_cast<size_t>(depth) * 2, ' ');

    if (node == nullptr) {
        std::cout << indent << "(empty)\n";
        return;
    }

    std::cout << indent;
    if (node->is_leaf()) {
        std::cout << "Leaf";
    } else {
        std::cout << "Internal";
    }

    std::cout << " bound=";
    for (size_t d = 0; d < state_.dim; d++) {
        if (d > 0) {
            std::cout << " x ";
        }
        std::cout << "[" << node->bound_.l(d) << ", " << node->bound_.r(d) << "]";
    }

    if (node->is_leaf()) {
        const BHLeafNode* leaf = static_cast<const BHLeafNode*>(node);
        std::cout << " particles=[";
        for (size_t i = 0; i < leaf->particle_indexes_.size(); i++) {
            if (i > 0) {
                std::cout << ", ";
            }
            const size_t ind = leaf->particle_indexes_[i];
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
    const BHInternalNode* internal = static_cast<const BHInternalNode*>(node);
    for (size_t i = 0; i < internal->children_.size(); i++) {
        std::cout << indent << "  child " << i << ":\n";
        printHelper(internal->children_[i].get(), depth + 2);
    }
}

void BHIndexedOrthoTree::calculateCenterOfMass(std::unique_ptr<BHNode>& node) {
    if (node == nullptr) {
        return;
    }

    node->center_.assign(state_.dim, 0.0);

    if (node->is_leaf()) {
        BHLeafNode* leaf = static_cast<BHLeafNode*>(node.get());
        double total_mass = 0.0;
        for (size_t index : leaf->particle_indexes_) {
            const Body& body = state_.bodies[index];
            total_mass += body.m;
        }
        leaf->mass_ = total_mass;
        if (total_mass == 0.0) {
            return;
        }
        for (size_t d = 0; d < state_.dim; d++) {
            double total_mass_distance = 0.0;
            for (size_t index : leaf->particle_indexes_) {
                const Body& body = state_.bodies[index];
                total_mass_distance += body.m * body.r[d];
            }
            leaf->center_[d] = total_mass_distance / total_mass;
        }
    }
    else {
        BHInternalNode* internal = static_cast<BHInternalNode*>(node.get());
        double total_mass = 0.0;
        for (size_t i = 0; i < internal->children_.size(); i++) {
            if (internal->children_[i] == nullptr) {
                continue;
            }
            calculateCenterOfMass(internal->children_[i]);
            total_mass += internal->children_[i]->mass_;
        }
        internal->mass_ = total_mass;
        if (total_mass == 0.0) {
            return;
        }
        for (size_t d = 0; d < state_.dim; d++) {
            double total_mass_distance = 0.0;
            for (size_t i = 0; i < internal->children_.size(); i++) {
                if (internal->children_[i] == nullptr) {
                    continue;
                }
                total_mass_distance +=
                    internal->children_[i]->center_[d]
                    * internal->children_[i]->mass_;
            }
            internal->center_[d] = total_mass_distance / total_mass;
        }
    }
};

BHIndexedOrthoTree::BHIndexedOrthoTree(
    const State& state,
    size_t max_particles_per_leaf
)
    : state_(state),
      max_particles_per_leaf_(max_particles_per_leaf) {

    std::vector<std::vector<double>> b(state_.dim, std::vector<double>(2));
    for (size_t d = 0; d < state_.bodies[0].dimensions(); d++) {
        b[d][0] = state_.bodies[0].r[d];
        b[d][1] = b[d][0];
    }

    for (const Body& body : state_.bodies) {
        for (size_t d = 0; d < body.dimensions(); d++) {
            b[d][0] = fmin(b[d][0], body.r[d]);
            b[d][1] = fmax(b[d][1], body.r[d]);
        }
    }

    double box_min = b[0][0];
    double box_max = b[0][1];
    for (size_t d = 1; d < state_.dim; d++) {
        box_min = fmin(box_min, b[d][0]);
        box_max = fmax(box_max, b[d][1]);
    }
    for (size_t d = 0; d < state_.dim; d++) {
        b[d][0] = box_min;
        b[d][1] = box_max;
    }

    root_ = std::make_unique<BHLeafNode>(state, BHBound(b));
    
    for (size_t i = 0; i < state.size(); i++) {
        insert(i);
    } 

    calculateCenterOfMass(root_);
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
    double gravitational_constant
)
    : Solver(softening, gravitational_constant) {}

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

    BHIndexedOrthoTree tree(state, MAX_PARTICLES_PER_LEAF);
// #ifdef _OPENMP
#pragma omp parallel for schedule(static)
// #endif
    for (size_t i = 0; i < state.size(); i++) {
        tree.calculateAcceleration(
            i,
            THETA,
            softening_,
            gravitational_constant_,
            acceleration
        );
    }
}
