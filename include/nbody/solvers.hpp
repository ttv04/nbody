#pragma once

#include "state.hpp"
#include <memory>

void accumulate_acceleration(
    AccelerationField& acceleration,
    size_t target,
    const std::vector<double>& target_pos,
    const std::vector<double>& source_pos,
    double source_mass,
    size_t dim,
    double softening,
    double grav
);

class Solver {
public:
    Solver(double softening = 0.0, double gravitational_constant = 6.67430e-11);
    virtual ~Solver() = default;

    virtual void solve(
        const State& state,
        AccelerationField& acceleration
    ) const = 0;

protected:
    double softening_;
    double gravitational_constant_;
};

class NaiveSolver : public Solver {
public:
    NaiveSolver(
        double softening,
        double gravitational_constant = 6.67430e-11
    );

    void solve(
        const State& state,
        AccelerationField& acceleration
    ) const;
};

class BHBound {
    friend class BHInternalNode;
    friend class BHIndexedOrthoTree;
public:
    BHBound() = default;
    BHBound(std::vector<std::vector<double>> bounds);

    double l(size_t dimension) const;
    double r(size_t dimension) const;
    double width() const;
private:
    std::vector<std::vector<double>> bounds;
    double width_ = 0.0;
};

class BHNode {
    friend class BHIndexedOrthoTree;
public:
    virtual bool is_leaf() const = 0;
    virtual void accelerationSum(
        size_t particle_index,
        double theta,
        double softening,
        double grav,
        AccelerationField& acceleration
    ) = 0;
    virtual ~BHNode() = default;
protected:
    BHNode(const BHBound& bound)
        : mass_(0.0),
          bound_(bound) {}

    std::vector<double> center_;
    double mass_;
    BHBound bound_;
};

class BHLeafNode : public BHNode {
    friend class BHIndexedOrthoTree;
public:
    BHLeafNode(const State& state, const BHBound& bound);
    bool is_leaf() const;
    void accelerationSum(
        size_t particle_index,
        double theta,
        double softening,
        double grav,
        AccelerationField& acceleration
    );
private:
    std::vector<size_t> particle_indexes_;
    const State& state_;
};

class BHInternalNode : public BHNode {
    friend class BHIndexedOrthoTree;
public:
    BHInternalNode(const State& state, const BHBound& bound);
    bool is_leaf() const;
    void accelerationSum(
        size_t particle_index,
        double theta,
        double softening,
        double grav,
        AccelerationField& acceleration
    );
private:
    std::vector<std::unique_ptr<BHNode>> children_;
    const State& state_;
};

class BHIndexedOrthoTree {
public:
    BHIndexedOrthoTree(const State& state, size_t max_particles_per_leaf);
    void insert(size_t index);
    void calculateAcceleration(
        size_t particle_index,
        double theta,
        double softening,
        double grav,
        AccelerationField& acceleration
    );
    void print() const;
private:
    void insertHelper(size_t index, std::unique_ptr<BHNode>& node);
    void printHelper(const BHNode* node, int depth) const;
    void calculateCenterOfMass(std::unique_ptr<BHNode>& node);
    const State& state_;
    size_t max_particles_per_leaf_;
    std::unique_ptr<BHNode> root_;
};

class BarnesHutSolver : public Solver {
public:
    BarnesHutSolver(
        double softening,
        double gravitational_constant = 6.67430e-11
    );
    
    void solve(
        const State& state,
        AccelerationField& acceleration
    ) const;
};
