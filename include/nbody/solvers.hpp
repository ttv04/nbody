#pragma once

#include "state.hpp"
#include <memory>
#include <vector>

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

class BHNode {
    friend class BHIndexedOrthoTree;
public:
    BHNode(const State& state, std::vector<double> lo, std::vector<double> hi);
    BHNode(const BHNode&) = delete;
    BHNode& operator=(const BHNode&) = delete;
    ~BHNode();

    void accelerationSum(
        size_t particle_index,
        double theta,
        double softening,
        double grav,
        AccelerationField& acceleration
    );
private:
    const State* state_;
    bool is_leaf_;
    double mass_;
    double width_;
    std::vector<double> lo_;
    std::vector<double> hi_;
    std::vector<double> center_;
    std::vector<size_t> particles_;
    std::vector<BHNode*> children_;
};

struct BHBound {
    std::vector<double> lo;
    std::vector<double> hi;
};

class BHIndexedOrthoTree {
public:
    BHIndexedOrthoTree(const State& state, size_t max_particles_per_leaf);
    BHIndexedOrthoTree(const BHIndexedOrthoTree&) = delete;
    BHIndexedOrthoTree& operator=(const BHIndexedOrthoTree&) = delete;
    ~BHIndexedOrthoTree();

    void insert(size_t index);
    void calculateAcceleration(
        size_t particle_index,
        double theta,
        double softening,
        double grav,
        AccelerationField& acceleration
    );
    void collect_bounds(std::vector<BHBound>& out) const;
    void print() const;
private:
    void insertHelper(size_t index, BHNode* node);
    void collectBoundsHelper(const BHNode* node, std::vector<BHBound>& out) const;
    void printHelper(const BHNode* node, int depth) const;
    void calculateCenterOfMass(BHNode* node);

    const State& state_;
    size_t max_particles_per_leaf_;
    BHNode* root_;
};

class BarnesHutSolver : public Solver {
public:
    BarnesHutSolver(
        double softening,
        double gravitational_constant = 6.67430e-11,
        double theta = 0.75,
        size_t max_particles_per_leaf = 12
    );
    
    void solve(
        const State& state,
        AccelerationField& acceleration
    ) const;

    const std::vector<BHBound>& last_bounds() const;

private:
    double theta_;
    size_t max_particles_per_leaf_;
    mutable std::vector<BHBound> last_bounds_;
};
