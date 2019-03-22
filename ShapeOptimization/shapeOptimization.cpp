#include "shapeOptimization.h"

tudat::ShapeOptimization::ShapeOptimization() = default;

vector_double tudat::ShapeOptimization::fitness(const std::vector<double> &x) const
{
    vector_double fitness;
    fitness.push_back(0.0);
    return fitness;
}

std::pair<vector_double, vector_double> tudat::ShapeOptimization::get_bounds() const
{
    std::pair<vector_double, vector_double> bounds;
    bounds.first = {0.0, 0.0};
    bounds.second = {1.0, 1.0};

    return bounds;
}