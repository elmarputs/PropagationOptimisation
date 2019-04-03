#include "shapeOptimization.h"

#include <cmath>

tudat::ShapeOptimization::ShapeOptimization() = default;

vector_double tudat::ShapeOptimization::fitness(const std::vector<double> &cv) const
{

	// Himmelblau test function
	double x = cv.at(0);
	double y = cv.at(1);

	// Compute fitness (Himmelblau)
	//double fitnessValue = std::pow(x*x + y - 11.0, 2.0) + std::pow(x + y*y - 7.0, 2.0);

	// Compute fitness (Rosenbrock)
	//double fitnessValue = 100.0*std::pow(y - x*x, 2.0) + std::pow(1 - x*x, 2.0);

	// Compute fitness (Booth)
	double fitnessValue = std::pow(x + 2.0*y - 7, 2.0) + std::pow(2.0*x + y - 5, 2.0);

    vector_double fitness;
    fitness.push_back(fitnessValue);
    return fitness;
}

std::pair<vector_double, vector_double> tudat::ShapeOptimization::get_bounds() const
{
    std::pair<vector_double, vector_double> bounds;
    bounds.first = {-5.0, -5.0};
    bounds.second = {5.0, 5.0};

    return bounds;
}