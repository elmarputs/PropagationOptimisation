#ifndef TUDAT_APPLICATION_PAGMO_PROBLEM_SHAPE_OPTIMIZATION_H
#define TUDAT_APPLICATION_PAGMO_PROBLEM_SHAPE_OPTIMIZATION_H

#include <tudat/SimulationSetup/tudatSimulationHeader.h>
#include <vector>

typedef std::vector<double> vector_double;

namespace tudat
{

    class ShapeOptimization
    {
    public:

        // Empty constructor
        ShapeOptimization();

        vector_double fitness(const vector_double &cv) const;

        std::pair<vector_double, vector_double> get_bounds() const;

    };

}

#endif