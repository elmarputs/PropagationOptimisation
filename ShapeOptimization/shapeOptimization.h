#ifndef TUDAT_APPLICATION_PAGMO_PROBLEM_SHAPE_OPTIMIZATION_H
#define TUDAT_APPLICATION_PAGMO_PROBLEM_SHAPE_OPTIMIZATION_H

#include <tudat/SimulationSetup/tudatSimulationHeader.h>
#include <vector>
#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <pagmo/problem.hpp>

#include "../applicationOutput.h"

#include <iostream>

#include <boost/filesystem.hpp>

#include "applicationOutput.h"

#include "Tudat/InputOutput/basicInputOutput.h"


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


    private:
        simulation_setup::NamedBodyMap bodyMap_;
        double simulationStartEpoch;
        double simulationEndEpoch;
        double vehicleDensity;
        Eigen::Vector6d capsuleSphericalEntryState;
    };

}

/*
void printPopulationToFile( const std::vector< std::vector< double > >& population,
                            const std::string fileSuffix,
                            const bool isFitness )
{

    Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
    for( unsigned int i = 0; i < population.size( ); i++ )
    {
        for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
        {
            matrixToPrint( i, j ) = population.at( i ).at( j );
        }
    }

    if( !isFitness )
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint, "population_" + fileSuffix + ".dat", 16,
                                               tudat_applications::getOutputPath( ) );
    }
    else
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint, "fitness_" + fileSuffix + ".dat", 16,
                                                tudat_applications::getOutputPath( )  );
    }
}
*/
#endif
