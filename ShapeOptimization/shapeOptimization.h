#ifndef TUDAT_APPLICATION_PAGMO_PROBLEM_SHAPE_OPTIMIZATION_H
#define TUDAT_APPLICATION_PAGMO_PROBLEM_SHAPE_OPTIMIZATION_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <pagmo/problem.hpp>

#include "../applicationOutput.h"

#include <boost/filesystem.hpp>


typedef std::vector<double> vector_double;

namespace tudat
{

    class ShapeOptimization
    {
    public:

    	// Static variables
	    static simulation_setup::NamedBodyMap bodyMap_;
	    static Eigen::Vector6d capsuleSphericalEntryState;
	    static double simulationStartEpoch;
	    static double simulationEndEpoch;
	    static bool isSpiceLoaded;

	    double heatLoadFunction( const double, const double, tudat::simulation_setup::NamedBodyMap& ) const;
        double flightRangeFunction( const double, const double, tudat::simulation_setup::NamedBodyMap& ) const;
	    // Empty constructor
        ShapeOptimization();

        vector_double::size_type get_nobj() const
        {
            return 2;
        }

        vector_double::size_type get_nec() const
		{
        	return 0;
		}

		vector_double::size_type get_nic() const
		{
        	return 0;
		}

        vector_double fitness( const vector_double& ) const;

        std::pair<vector_double, vector_double> get_bounds() const;

        template <typename Archive>
        void serialize(Archive &ar)
        {
        	//ar(simulationStartEpoch, simulationEndEpoch);
        }


    private:

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

PAGMO_REGISTER_PROBLEM(tudat::ShapeOptimization)

#endif
