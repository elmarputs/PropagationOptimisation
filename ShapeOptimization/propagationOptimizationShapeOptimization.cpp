/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <pagmo/problem.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/simulated_annealing.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include <tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h>

#include "../applicationOutput.h"
#include "shapeOptimization.h"

#include <iostream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::aerodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::mathematical_constants;
using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::statistics;
using namespace pagmo;


//! Execute propagation of orbits of Capsule during entry.
int main( )
{

    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimisationGroup" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    problem prob{ShapeOptimization( ) };


    // Instantiate a pagmo algorithm
    algorithm algo{de1220( )};

    pagmo::population::size_type populationSize = 128;

    island isl{algo, prob, populationSize};

    // Evolve for 25 generations
    for( int i = 0; i < 25; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( isl.get_population( ).get_x( ), "targetingPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
        printPopulationToFile( isl.get_population( ).get_f( ), "targetingPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

        std::cout<<i<<std::endl;
    }

	// The exit code EXIT_SUCCESS indicates that the program was successfully executed.
	return EXIT_SUCCESS;
}
