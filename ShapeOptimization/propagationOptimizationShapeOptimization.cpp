/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "shapeOptimization.h"

#include <pagmo/algorithms/sade.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/moead.hpp>
#include <pagmo/algorithms/nsga2.hpp>
#include <pagmo/algorithms/ihs.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/simulated_annealing.hpp>
#include "getAlgorithm.h"
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include "saveOptimizationResults.h"


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
using namespace tudat_pagmo_applications;


double getVariableSettingFromCodedValue(const double codedValue, const double min, const double max)
{
	// Find midpoint and spread
	double midpoint = (min + max) / 2.0;
	double diff		= max - midpoint;

	return midpoint + diff * codedValue;
}

void performCompositeDesign(const pagmo::problem &prob, const std::string &fileName, const std::string &filePath)
{
	// Load CCD matrix from file
	Eigen::MatrixXd ccdMatrix = tudat::input_output::readMatrixFromFile(filePath + fileName, ",");
	std::cout << "Successfully loaded coded variables matrix\n";

	// Convert CCD matrix with coded values to actual control variable values
	for( int i = 0; i < ccdMatrix.rows(); i++ )
	{
		for( int j = 0; j < ccdMatrix.cols(); j++ )
		{
			ccdMatrix(i, j) = getVariableSettingFromCodedValue(ccdMatrix(i, j), prob.get_bounds().first.at(j), prob.get_bounds().second.at(j));
		}
	}
	std::cout << "Converted coded variable values into real values\n";

	std::map<int, Eigen::VectorXd> ccdFitnessMap;

	// Loop through all designs
	for( int design = 0; design < ccdMatrix.rows(); design++ )
	{
		std::cout << "Running design " << design << "...\n";

		Eigen::VectorXd cv = ccdMatrix.row(design);
		std::cout << cv << std::endl;
		// Convert Eigen vector to std::vector for use with fitness function
		const std::vector<double> controlVars{&cv[0], cv.data() + cv.size()};

		// Run propagation with control variables for current design iteration and save to map
		std::vector<double> fitness = prob.fitness(controlVars);
		Eigen::VectorXd objectives(2);
		objectives << fitness.at(0), fitness.at(1);
		std::cout << "Fitness: " << objectives << std::endl;
		Eigen::VectorXd outputVector(cv.size() + objectives.size());
		outputVector << cv, objectives;
		ccdFitnessMap[design] = outputVector;
	}

	// Print map with fitness values to file
	tudat::input_output::writeDataMapToTextFile(ccdFitnessMap, "ccdFitness.txt", filePath);

}

int main( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed(255);

	bool runOptimisation;

	std::cout << "Please enter a 1 if you want to run the optimisation; for CCD enter a 0.\n";
    runOptimisation = true;

    //std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimisationGroup" );
    std::string outputPath{__FILE__};
    outputPath = outputPath.substr(0, outputPath.find_last_of("/\\") + 1);



    problem prob{ShapeOptimization( ) };
    std::cout<<"Created Problem \n";

	// Perform central composite design
	if(!runOptimisation)
	{
		performCompositeDesign(prob, "variableSettings.txt", outputPath);
	}


	if(runOptimisation)
    {
        unsigned int numberOfOptimizers = 1;
        unsigned int cases = 1;

        for( unsigned int currentQuestion = 1; currentQuestion <= numberOfOptimizers; currentQuestion++ )
        {
            if( currentQuestion == 1 )
            {
                // Define number of cases for the single optimizer question here
                cases = 1;
            }
            else
            {
                // Define cases for multi opmizer question here
                cases = 1;
            }

            for( unsigned int currentCase = 1; currentCase <= cases; currentCase++ )
            {
            std::cout<<"Initializing case:"<<currentCase<<"\n";

            if( currentQuestion != 1 )
            {
            currentCase = 2;
            }
            // Instantiate a pagmo algorithm
            // Create optimization algorithm

            // Start moead optimizer
            algorithm algo = getMultiObjectiveAlgorithm( currentQuestion-1 );
            algo.set_seed(255);
            //algorithm algo{moead( )};


            std::cout << "Created pagmo algorithm \n";

            pagmo::population::size_type populationSize = 32;
            if( currentCase == 1 )
            {
                populationSize = 21;
            }
            if( currentCase == 2 )
            {
                populationSize = 32;
            }
            else if( currentCase == 3 )
            {
                populationSize = 64;
            }
            else if( currentCase == 4 )
            {
                populationSize = 128;
            }


            std::cout << "Created populationSize \n";



            std::cout << "Created archipelago \n";

            int generations = 25;
            if( currentCase == 5 )
            {
                generations = 10;
            }
            if( currentCase == 6 )
            {
                generations = 25;
            }
            if( currentCase == 7)
            {
                generations = 50;
            }
            if( currentCase == 8 )
            {
                generations = 100;
            }
            /*
            //archipelago arch(4, algo, prob, populationSize);
            island arch(algo, prob, populationSize);

            // Evolve for 25 generations
            for (int i = 0; i < generations; i++)
            {
                std::cout << "Iteration " << i << " started \n";
                arch.evolve();
                while (arch.status() != pagmo::evolve_status::idle &&
                       arch.status() != pagmo::evolve_status::idle_error)
                {
                    arch.wait();
                }
                std::cout << "Waiting for islands to finish evolution...\n";
                try
                {
                    arch.wait_check(); // Raises errors
                }
                catch(std::exception& e)
                {
                    std::cout << "Something went wrong during evolution. Quitting application...\n";
                    std::cout << "Exception: " << e.what() << std::endl;
                    return 1;
                }
                std::cout << "Writing champions to file...\n";
                // Write current iteration results to file
                //std::cout << arch.get_champions_f()[0][0] << std::endl;
                //printPopulationToFile(arch.get_champions_f(), "targetingPropagation_" + std::to_string(i) + "_" + std::to_string(i), false);
                printPopulationToFile(arch.get_population().get_f(), "targetingPropagation_" + std::to_string(currentQuestion ) + "_" + std::to_string(currentCase) + "_" + std::to_string(i), true);
            }
            */
            // Create an island with 1024 individuals
            island isl{algo, prob, populationSize};

            // Evolve for 100 generations
            for( int i = 0 ; i < generations; i++ )
            {
                isl.evolve( );
                while( isl.status( ) != pagmo::evolve_status::idle &&
                       isl.status( ) != pagmo::evolve_status::idle_error )
                {
                    isl.wait( );
                }
                isl.wait_check( ); // Raises errors

                // Write current iteration results to file
                //printPopulationToFile( isl.get_population( ).get_x( ), "earthMarsLambert_" + std::to_string( j ) + "_" + std::to_string( i ) , false );
               // printPopulationToFile( isl.get_population( ).get_f( ), "earthMarsLambert_" + std::to_string( j ) + "_" + std::to_string( i ) , true );
                printPopulationToFile( isl.get_population().get_f(), "targetingPropagation_" + std::to_string(currentQuestion ) + "_" + std::to_string(currentCase) + "_" + std::to_string(i), true);
            }
                //std::cout<<i<<" "<<algorithmIndex<<std::endl;


            }





        }
	}
	// The exit code EXIT_SUCCESS indicates that the program was successfully executed.
	return EXIT_SUCCESS;
}
