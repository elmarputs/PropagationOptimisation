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



/*
namespace tudat
{

//! Class to set the aerodynamic angles of the capsule (default: all angles 0)
class CapsuleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor
    CapsuleAerodynamicGuidance(
            const NamedBodyMap bodyMap,
            const double fixedAngleOfAttack ):bodyMap_( bodyMap ), fixedAngleOfAttack_( fixedAngleOfAttack )
    {
    }

    //! The aerodynamic angles are to be computed here
    void updateGuidance( const double time ) override
    {
        currentAngleOfAttack_ = fixedAngleOfAttack_;
        currentAngleOfSideslip_ = 0.0;
        currentBankAngle_ = 0.0;

    }

private:

    //! List of body objects that constitute the environment
    NamedBodyMap bodyMap_;

    //! Fixed angle of attack that is to be used by vehicle
    double fixedAngleOfAttack_;
};

}
*/

/*!
 *   This function computes the entry trajectory of a capsule, where the shape of the capsule is used to determine the vehicle's
 *   aerodynamic force and moment coefficients. The aerodynamic coefficients are based on local inclination methods, and computed
 *   by an object of the HypersonicLocalInclinationAnalysis class.
 *
 *   The present code uses only an aerodynamic force, and a point-mass Earth gravitational acceleration. A class
 *   CapsuleAerodynamicGuidance has been provided, which is currently has no direct functionality: it sets the aerodynamic angles
 *   (attack, sideslip, bank) to 0 degrees. This can (and should) be overridden by the user in favor of something more realistic
 *
 *   Key outputs:
 *
 *   propagatedStateHistory Numerically propagated Cartesian state
 *   dependentVariableHistory Dependent variables (default none) saved during the state propagation of the entry capsule
 *
 *   Input parameters:
 *
 *   shapeParameters: A vector defining the problem as follows (first five parameters describe shape; sixth paramater related
 *      to guidance): Nose radius, Middle radius, Rear length, Rear angle
 *      Side radius, Constant Angle of Attack (see Dirkx and Mooij, 2018 for more details)
 *
 */

//! Execute propagation of orbits of Capsule during entry.

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
		// Convert Eigen vector to std::vector for use with fitness function
		const std::vector<double> controlVars{&cv[0], cv.data() + cv.size()};

		// Run propagation with control variables for current design iteration and save to map
		std::vector<double> fitness = prob.fitness(controlVars);
		Eigen::VectorXd objectives(1);
		objectives << fitness.at(0);
		Eigen::VectorXd outputVector(cv.size() + objectives.size());
		outputVector << cv, objectives;
		ccdFitnessMap[design] = outputVector;
	}

	// Print map with fitness values to file
	tudat::input_output::writeDataMapToTextFile(ccdFitnessMap, "ccdFitness.txt", filePath);

}

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
