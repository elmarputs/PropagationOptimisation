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

#include "../applicationOutput.h"
#include "shapeOptimization.h"

#include <iostream>

using namespace tudat;
using namespace tudat::aerodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::mathematical_constants;

/*!
 *  Function that creates an aerodynamic database for a capsule, based on a set of shape parameters
 *  The Capsule shape consists of four separate geometrical components: a sphere segment for the nose, a torus segment for the
 *  shoulder/edge, a conical frustum for the rear body, and a sphere segment for the rear cap (see Dirkx and Mooij, 2016).
 *  The code used in this function discretizes these surfaces into a structured mesh of quadrilateral panels. The parameters
 *  numberOfPoints and numberOfLines define the number of discretization points (for each part) in both independent directions
 *  (lengthwise and circumferential).
 */
std::shared_ptr< HypersonicLocalInclinationAnalysis > getCapsuleCoefficientInterface(
        const std::shared_ptr< geometric_shapes::Capsule > capsule,
        const std::string directory,
        const std::string filePrefix,
        const bool useNewtonianMethodForAllPanels = true )
{

    // Define settings for surface discretization of capsule
    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 31;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;

    // DO NOT CHANGE THESE (setting to true will turn parts of the vehicle 'inside out')
    std::vector< bool > invertOrders;
    invertOrders.resize( 4 );
    invertOrders[ 0 ] = 0;
    invertOrders[ 1 ] = 0;
    invertOrders[ 2 ] = 0;
    invertOrders[ 3 ] = 0;

    // Define moment reference point
    Eigen::Vector3d momentReference;
    momentReference( 0 ) = -0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;

    // Define independent variable values
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 15 );
    for ( int i = 0; i < 15; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }
    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

    // Define local inclination methods to use
    std::vector< std::vector< int > > selectedMethods;
    selectedMethods.resize( 2 );
    selectedMethods[ 0 ].resize( 4 );
    selectedMethods[ 1 ].resize( 4 );
    if( !useNewtonianMethodForAllPanels )
    {
        selectedMethods[ 0 ][ 0 ] = 1;
        selectedMethods[ 0 ][ 1 ] = 5;
        selectedMethods[ 0 ][ 2 ] = 5;
        selectedMethods[ 0 ][ 3 ] = 1;
        selectedMethods[ 1 ][ 0 ] = 6;
        selectedMethods[ 1 ][ 1 ] = 3;
        selectedMethods[ 1 ][ 2 ] = 3;
        selectedMethods[ 1 ][ 3 ] = 3;
    }
    else
    {
        selectedMethods[ 0 ][ 0 ] = 0;
        selectedMethods[ 0 ][ 1 ] = 0;
        selectedMethods[ 0 ][ 2 ] = 0;
        selectedMethods[ 0 ][ 3 ] = 0;
        selectedMethods[ 1 ][ 0 ] = 0;
        selectedMethods[ 1 ][ 1 ] = 0;
        selectedMethods[ 1 ][ 2 ] = 0;
        selectedMethods[ 1 ][ 3 ] = 0;
    }

    // Create aerodynamic database
    std::shared_ptr< HypersonicLocalInclinationAnalysis > hypersonicLocalInclinationAnalysis =
            std::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * std::pow( capsule->getMiddleRadius( ), 2.0 ),
                capsule->getMiddleRadius( ), momentReference, false );

    // Save vehicle mesh to a file
    aerodynamics::saveVehicleMeshToFile(
                hypersonicLocalInclinationAnalysis, directory, filePrefix );

    // Create analysis object and capsule database.
    return  hypersonicLocalInclinationAnalysis;
}

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
int main( )
{
	tudat::spice_interface::loadStandardSpiceKernels();

	std::vector<std::string> bodiesToCreate;
	bodiesToCreate.push_back("Earth");
	bodiesToCreate.push_back("Sun");
	bodiesToCreate.push_back("Moon");
	std::map<std::string, std::shared_ptr<BodySettings>> bodySettings = getDefaultBodySettings(bodiesToCreate, 0.0-300.0, 24.0*3600.0);
	for(unsigned int i = 0; i < bodiesToCreate.size(); i++)
	{
		bodySettings[bodiesToCreate.at(i)]->rotationModelSettings->resetOriginalFrame("J2000");
		bodySettings[bodiesToCreate.at(i)]->ephemerisSettings->resetFrameOrientation("J2000");
	}
	double doubleValue;
	std::shared_ptr<double> doublePtr = std::make_shared<double>(doubleValue);

	std::cout << *doublePtr << std::endl;
	std::cout << doublePtr << std::endl;

	bodySettings["Earth"]->shapeModelSettings = std::make_shared<OblateSphericalBodyShapeSettings>(6378e3, 1.0/300.0);
	// The exit code EXIT_SUCCESS indicates that the program was successfully executed.
	return EXIT_SUCCESS;
}
