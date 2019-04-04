#include "shapeOptimization.h"
#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <pagmo/problem.hpp>

#include "../applicationOutput.h"

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


tudat::ShapeOptimization::ShapeOptimization()
{




}

vector_double tudat::ShapeOptimization::fitness(const std::vector<double> &decisionVariables) const
{
    // Set spherical elements for Capsule.
    Eigen::Vector6d capsuleSphericalEntryState;
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 0.0 );
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 68.75 );
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.83E03;
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            unit_conversions::convertDegreesToRadians( -1.5 );
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) =
            unit_conversions::convertDegreesToRadians( 34.37 );

    // Vehicle properties
    double vehicleDensity = 250.0;

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > shapeParameters =
            { 7.353490006527863,	 2.99718480813317, 	 1.006902199215256,	 -0.9043838974848388,	 0.4513045128015801,	 0.02889224366077636 };

    // Set simulation start epoch.
    double simulationStartEpoch = 0.0;

    // Set simulation start epoch.
    double simulationEndEpoch = 2*86400;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back("Earth");


    // Set Earth gravity field settings.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );



    // Change default parameters based on earlier analyses
    bodySettings[ "Earth" ]->gravityFieldSettings =
            std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02s );

    bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< AtmosphereSettings >( nrlmsise00 );

    double bodyRadius = 6378.0E3;
    double bodyFlattening = 1.0 / 300.0;
    bodySettings[ "Earth" ]->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >( bodyRadius, bodyFlattening );

    bodySettings[ "Earth"]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth"]->ephemerisSettings->resetFrameOrientation( "J2000" );


    // Create Earth object
    bodyMap_ = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap_[ "Capsule" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap_, "Earth", "J2000" );

    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimisationGroup" );
	std::cout << "Fitness function called!\n";

    std::cout<<"Starting fitness \n";

    vector_double fitness;

    shapeParameters = decisionVariables;
    double limitLength = ( shapeParameters[ 1 ] - shapeParameters[ 4 ] * ( 1.0 - std::cos( shapeParameters[ 3 ] ) ) ) /
            std::tan( -shapeParameters[ 3 ] );
    if( shapeParameters[ 2 ] >= limitLength - 0.01 )
    {
        shapeParameters[ 2 ] = limitLength -0.01;
    }

    std::cout << "Creating capsule...\n";

    // Create capsule.
    std::shared_ptr< geometric_shapes::Capsule > capsule
            = std::make_shared< geometric_shapes::Capsule >(
                shapeParameters[ 0 ], shapeParameters[ 1 ], shapeParameters[ 2 ],
            shapeParameters[ 3 ], shapeParameters[ 4 ] );

    bodyMap_.at("Capsule")->setConstantBodyMass(
                capsule->getVolume( ) * vehicleDensity );

    std::cout<<"Creating Aerodynamics \n";


    // Create vehicle aerodynamic coefficients
    bodyMap_.at("Capsule")->setAerodynamicCoefficientInterface(
                getCapsuleCoefficientInterface( capsule, outputPath, "output_", true ) );

	std::cout << "Aerodynamic coefficient interface set\n";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////            CREATE ACCELERATIONS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    std::cout<<"Creating Accelerations \n";
    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCapsule;

    accelerationsOfCapsule["Earth"].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(2,0) );
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

    accelerationMap[ "Capsule" ] = accelerationsOfCapsule;

    bodiesToPropagate.push_back( "Capsule" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );

    std::cout << "Created acceleration map.\n";

    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap_, shapeParameters.at( 5 ) );
    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap_.at( "Capsule" ) );

	std::cout << "Guidance angle function set\n";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////            CREATE PROPAGATION SETTINGS            /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Creating Propagator Settings \n";
    // Set Integrator
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
            (simulationStartEpoch, 0.5, RungeKuttaCoefficients::rungeKuttaFehlberg56, std::numeric_limits<double>::epsilon(), 4.0, 1e-5, 1e-5 );


    // Convert capsule state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                capsuleSphericalEntryState );

    // Convert the state to the global (inertial) frame.
    systemInitialState = transformStateToGlobalFrame(
                systemInitialState, simulationStartEpoch, bodyMap_.at( "Earth" )->getRotationalEphemeris( ) );

    // Define termination conditions
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Capsule", "Earth" );
    terminationSettingsList.push_back(
                std::make_shared< PropagationDependentVariableTerminationSettings >(
                    terminationDependentVariable, 25.0E3, true ) );
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           24.0 * 3600.0 ) );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = std::make_shared<
            PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          airspeed_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          aerodynamic_force_coefficients_dependent_variable, "Capsule" ) );
    // Add variables to get latitude and longitude
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::latitude_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::longitude_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::heading_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          total_aerodynamic_g_load_variable, "Capsule", "Earth" ) );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


    // Create propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings<  > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                terminationSettings, cowell, dependentVariablesToSave );

    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap_, integratorSettings, propagatorSettings );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            PERFORM PROPAGATION            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Propagate Entry \n";

    // Create simulation object and propagate dynamics.

    std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

	std::shared_ptr<PropagationTerminationDetails > propagationTerminationDetails = dynamicsSimulator.getPropagationTerminationReason();
	std::cout << "Termination reason: " << propagationTerminationDetails->getPropagationTerminationReason() << std::endl;

    fitness.push_back(0.0);
    return fitness;
}

std::pair<vector_double, vector_double> tudat::ShapeOptimization::get_bounds() const
{

    std::pair<vector_double, vector_double> bounds;
    bounds.first = {7.0, 2.5, 0.5, -1.5, 0.0, 0.01};
    bounds.second = {8.0, 3.5, 1.5, -0.5, 1.0, 0.04};

    return bounds;
}



