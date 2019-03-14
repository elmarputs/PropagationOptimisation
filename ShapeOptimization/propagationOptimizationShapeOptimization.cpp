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
using namespace tudat::reference_frames;
using namespace tudat;

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
    void updateGuidance( const double time )
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
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimization" );

    bool runMonteCarlo = true;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::function< double( ) > altitudeErrorFunction = tudat::statistics::createBoostContinuousRandomVariableGeneratorFunction(
            tudat::statistics::normal_boost_distribution, {0, 1000}, 0 );
    std::function< double( ) > velocityErrorFunction = tudat::statistics::createBoostContinuousRandomVariableGeneratorFunction(
            tudat::statistics::normal_boost_distribution, {0, 250}, 0 );
    std::function< double( ) > gammaErrorFunction = tudat::statistics::createBoostContinuousRandomVariableGeneratorFunction(
            tudat::statistics::normal_boost_distribution, {0, 0.2}, 0 );

    std::function< double( ) > vehicleDensityErrorFunction = tudat::statistics::createBoostContinuousRandomVariableGeneratorFunction(
            tudat::statistics::uniform_boost_distribution, {-50.0, 50.0}, 0 );

    // Vehicle properties
    double vehicleDensity = 250.0;
    double area = 22.0;
    double radiationPressureCoefficient = 1.0;

    //std::cout << "Enter a constant angle of attack: \n";
    //std::cin >> angleOfAttack;

    std::vector< double > shapeParameters =
    { 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };

    std::vector<std::map<double, Eigen::VectorXd> > propagationHistories;
    std::vector<std::map<double, Eigen::VectorXd> > dependentVariableHistories;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int cases = 10;

    for( unsigned int currentCase = 9; currentCase < cases; currentCase++ )
    {
        std::cout << "Running case " << currentCase << ".\n";
        // Create environment
        // Define simulation body settings.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Moon" );
        bodiesToCreate.push_back( "Jupiter" );
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 3600.0, simulationStartEpoch + 24.0 * 3600.0 );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        }

        std::shared_ptr< RadiationPressureInterfaceSettings > capsuleRadiationPressureSettings;

        switch( currentCase )
        {
            case 2:
            {
                // NRLMSISE00 atmosphere model
                std::cout << "NRLMSISE00 atmosphere model\n";
                bodySettings["Earth"]->atmosphereSettings = std::make_shared<AtmosphereSettings>(nrlmsise00);
                break;
            }
            case 5:
            {
                std::cout << "Radiation pressure of the Sun\n";
                std::vector<std::string> occultingBodies;bodySettings["Earth"]->atmosphereSettings = std::make_shared<AtmosphereSettings>(nrlmsise00);
                occultingBodies.push_back("Earth");
                capsuleRadiationPressureSettings = std::make_shared<
                        CannonBallRadiationPressureInterfaceSettings>("Sun", area, radiationPressureCoefficient, occultingBodies);
                break;
            }
            case 6:
            {
                // Oblate Earth as shape model
                bodySettings["Earth"]->shapeModelSettings = std::make_shared<OblateSphericalBodyShapeSettings>(6378e3, 1.0/300.0);
                break;
            }
            case 9:
            {
                bodySettings["Earth"]->atmosphereSettings = std::make_shared<AtmosphereSettings>(nrlmsise00);
                bodySettings["Earth"]->shapeModelSettings = std::make_shared<OblateSphericalBodyShapeSettings>(6378e3, 1.0/300.0);
                break;
            }
            default:
                break;
        }

        // Create Earth object
        simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );
        // Create vehicle objects.
        bodyMap["Capsule"] = std::make_shared<simulation_setup::Body>();

        // Finalize body creation.
        setGlobalFrameBodyEphemerides(bodyMap, "Earth", "J2000");

        double limitLength = (shapeParameters[1] - shapeParameters[4] * (1.0 - std::cos(shapeParameters[3]))) /
                             std::tan(-shapeParameters[3]);
        if (shapeParameters[2] >= limitLength - 0.01)
        {
            shapeParameters[2] = limitLength - 0.01;
        }

        // Create capsule.
        std::shared_ptr<geometric_shapes::Capsule> capsule
                = std::make_shared<geometric_shapes::Capsule>(
                        shapeParameters[0], shapeParameters[1], shapeParameters[2],
                        shapeParameters[3], shapeParameters[4]);

        bodyMap["Capsule"]->setConstantBodyMass(capsule->getVolume() * vehicleDensity);

        // Create vehicle aerodynamic coefficients
        bodyMap["Capsule"]->setAerodynamicCoefficientInterface(
                getCapsuleCoefficientInterface(capsule, outputPath, "output_", true));

        if( currentCase == 5 )
        {
            bodyMap["Capsule"]->setRadiationPressureInterface("Sun", createRadiationPressureInterface(capsuleRadiationPressureSettings, "Capsule",
                    bodyMap));
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector<std::string> bodiesToPropagate;
        std::vector<std::string> centralBodies;

        // Define acceleration model settings.
        std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfCapsule;

        switch( currentCase )
        {
            case 0:
                accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared<AccelerationSettings>(central_gravity) );
                break;
            case 1:
                accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
                break;
            case 2:
                // NRLMSISE00 atmosphere model
                accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                break;
            case 3:
                // 3rd body effect of Moon on Capsule
                accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                accelerationsOfCapsule[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                break;
            case 4:
                accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                accelerationsOfCapsule["Sun"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                break;
            case 5:
                accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                accelerationsOfCapsule["Sun"].push_back(std::make_shared<AccelerationSettings>(cannon_ball_radiation_pressure));
                break;
            case 6:
                accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                break;
            case 7:
                accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                accelerationsOfCapsule["Jupiter"].push_back(std::make_shared<AccelerationSettings>(central_gravity));
                break;
            case 8:
                accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                break;
            case 9:
                accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
            default:
                break;
        }

        accelerationsOfCapsule["Earth"].push_back(std::make_shared<AccelerationSettings>(aerodynamic));

        accelerationMap["Capsule"] = accelerationsOfCapsule;

        bodiesToPropagate.push_back("Capsule");
        centralBodies.push_back("Earth");

        // Create acceleration models
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies);

        std::shared_ptr<CapsuleAerodynamicGuidance> capsuleGuidance =
                std::make_shared<CapsuleAerodynamicGuidance>(bodyMap, shapeParameters.at(5));
        setGuidanceAnglesFunctions(capsuleGuidance, bodyMap.at("Capsule"));

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



        // Define termination conditions
        std::vector<std::shared_ptr<PropagationTerminationSettings> > terminationSettingsList;
        std::shared_ptr<SingleDependentVariableSaveSettings> terminationDependentVariable =
                std::make_shared<SingleDependentVariableSaveSettings>(altitude_dependent_variable, "Capsule", "Earth");
        terminationSettingsList.push_back(
                std::make_shared<PropagationDependentVariableTerminationSettings>(
                        terminationDependentVariable, 25.0E3, true));
        terminationSettingsList.push_back(std::make_shared<PropagationTimeTerminationSettings>(
                24.0 * 3600.0));
        std::shared_ptr<PropagationTerminationSettings> terminationSettings = std::make_shared<
                PropagationHybridTerminationSettings>(terminationSettingsList, true);

        // Define dependent variables
        std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariablesList;
        dependentVariablesList.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
                altitude_dependent_variable, "Capsule", "Earth"));
        dependentVariablesList.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
                airspeed_dependent_variable, "Capsule", "Earth"));
        dependentVariablesList.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
                aerodynamic_force_coefficients_dependent_variable, "Capsule"));
        dependentVariablesList.push_back(std::make_shared<BodyAerodynamicAngleVariableSaveSettings>("Capsule", latitude_angle));
        dependentVariablesList.push_back(std::make_shared<BodyAerodynamicAngleVariableSaveSettings>("Capsule", longitude_angle));
        dependentVariablesList.push_back(std::make_shared<SingleAccelerationDependentVariableSaveSettings>(aerodynamic, "Capsule", "Earth", false));
        std::shared_ptr<DependentVariableSaveSettings> dependentVariablesToSave =
                std::make_shared<DependentVariableSaveSettings>(dependentVariablesList);

        // Define propagator type
        TranslationalPropagatorType propagatorType;
        std::shared_ptr<IntegratorSettings<> > integratorSettings;

        std::shared_ptr<interpolators::OneDimensionalInterpolator<double, Eigen::VectorXd> > interpolator;
        std::shared_ptr<interpolators::OneDimensionalInterpolator<double, Eigen::VectorXd> > depVarInterpolator;

        // Benchmark
        if(currentCase == -1)
        {
            integratorSettings = std::make_shared<RungeKuttaVariableStepSizeSettings<> >(simulationStartEpoch, 0.04,
                    RungeKuttaCoefficients::rungeKuttaFehlberg78, 0.005, 0.5, 1e-14, 1e-14);
        }
        else
        {
            integratorSettings = std::make_shared<RungeKuttaVariableStepSizeSettings<> >(simulationStartEpoch, 0.04,
                    RungeKuttaCoefficients::rungeKuttaFehlberg56, 0.005, 0.5, 1e-5, 1e-5);
        }

        propagatorType = cowell;

        // Double for loop for Monte Carlo
        unsigned int numberOfInputs = 1;
        if( runMonteCarlo )
        {
            numberOfInputs = 4; // How many initial conditions are to be analysed
        }

        for( unsigned int i = 0; i < numberOfInputs; i++ )
        {
            Eigen::Vector6d capsuleSphericalEntryState;
            capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
                    spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
            capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
                    unit_conversions::convertDegreesToRadians( 0.0 );
            capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) =
                    unit_conversions::convertDegreesToRadians( 68.75 );
            capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 6.5E3;
            capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
                    unit_conversions::convertDegreesToRadians( -1.5 );
            capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) =
                    unit_conversions::convertDegreesToRadians( 34.37 );

            std::map< int, Eigen::VectorXd > finalStates;
            std::map< int, Eigen::VectorXd > finalDependentVariables;

            unsigned int numberOfSamples = 250;
            if( i == 0 )
            {
                // Only one run to obtain reference trajectory
                numberOfSamples = 1;
            }

            for( unsigned int j = 0; j < numberOfSamples; j++ )
            {
                std::cout << "MC run " << i << " : " << j << "\n";
                // Initialise initial conditions (= systemInitialState)
                // Do nothing if i == 0 (reference trajectory)
                if( i == 1 )
                {
                    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) = spice_interface::getAverageRadius( "Earth" ) +
                            120.0E3 + altitudeErrorFunction();
                }
                else if( i == 2 )
                {
                    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 6.5E3 + velocityErrorFunction();
                }
                else if( i == 3 )
                {
                    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) = unit_conversions::convertDegreesToRadians(
                            -1.5 + gammaErrorFunction( ) );
                }
                else if( i == 4 )
                {

                }

                // Convert capsule state from spherical elements to Cartesian elements.
                Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                        capsuleSphericalEntryState);

                // Convert the state to the global (inertial) frame.
                systemInitialState = transformStateToGlobalFrame(
                        systemInitialState, simulationStartEpoch, bodyMap.at("Earth")->getRotationalEphemeris());


                // Set up propagator settings with new initial conditions
                std::shared_ptr<TranslationalStatePropagatorSettings<> > propagatorSettings =
                        std::make_shared<TranslationalStatePropagatorSettings<> >(
                                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                                terminationSettings, propagatorType, dependentVariablesToSave);

                // Perform propagation and extract solution map
                SingleArcDynamicsSimulator<> dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings);
                std::map<double, Eigen::VectorXd> propagationHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
                std::map<double, Eigen::VectorXd> dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory();

                // Extract last state (which is relevant for error analysis)
                finalStates[ j ] = propagationHistory.rbegin()->second;
                finalDependentVariables[ j ] = dependentVariableHistory.rbegin()->second;
            }

            // Write map with final states to file
            input_output::writeDataMapToTextFile( finalStates, "MC_" + std::to_string(i) + ".dat", outputPath );
            input_output::writeDataMapToTextFile( finalDependentVariables, "MC_" + std::to_string(i) + "_dep.dat", outputPath );
        }

//        std::map<double, Eigen::VectorXd> propagationHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
//        std::cout << "Number of function evaluations: " << dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations().rbegin()->second << std::endl;
//        propagationHistories.push_back(propagationHistory);
//        std::map<double, Eigen::VectorXd> dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory();
//        dependentVariableHistories.push_back(dependentVariableHistory);

//        input_output::writeDataMapToTextFile(propagationHistory, "stateHistory_" + std::to_string(currentCase) + ".dat", outputPath);
//        input_output::writeDataMapToTextFile(dependentVariableHistory, "dependentVariables_" + std::to_string(currentCase) + ".dat", outputPath);

        // Create interpolator for map of current run
//        std::shared_ptr<interpolators::InterpolatorSettings> interpolatorSettings = std::make_shared<
//                interpolators::LagrangeInterpolatorSettings>(10);
//        interpolator = interpolators::createOneDimensionalInterpolator(propagationHistory, interpolatorSettings);
//        depVarInterpolator = interpolators::createOneDimensionalInterpolator(dependentVariableHistory, interpolatorSettings);
//
//        // Interpolate to time steps of first (i.e., benchmark) run
//        std::map<double, Eigen::VectorXd> interpolatedState;
//        for (const auto &stateIterator : propagationHistories.at(0))
//        {
//            interpolatedState[stateIterator.first] = interpolator->interpolate(stateIterator.first);
//        }
//        std::map<double, Eigen::VectorXd> interpolatedDepVar;
//        for (const auto &stateIterator : dependentVariableHistories.at(0))
//        {
//            interpolatedDepVar[stateIterator.first] = depVarInterpolator->interpolate(stateIterator.first);
//        }
//
//        input_output::writeDataMapToTextFile(interpolatedState, "stateHistory_interpolated_" + std::to_string(currentCase) + ".dat", outputPath);
//        input_output::writeDataMapToTextFile(interpolatedDepVar, "dependentVariables_interpolated_" + std::to_string(currentCase) + ".dat", outputPath);

    }
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
