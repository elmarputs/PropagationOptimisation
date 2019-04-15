#include "shapeOptimization.h"

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

// Static member variables
simulation_setup::NamedBodyMap ShapeOptimization::bodyMap_;
Eigen::Vector6d ShapeOptimization::capsuleSphericalEntryState;
double ShapeOptimization::simulationStartEpoch = 0.0;
double ShapeOptimization::simulationEndEpoch = 2*86400;
bool ShapeOptimization::isSpiceLoaded = false;

// Heat load derivative function
inline double tudat::ShapeOptimization::heatLoadFunction( const double time, const double state, NamedBodyMap& bodyMap_ ) const
{
	std::shared_ptr<aerodynamics::AtmosphericFlightConditions> atmosphericFlightConditions;
	atmosphericFlightConditions = std::dynamic_pointer_cast<aerodynamics::AtmosphericFlightConditions>(bodyMap_.at("Capsule")->getFlightConditions());

	double heatRate = atmosphericFlightConditions->getCurrentAerodynamicHeatRate();

	//std::cout << "Current heat rate: " << heatRate << " W/m^2\n";

	return heatRate;
}


tudat::ShapeOptimization::ShapeOptimization()
{
	std::cout << "ShapeOptimization problem constructed.\n";
	//std::cout << "Setting up Tudat...\n";

	if( !isSpiceLoaded )
	{
		spice_interface::loadStandardSpiceKernels( );
		std::cout<<"Created Spice Kernels \n";
		isSpiceLoaded = true;
	}

	// Set spherical elements for Capsule
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


	std::vector<std::string> bodiesToCreate;
	bodiesToCreate.push_back("Earth");

	std::map<std::string, std::shared_ptr<simulation_setup::BodySettings> > bodySettings_;

	// Set Earth gravity field settings.
	bodySettings_ =
			getDefaultBodySettings(bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0);

	// Change default parameters based on earlier analyses
	bodySettings_["Earth"]->gravityFieldSettings =
			std::make_shared<FromFileSphericalHarmonicsGravityFieldSettings>(ggm02s);

	bodySettings_["Earth"]->atmosphereSettings = std::make_shared<AtmosphereSettings>(nrlmsise00);

	double bodyRadius = 6378.0E3;
	double bodyFlattening = 1.0 / 300.0;
	bodySettings_["Earth"]->shapeModelSettings = std::make_shared<OblateSphericalBodyShapeSettings>(bodyRadius, bodyFlattening);

	bodySettings_["Earth"]->rotationModelSettings->resetOriginalFrame("J2000");
	bodySettings_["Earth"]->ephemerisSettings->resetFrameOrientation("J2000");

	// Create Earth object
	bodyMap_ = simulation_setup::createBodies(bodySettings_);

}

vector_double tudat::ShapeOptimization::fitness(const vector_double& decisionVariables) const
{
	//std::cout << "Fitness\n";

	// Vehicle properties
	double vehicleDensity = 250.0;

	bodyMap_[ "Capsule" ] = std::make_shared< simulation_setup::Body >( );
	std::shared_ptr<Body> capsulePtr = bodyMap_["Capsule"];

	// Finalize body creation.
	setGlobalFrameBodyEphemerides( bodyMap_, "Earth", "J2000" );

    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimisationGroup" );

    vector_double fitness;

    std::vector<double> shapeParameters = decisionVariables;
	//std::vector< double > shapeParameters =
	//		{ 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };
    double limitLength = ( shapeParameters[ 1 ] - shapeParameters[ 4 ] * ( 1.0 - std::cos( shapeParameters[ 3 ] ) ) ) /
            std::tan( -shapeParameters[ 3 ] );
    if( shapeParameters[ 2 ] >= limitLength - 0.01 )
    {
        shapeParameters[ 2 ] = limitLength -0.01;
    }

    // Create capsule.
    std::shared_ptr< geometric_shapes::Capsule > capsule
            = std::make_shared< geometric_shapes::Capsule >(
                shapeParameters[ 0 ], shapeParameters[ 1 ], shapeParameters[ 2 ],
            shapeParameters[ 3 ], shapeParameters[ 4 ] );

    bodyMap_.at("Capsule")->setConstantBodyMass(
                capsule->getVolume( ) * vehicleDensity );

    // Create vehicle aerodynamic coefficients
    bodyMap_.at("Capsule")->setAerodynamicCoefficientInterface(
                getCapsuleCoefficientInterface( capsule, outputPath, "output_", true ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////            CREATE ACCELERATIONS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap_, shapeParameters.at( 5 ) );
    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap_.at( "Capsule" ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////            CREATE PROPAGATION SETTINGS            /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
	dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
			stagnation_point_heat_flux_dependent_variable, "Capsule" ) );
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
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );


    // Create propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                terminationSettings, cowell, dependentVariablesToSave );

    double initialHeatLoad = 0;

    std::function<double (double, double)> heatLoadDerivativeFunction = std::bind(&tudat::ShapeOptimization::heatLoadFunction, this,
    		std::placeholders::_1, std::placeholders::_2, bodyMap_);

    // Heat load propagator
    std::shared_ptr< CustomStatePropagatorSettings<double, double> > heatLoadPropagatorSettings =
    		std::make_shared< CustomStatePropagatorSettings<double, double> >(heatLoadDerivativeFunction, initialHeatLoad,
    				terminationSettings);

	std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
	propagatorSettingsVector.push_back( translationalPropagatorSettings );
	propagatorSettingsVector.push_back( heatLoadPropagatorSettings );

	std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
			std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            PERFORM PROPAGATION            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
	SingleArcDynamicsSimulator< > dynamicsSimulator(
			bodyMap_, integratorSettings, propagatorSettings, true, true );

    std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

	double propagationTime = propagatedStateHistory.rbegin()->first;
	Eigen::Vector7d finalState = propagatedStateHistory.rbegin()->second;
	Eigen::Vector7d initialState = propagatedStateHistory.begin()->second;
	double integratedHeatRateFinal = finalState(6);
	double integratedHeatRateInit = initialState(6);

	//std::cout << "Initial heat load: " << integratedHeatRateInit << ", final: " << integratedHeatRateFinal << "\n";

    fitness.push_back(integratedHeatRateFinal);
    return fitness;
}

std::pair<vector_double, vector_double> tudat::ShapeOptimization::get_bounds() const
{

    std::pair<vector_double, vector_double> bounds;
    bounds.first = {7.9, 2.5, 0.15, -1.0, 0.15, 0.03};
    bounds.second = {8.4, 3.5, 0.80, -0.3, 0.35, 0.05};

    return bounds;
}

pagmo::thread_safety tudat::ShapeOptimization::get_thread_safety() const
{
	return pagmo::thread_safety::none;
}
