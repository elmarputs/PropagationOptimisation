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
        ////                                                                                ////
        ////    Retrieve relevant objects from environment and set as member variables      ////
        ////                                                                                ////

        vehicleFlightConditions_ = std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( "Capsule" )->getFlightConditions( ) );

        earthBody_ = bodyMap.at("Earth");
        //rotationalVelocityEarth_ = bodyMap.at( "Earth" )->getCurrentAngularVelocityVectorInGlobalFrame( )(2);
        //rotationalVelocityEarth_ = rotationalVelocityVector.norm( );

        vehicleMass_ = bodyMap.at( "Capsule" )->getBodyMass( );
        mu_ = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );


    }

    //! The aerodynamic angles are to be computed here
    void updateGuidance( const double time ) override
    {
        // TEMPORARY FIX TO OBTAIN OUTPUTS
        //fixedAngleOfAttack_ = 0.048619431524615;
        rotationalVelocityEarth_ = earthBody_->getCurrentAngularVelocityVectorInGlobalFrame()(2);

        // Algorithm for angle of attack
        double machNumber = vehicleFlightConditions_->getCurrentMachNumber( );

        if( vehicleFlightConditions_->getCurrentAltitude( ) > 150.0e3 ||  machNumber > 12.0 )
        {
            currentAngleOfAttack_ = fixedAngleOfAttack_; // * mathematical_constants::PI / 180.0;
        }
        else if( machNumber < 6.0 )
        {
            currentAngleOfAttack_ = 0 * mathematical_constants::PI / 180.0;
        }
        else
        {
            double B = (fixedAngleOfAttack_*180/mathematical_constants::PI - 0 - pow(6,3)*0)/(12*(1-pow(6,3)*12/(12*108)));
            double A = (0 - 12*B)/108;
            double aoa = -1*(A*pow((machNumber - 6.0),3) + B*pow((machNumber-6.0),2) + 0);
            currentAngleOfAttack_ = aoa * mathematical_constants::PI / 180.0;
        }

        //    currentAngleOfAttack_ = fixedAngleOfAttack_;


        double airspeed = vehicleFlightConditions_->getCurrentAirspeed();

        // Calculate lift
        Eigen::Vector3d forceCoefficients = vehicleFlightConditions_->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients( );
        double referenceArea = vehicleFlightConditions_->getAerodynamicCoefficientInterface( )->getReferenceArea( );
        double dynamicPressure = vehicleFlightConditions_->getCurrentDynamicPressure( );

        // forceCoefficients(0) -> C_L?
        double liftForce = forceCoefficients(2) * dynamicPressure * referenceArea;


        // Calculate downward component of g
        double r = vehicleFlightConditions_->getCurrentBodyCenteredBodyFixedState().head(3).norm();
        double gd = mu_ / (r * r);


        double latitude = vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        double heading = vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        double flightpath = vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );

        // Solve equation 3 for bank angle (sigma) to minimise gamma dot
        //double sigmaMax = tudat::mathematical_constants::PI/2;

        double sigma = 0.0;


        double coriolis = 2 * rotationalVelocityEarth_ * airspeed * cos(latitude) * sin(heading);
        double centripetal = rotationalVelocityEarth_ * rotationalVelocityEarth_ * r * cos(latitude) * (cos(latitude)*cos(flightpath)+sin(flightpath)*sin(latitude)*cos(heading));
        double weightcomp = - ( gd - ( airspeed * airspeed / r ) ) * cos(flightpath);

        double cosCalculations = ( weightcomp + coriolis + centripetal ) / ( liftForce/vehicleMass_ );
        // Loop defining how we choose the bank angle
        if( question_ == 2 && case_ == 1 )
        {
            sigma = 0.0;
        }
        else if( cosCalculations > 1 )
        {
            sigma = tudat::mathematical_constants::PI;
        }
        else if( cosCalculations < -1 )
        {
            sigma = 0.0;
        }
        else
        {
            sigma = std::acos( - cosCalculations );
        }


        currentBankAngle_ = sigma;
        //currentBankAngle_ = 0.0;




        currentAngleOfSideslip_ = 0.0;
    }

private:

    //! List of body objects that constitute the environment
    NamedBodyMap bodyMap_;

    //! Fixed angle of attack that is to be used by vehicle
    double fixedAngleOfAttack_;

    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > vehicleFlightConditions_;
    double rotationalVelocityEarth_;
    double vehicleMass_;
    double mu_;
    std::shared_ptr< simulation_setup::Body >earthBody_;
    unsigned int question_;
    unsigned int case_;


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

    //double heatRate = atmosphericFlightConditions->getCurrentAerodynamicHeatRate();

    double heatRateRF = propagators::computeEquilibriumFayRiddellHeatFluxFromProperties(atmosphericFlightConditions, bodyMap_.at("Capsule")
    ->getVehicleSystems());

	//std::cout << "RF heat rate: " << heatRateRF << " W/m2. Via flight conditions: " << heatRate << " W/m2 \n";

    return heatRateRF;

	//std::cout << "Current heat rate: " << heatRate << " W/m^2\n";

    //return heatRate;
}

inline double tudat::ShapeOptimization::flightRangeFunction( const double time, const double state, NamedBodyMap& bodyMap_ ) const
{
    std::shared_ptr<aerodynamics::AtmosphericFlightConditions> atmosphericFlightConditions;
    atmosphericFlightConditions = std::dynamic_pointer_cast<aerodynamics::AtmosphericFlightConditions>(bodyMap_.at("Capsule")->getFlightConditions());

    double rangeVelocity = atmosphericFlightConditions->getCurrentAirspeedBasedVelocity().norm();
    //std::cout<<"Current Velocity:"<<rangeVelocity<<"\n";
    return rangeVelocity;
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

    //Temporary values to obtain the flight parameters

    //std::vector< double > shapeParameters =
    //        { 8.395074713222693, 2.718209058425671, 0.200568003033042, -0.302674762298538, 0.349473115033317, 0.048619431524615 };
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

    std::shared_ptr<system_models::VehicleSystems> vehicleSystems = std::make_shared<system_models::VehicleSystems>(bodyMap_.at("Capsule")
    		->getBodyMass());
    vehicleSystems->setWallEmissivity(0.40);
    vehicleSystems->setNoseRadius(shapeParameters[0]);

    bodyMap_.at("Capsule")->setVehicleSystems(vehicleSystems);

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

    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          airspeed_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          aerodynamic_force_coefficients_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            stagnation_point_heat_flux_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            total_aerodynamic_g_load_variable, "Capsule", "Earth") );


    // Add variables to get latitude and longitude
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::latitude_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::longitude_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::heading_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule", tudat::reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle ,"Earth" ) );

    // Add aerodynamic Angles
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule",tudat::reference_frames::AerodynamicsReferenceFrameAngles::bank_angle ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Capsule",tudat::reference_frames::AerodynamicsReferenceFrameAngles::angle_of_attack ,"Earth" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                "Capsule",tudat::reference_frames::AerodynamicsReferenceFrameAngles::angle_of_sideslip , "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            mach_number_dependent_variable, "Capsule", "Earth") );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


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



    // Create propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                terminationSettings, cowell );

    double initialHeatLoad = 0;
    double initialFlightRange = 0;

    std::function<double (double, double)> heatLoadDerivativeFunction = std::bind(&tudat::ShapeOptimization::heatLoadFunction, this,
    		std::placeholders::_1, std::placeholders::_2, bodyMap_);

    std::function<double (double, double)> flightRangeDerivativeFunction = std::bind(&tudat::ShapeOptimization::flightRangeFunction, this,
            std::placeholders::_1, std::placeholders::_2, bodyMap_);

    // Heat load propagator
    std::shared_ptr< CustomStatePropagatorSettings<double, double> > heatLoadPropagatorSettings =
    		std::make_shared< CustomStatePropagatorSettings<double, double> >(heatLoadDerivativeFunction, initialHeatLoad,
                    terminationSettings);

    // Flight range propagator
    std::shared_ptr< CustomStatePropagatorSettings<double, double> > flightRangePropagatorSettings =
            std::make_shared< CustomStatePropagatorSettings<double, double> >(flightRangeDerivativeFunction, initialFlightRange,
                    terminationSettings);

	std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
	propagatorSettingsVector.push_back( translationalPropagatorSettings );
	propagatorSettingsVector.push_back( heatLoadPropagatorSettings );
    propagatorSettingsVector.push_back( flightRangePropagatorSettings );

	std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            PERFORM PROPAGATION            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
	SingleArcDynamicsSimulator< > dynamicsSimulator(
			bodyMap_, integratorSettings, propagatorSettings, true, true );

    std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );


    writeDataMapToTextFile( propagatedStateHistory, "Optimal_Capsule_stateHistory.dat", outputPath );
    writeDataMapToTextFile( dependentVariableHistory, "Optimal_Capsule_dependentVariableHistory.dat", outputPath );


    // Find Propagation Termination Reason for debugging purposes
    // std::shared_ptr< tudat::propagators::PropagationTerminationDetails > reason = dynamicsSimulator.getPropagationTerminationReason() ;
    // std::cout<<"Termination reason:"<< reason->getPropagationTerminationReason() <<"\n";

	double propagationTime = propagatedStateHistory.rbegin()->first;
    Eigen::VectorXd finalState = propagatedStateHistory.rbegin()->second;
    Eigen::VectorXd initialState = propagatedStateHistory.begin()->second;
	double integratedHeatRateFinal = finalState(6);

    double integratedFlightRangeFinal = finalState(7);


    double maxStagHeatFlux = 0;
    double stagnationHeatFlux = 0;
    double time = 0;
    double gLoad;
    double maxgLoad = 0;
    double timeMaxgLoad;
    double timeMaxStagFlux;


    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = dependentVariableHistory.begin( );
         stateIterator != dependentVariableHistory.end( ); stateIterator++ )
    {
        double time2 = stateIterator->first;
        Eigen::MatrixXd depVariables = stateIterator->second;
        stagnationHeatFlux = depVariables(5);
        std::cout<<"Stagnation Heat flux:"<<stagnationHeatFlux<<"\n";

        gLoad = depVariables(6);
        std::cout<<"G load:"<<gLoad<<"\n";


        if( maxgLoad < gLoad )
        {
            maxgLoad = gLoad;
            timeMaxgLoad = time2;
        }

        if( maxStagHeatFlux < stagnationHeatFlux )
        {
            maxStagHeatFlux = stagnationHeatFlux;
            timeMaxStagFlux = time2;
        }

        time = stateIterator->first;

    }



    if( maxStagHeatFlux > 1.0E7 || maxgLoad > 8.0 )
    {
        fitness.push_back( integratedHeatRateFinal*1000 );
        fitness.push_back( -1*integratedFlightRangeFinal*0.001 );


    }
    else
    {
    fitness.push_back( integratedHeatRateFinal );
    fitness.push_back( -1*integratedFlightRangeFinal );
    }



	//std::cout << "Initial heat load: " << integratedHeatRateInit << ", final: " << integratedHeatRateFinal << "\n";

    // Objective values OLD, before addition of constraints/penalties
   // fitness.push_back( integratedHeatRateFinal );
   // fitness.push_back( -1*integratedFlightRangeFinal );

    // Equality constraints

    // Inequality constraints
    //fitness.push_back( bodyMap_["Capsule"]->getBodyMass() - 64000 );

    return fitness;
}

std::pair<vector_double, vector_double> tudat::ShapeOptimization::get_bounds() const
{

    std::pair<vector_double, vector_double> bounds;
    bounds.first = {8.35, 2.5, 0.15, -0.35, 0.30, 0.045};
    bounds.second = {8.4, 2.9, 0.25, -0.23, 0.45, 0.05};

    return bounds;
}
