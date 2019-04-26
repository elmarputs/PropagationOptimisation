#include <iostream>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/algorithms/sea.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/nlopt.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "pagmo/algorithms/xnes.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/ihs.hpp"

pagmo::algorithm getMultiObjectiveAlgorithm( const int index )
{
    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::nsga2( ) };
        return algo;
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::moead( ) };
        return algo;
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        return algo;
        break;
    }
    // Change the crossover rate
    case 3:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.66, 10., 0.01,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 4:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.33, 10., 0.01,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 5:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.0, 10., 0.01,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    //Change the mutation rate
    case 6:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.05,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 7:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.1,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 8:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.5,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    // Change the crossover distribution index
    case 9:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 1, 0.01,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 10:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 50, 0.01,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 11:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 90., 0.01,  50.,pagmo::random_device::next()) };
        return algo;
        break;
    }
        // Change the mutation distribution index
    case 12:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.01,  1., pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 13:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.01,  10., pagmo::random_device::next()) };
        return algo;
        break;
    }
    case 14:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.01,  90., pagmo::random_device::next()) };
        return algo;
        break;
    }
        // Optimal Algorithm settings
    case 15:
    case 16:
    {
        pagmo::algorithm algo{ pagmo::nsga2( 1u, 0.95, 10., 0.01,  50., pagmo::random_device::next()) };
        return algo;
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, multi-objective pagmo algorithm " + std::to_string( index ) + " was not found." );
    }
    }
}
pagmo::algorithm getAlgorithm( const int index )
{
    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::de( ) };
        return algo;
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::sade( ) };
        return algo;
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::de1220( ) };
        return algo;
        break;
    }
    case 3:
    {
        pagmo::algorithm algo{ pagmo::pso( ) };
        return algo;
        break;
    }
    case 4:
    {
        pagmo::algorithm algo{ pagmo::sea( ) };
        return algo;
        break;
    }
    case 5:
    {
        pagmo::algorithm algo{ pagmo::sga( ) };
        return algo;
        break;
    }
    case 6:
    {
        pagmo::algorithm algo{ pagmo::simulated_annealing( ) };
        return algo;
        break;
    }
    case 7:
    {
        pagmo::algorithm algo{ pagmo::bee_colony( ) };
        return algo;
        break;
    }
    case 8:
    {
        pagmo::algorithm algo{ pagmo::cmaes( ) };
        return algo;
        break;
    }
    case 9:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        return algo;
        break;
    }
    case 10:
    {
        pagmo::algorithm algo{ pagmo::xnes( ) };
        return algo;
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, sinlge-objective pagmo algorithm " + std::to_string( index ) + " was not found." );
    }
    }
}

