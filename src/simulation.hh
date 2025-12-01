/*!
\file simulation.hh
\brief Declares application level methods to perform a complete simulation from start to finish
\author Juan G Alonso Guzman
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMULATION_HH
#define SPECTRUM_SIMULATION_HH

#include "simulation_master.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Top level methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 01/07/2025
\param[in] argc Number of command line arguments
\param[in] argv Command line arguments
*/
template <typename HConfig, typename Trajectory>
std::unique_ptr<SimulationWorker<HConfig, Trajectory>> CreateSimulation(int argc, char** argv)
{
   using MPI = MPI<HConfig, HConfig::MPI_enabled>;
// Initialize a single instance of "MPI" that will persist until the program terminates
   static MPI mpicfg(argc, argv);

   if (MPI::is_master) return std::make_unique<SimulationMaster<HConfig, Trajectory>>();
   else if (MPI::is_server) return std::make_unique<SimulationServer<HConfig, Trajectory>>();
   else return std::make_unique<SimulationWorker<HConfig, Trajectory>>();
};

};


#endif
