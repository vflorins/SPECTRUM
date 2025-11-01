#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_uniform.hh"
#include "src/boundary_time.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{

   DataContainer container;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a simulation object
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle type
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int specie = SPECIES_PROTON_BEAM;
   simulation->SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   double t0 = 0.0;
   container.Insert(t0);

// Origin
   container.Insert(gv_zeros);

// Velocity
   container.Insert(gv_zeros);

// Magnetic field
   double Bmag = 1.0 / unit_magnetic_fluid;
   GeoVector B0(0.0, 0.0, Bmag);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax = 0.1;
   container.Insert(dmax);

   simulation->AddBackground(BackgroundUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Start time for interval
   double start_t = -5.0 / unit_time_fluid;
   container.Insert(start_t);

// End time for interval
   double end_t = 5.0 / unit_time_fluid;
   container.Insert(end_t);

// Number of subintervals (0 for random points)
   int n_intervals = 0;
   container.Insert(n_intervals);

   simulation->AddInitial(InitialTimeInterval(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   container.Insert(gv_zeros);

   double radius = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(radius);

   simulation->AddInitial(InitialSpaceSphere(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double momentum = Mom(100.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum);

   double theta = DegToRad(90.0);
   container.Insert(theta);

   simulation->AddInitial(InitialMomentumRing(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition (terminal)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(0);
   actions_time.push_back(0);
   container.Insert(actions_time);
   
// Duration of the trajectory
   double maxtime = 10.0 / unit_time_fluid;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1 (time)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins1(100, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(start_t, 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(end_t, 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(0, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro variable
   double unit_distro1 = 1.0;
   container.Insert(unit_distro1);

// Physical units of the bin variable
   GeoVector unit_val1 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

// Value for the "hot" condition
   double val_hot1 = 1.0;
   container.Insert(val_hot1);

// Value for the "cold" condition
   double val_cold1 = 0.0;
   container.Insert(val_cold1);

// Coordinates to use (initial or final)
   int val_time1 = 0;
   container.Insert(val_time1);

   simulation->AddDistribution(DistributionTimeUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2 (position)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins2(0, 0, 0);
   container.Insert(n_bins2);
   
// Smallest value
   GeoVector minval2(0.0, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(0.0, 0.0, 0.0);
   container.Insert(maxval2);

// Linear or logarithmic bins
   MultiIndex log_bins2(0, 0, 0);
   container.Insert(log_bins2);

// Add outlying events to the end bins
   MultiIndex bin_outside2(0, 0, 0);
   container.Insert(bin_outside2);

// Physical units of the distro variable
   double unit_distro2 = 1.0;
   container.Insert(unit_distro2);

// Physical units of the bin variable
   GeoVector unit_val2 = {unit_length_fluid, unit_length_fluid, unit_length_fluid};
   container.Insert(unit_val2);

// Don't keep records
   bool keep_records2 = true;
   container.Insert(keep_records2);

// Value for the "hot" condition
   double val_hot2 = 1.0;
   container.Insert(val_hot2);

// Value for the "cold" condition
   double val_cold2 = 0.0;
   container.Insert(val_cold2);

// Coordinates to use (initial or final)
   int val_time2 = 0;
   container.Insert(val_time2);

// Coordinate representation to use
   int val_coord2 = 0;
   container.Insert(val_coord2);

   simulation->AddDistribution(DistributionPositionUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "main_test_init_cond_records_" + simulation->GetTrajectoryName() + "_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "init_time.dat", true);
   simulation->PrintRecords(1, simulation_files_prefix + "pos_records.dat", false);

   if(MPI_Config::is_master) {
      std::cout << std::endl;
      std::cout << "INITIAL CONDITION RECORDS" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Trajectory type: " << simulation->GetTrajectoryName() << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Distribution files outputed to " << simulation_files_prefix << std::endl;
      std::cout << std::endl;
   };
   
   return 0;
};
