#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_uniform.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
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

// Initial time
   double init_t = 0.0;
   container.Insert(init_t);

   simulation->AddInitial(InitialTimeFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   container.Insert(gv_zeros);

   simulation->AddInitial(InitialSpaceFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double momentum = Mom(1.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum);

   double theta = DegToRad(45.0);
   container.Insert(theta);

   simulation->AddInitial(InitialMomentumRing(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Scattering frequency
   double D0 = 0.01 * CyclotronFrequency(Vel(momentum),Bmag,specie);
   container.Insert(D0);

// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(DiffusionIsotropicConstant(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition 1 (pass)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   int max_crossings_time = -1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(-1);
   actions_time.push_back(0);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// First timemark
   double timemark1 = 0.5 / D0;
   container.Insert(timemark1);

   simulation->AddBoundary(BoundaryTimePass(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition 2 (pass)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   container.Insert(max_crossings_time);

// Action
   actions_time.clear();
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(0);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Second timemark
   double timemark2 = 2.0 * timemark1;
   container.Insert(timemark2);

   simulation->AddBoundary(BoundaryTimePass(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition 3 (pass)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   container.Insert(max_crossings_time);

// Action
   actions_time.clear();
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(0);
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Third timemark
   double timemark3 = 2.0 * timemark2;
   container.Insert(timemark3);

   simulation->AddBoundary(BoundaryTimePass(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition 4 (cutoff)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1
   max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   actions_time.clear();
   actions_time.push_back(0);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(0);
   container.Insert(actions_time);
   
// Duration of the trajectory
   double maxtime = 2.0 * timemark3;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins(1, 100, 0);
   container.Insert(n_bins);
   
// Smallest value
   GeoVector minval(0.99 * momentum, -1.0, 0.0);
   container.Insert(minval);

// Largest value
   GeoVector maxval(1.01 * momentum, 1.0, 0.0);
   container.Insert(maxval);

// Linear or logarithmic bins
   MultiIndex log_bins(0, 0, 0);
   container.Insert(log_bins);

// Add outlying events to the end bins
   MultiIndex bin_outside(0, 0, 1);
   container.Insert(bin_outside);

// Physical units of the distro variable
   double unit_distro = 1.0;
   container.Insert(unit_distro);

// Physical units of the bin variable
   GeoVector unit_val = {unit_momentum_particle, 1.0, 1.0};
   container.Insert(unit_val);

// Don't keep records
   bool keep_records = false;
   container.Insert(keep_records);

// Value for the "hot" condition
   double val_hot = 1.0;
   container.Insert(val_hot);

// Value for the "cold" condition
   double val_cold = 0.0;
   container.Insert(val_cold);

// Coordinates to use (initial or final)
   int val_time = 0;
   container.Insert(val_time);

// Coordinate representation to use (native or locally spherical)
   int val_coord = 1;
   container.Insert(val_coord);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   container.Insert(n_bins);
   
// Smallest value
   container.Insert(minval);

// Largest value
   container.Insert(maxval);

// Linear or logarithmic bins
   container.Insert(log_bins);

// Add outlying events to the end bins
   container.Insert(bin_outside);

// Physical units of the distro variable
   container.Insert(unit_distro);

// Physical units of the bin variable
   container.Insert(unit_val);

// Don't keep records
   container.Insert(keep_records);

// Value for the "hot" condition
   container.Insert(val_hot);

// Value for the "cold" condition
   container.Insert(val_cold);

// Coordinates to use (initial or final)
   val_time = 1;
   container.Insert(val_time);

// Coordinate representation to use (native or locally spherical)
   container.Insert(val_coord);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 3
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   container.Insert(n_bins);
   
// Smallest value
   container.Insert(minval);

// Largest value
   container.Insert(maxval);

// Linear or logarithmic bins
   container.Insert(log_bins);

// Add outlying events to the end bins
   container.Insert(bin_outside);

// Physical units of the distro variable
   container.Insert(unit_distro);

// Physical units of the bin variable
   container.Insert(unit_val);

// Don't keep records
   container.Insert(keep_records);

// Value for the "hot" condition
   container.Insert(val_hot);

// Value for the "cold" condition
   container.Insert(val_cold);

// Coordinates to use (initial or final)
   container.Insert(val_time);

// Coordinate representation to use (native or locally spherical)
   container.Insert(val_coord);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 4
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   container.Insert(n_bins);
   
// Smallest value
   container.Insert(minval);

// Largest value
   container.Insert(maxval);

// Linear or logarithmic bins
   container.Insert(log_bins);

// Add outlying events to the end bins
   container.Insert(bin_outside);

// Physical units of the distro variable
   container.Insert(unit_distro);

// Physical units of the bin variable
   container.Insert(unit_val);

// Don't keep records
   container.Insert(keep_records);

// Value for the "hot" condition
   container.Insert(val_hot);

// Value for the "cold" condition
   container.Insert(val_cold);

// Coordinates to use (initial or final)
   container.Insert(val_time);

// Coordinate representation to use (native or locally spherical)
   container.Insert(val_coord);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 5
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   container.Insert(n_bins);
   
// Smallest value
   container.Insert(minval);

// Largest value
   container.Insert(maxval);

// Linear or logarithmic bins
   container.Insert(log_bins);

// Add outlying events to the end bins
   container.Insert(bin_outside);

// Physical units of the distro variable
   container.Insert(unit_distro);

// Physical units of the bin variable
   container.Insert(unit_val);

// Don't keep records
   container.Insert(keep_records);

// Value for the "hot" condition
   container.Insert(val_hot);

// Value for the "cold" condition
   container.Insert(val_cold);

// Coordinates to use (initial or final)
   container.Insert(val_time);

// Coordinate representation to use (native or locally spherical)
   container.Insert(val_coord);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "main_test_pa_distro_isotrop_" + simulation->GetTrajectoryName() + "_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 1, simulation_files_prefix + "initial_distro.dat", true);
   simulation->PrintDistro1D(1, 1, simulation_files_prefix + "timemark1_distro.dat", true);
   simulation->PrintDistro1D(2, 1, simulation_files_prefix + "timemark2_distro.dat", true);
   simulation->PrintDistro1D(3, 1, simulation_files_prefix + "timemark3_distro.dat", true);
   simulation->PrintDistro1D(4, 1, simulation_files_prefix + "final_distro.dat", true);

   if(MPI_Config::is_master) {
      std::cout << std::endl;
      std::cout << "PITCH ANGLE DISTRIBUTION ISOTROPIZATION" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Trajectory type: " << simulation->GetTrajectoryName() << std::endl;
      std::cout << "D0 = " << D0 / unit_time_fluid << std::endl;
      std::cout << "distros outputed at times:" << std::endl;
      std::cout << "\t initial   = " << 0.0 << " s" << std::endl;
      std::cout << "\t timemark1 = " << timemark1 << " s" << std::endl;
      std::cout << "\t timemark2 = " << timemark2 << " s" << std::endl;
      std::cout << "\t timemark3 = " << timemark3 << " s" << std::endl;
      std::cout << "\t final     = " << maxtime << " s" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Distribution files outputed to " << simulation_files_prefix << std::endl;
      std::cout << std::endl;
   };

   return 0;
};
