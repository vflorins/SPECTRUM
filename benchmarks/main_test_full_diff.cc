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
   double momentum = Mom(100.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum);

   simulation->AddInitial(InitialMomentumShell(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Perpendicular diffusion coefficient
   double D0 = Sqr(LarmorRadius(momentum,Bmag,specie)) * CyclotronFrequency(Vel(momentum),Bmag,specie);
   container.Insert(D0);

// Parallel diffusion coefficient
   container.Insert(D0);

// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(DiffusionFullConstant(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition 1 (recurrent)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   int max_crossings_time = -1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(0);
   actions_time.push_back(0);
   container.Insert(actions_time);
   
// Spacing between dumps
   double timemark = 0.1 * D0 / Sqr(Vel(momentum));
   container.Insert(timemark);

   simulation->AddBoundary(BoundaryTimeRecurrent(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition 2 (cutoff)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1
   max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   actions_time.clear();
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Duration of the trajectory
   double maxtime = 1001.0 * timemark;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins1(maxtime / timemark - 1, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(0.5 * timemark, 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1((n_bins1[0] + 0.5) * timemark, 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(0, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro variable
   GeoVector unit_distro1 = unit_length_fluid * gv_ones;
   container.Insert(unit_distro1);

// Physical units of the bin variable
   GeoVector unit_val1 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

   simulation->AddDistribution(DistributionPositionCumulativeOrder1(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   container.Insert(n_bins1);
   
// Smallest value
   container.Insert(minval1);

// Largest value
   container.Insert(maxval1);

// Linear or logarithmic bins
   container.Insert(log_bins1);

// Add outlying events to the end bins
   container.Insert(bin_outside1);

// Physical units of the distro variable
   GeoMatrix unit_distro2 = Sqr(unit_length_fluid) * gm_ones;
   container.Insert(unit_distro2);

// Physical units of the bin variable
   container.Insert(unit_val1);

// Don't keep records
   container.Insert(keep_records1);

   simulation->AddDistribution(DistributionPositionCumulativeOrder2(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "main_test_full_diff_" + simulation->GetTrajectoryName() + "_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "cumulative_distro1.dat", true);
   simulation->PrintDistro1D(1, 0, simulation_files_prefix + "cumulative_distro2.dat", true);

   if(MPI_Config::is_master) {
      std::cout << std::endl;
      std::cout << "FULL DIFFUSION" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Trajectory type: " << simulation->GetTrajectoryName() << std::endl;
      std::cout << "D0 = " << D0 * Sqr(unit_length_fluid) / unit_time_fluid << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Distribution files outputed to " << simulation_files_prefix << std::endl;
      std::cout << std::endl;
   };
   
   return 0;
};
