#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_server_cartesian.hh"
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
   container.Insert(gv_zeros);

// Effective "mesh" resolution
   double dmax = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

   std::string fname_pattern = "../cartesian_backgrounds/parker_25_25_25";
   simulation->AddBackground(BackgroundServerCartesian(), container, fname_pattern);

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

   GeoVector init_pos(1.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0, 0.0);
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpaceFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Lower bound for momentum
   double momentum1 = Mom(10.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum1);

// Upper bound for momentum
   double momentum2 = Mom(5000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum2);

// Log bias
   bool log_bias = true;
   container.Insert(log_bias);

   simulation->AddInitial(InitialMomentumThickShell(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Inner boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_Sun = 1;
   container.Insert(max_crossings_Sun);

// Action
   std::vector<int> actions_Sun;
   actions_Sun.push_back(1);
   container.Insert(actions_Sun);

// Origin
   container.Insert(gv_zeros);

// Radius
   double inner_boundary = 0.05 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(inner_boundary);

   simulation->AddBoundary(BoundarySphereAbsorb(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Outer boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_outer = 1;
   container.Insert(max_crossings_outer);

// Action
   std::vector<int> actions_outer;
   actions_outer.push_back(0);
   container.Insert(actions_outer);

// Origin
   container.Insert(gv_zeros);

// Radius
   double outer_boundary = 80.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(outer_boundary);

   simulation->AddBoundary(BoundarySphereAbsorb(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time limit
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Max duration of the trajectory
   double maxtime = -60.0 * 60.0 * 24.0 * 365.0 / unit_time_fluid;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Diffusion coefficient normalization factor
   double kap0 = 1.5e22 / unit_diffusion_fluid;
   container.Insert(kap0);

// Rigidity normalization factor
   double T0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0);

// Magnetic field normalization factor
   double r0 = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(r0);

// Power law slope for rigidity
   double pow_law_T = 1.0;
   container.Insert(pow_law_T);

// Power law slope for magnetic field
   double pow_law_r = 2.0;
   container.Insert(pow_law_r);

// Ratio of kappa_perp to kappa_para
   double kap_rat = 1.00;
   container.Insert(kap_rat);

// Stream dependance index
   int stream_dep_idx = 0;
   container.Insert(stream_dep_idx);

// Upstream flow speed (unused)
   int u_up = 1.0;
   container.Insert(u_up);

// Shock width (unused)
   int w_sh = 0.0;
   container.Insert(w_sh);

// Shock strength (unused)
   int s_sh = 1.0;
   container.Insert(s_sh);

// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(DiffusionKineticEnergyRadialDistancePowerLaw(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins(100, 0, 0);
   container.Insert(n_bins);
   
// Smallest value
   GeoVector minval(EnrKin(momentum1, specie), 0.0, 0.0);
   container.Insert(minval);

// Largest value
   GeoVector maxval(EnrKin(momentum2, specie), 0.0, 0.0);
   container.Insert(maxval);

// Linear or logarithmic bins
   MultiIndex log_bins(1, 0, 0);
   container.Insert(log_bins);

// Add outlying events to the end bins
   MultiIndex bin_outside(0, 0, 0);
   container.Insert(bin_outside);

// Physical units of the distro variable
   double unit_distro = 1.0 / (Sqr(unit_length_fluid) * unit_time_fluid * M_4PI * unit_energy_particle);
   container.Insert(unit_distro);

// Physical units of the bin variable
   GeoVector unit_val = {unit_energy_particle, 1.0, 1.0};
   container.Insert(unit_val);

// Don't keep records
   bool keep_records = false;
   container.Insert(keep_records);

//! Normalization for the "hot" boundary
   double J0 = 1.0 / unit_distro;
   container.Insert(J0);

//! Characteristic energy
   T0 = 1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0);

//! Spectral power law
   pow_law_T = -1.8;
   container.Insert(pow_law_T);

//! Constant value for the "cold" condition
   double val_cold = 0.0;
   container.Insert(val_cold);

   simulation->AddDistribution(DistributionSpectrumKineticEnergyPowerLaw(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "main_test_modulation_cartesian_" + simulation->GetTrajectoryName() + "_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "spectrum.dat", true);

   if(MPI_Config::is_master) {
      std::cout << std::endl;
      std::cout << "MODULATION CARTESIAN PARKER FIELD" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Trajectory type: " << simulation->GetTrajectoryName() << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Distribution files outputed to " << simulation_files_prefix << std::endl;
      std::cout << std::endl;
   };
   
   return 0;
};
