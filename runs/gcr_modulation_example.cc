#include "common/fields.hh"
#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_solarwind.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
// todo when all trajectories are updated
//#include "src/trajectory.hh"
#include "src/trajectory_parker.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{

   DataContainer container;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Set the types
//----------------------------------------------------------------------------------------------------------------------------------------------------

   using Fields = Fields<Vel_t, Mag_t, HatMag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DelVel_t>;
   using Trajectory = TrajectoryParker<Fields>;
   using Background = BackgroundSolarWind<Fields>;

   using Simulation = SimulationWorker<Trajectory>;
   using InitialTime = InitialTimeFixed<Trajectory>;
   using InitialSpace = InitialSpaceFixed<Trajectory>;
   using InitialMomentum = InitialMomentumThickShell<Trajectory>;
   using BoundarySphere = BoundarySphereAbsorb<Trajectory>;
   using BoundaryTime = BoundaryTimeExpire<Trajectory>;
   using Diffusion = DiffusionRigidityMagneticFieldPowerLaw<Trajectory>;

   using Distribution1 = DistributionSpectrumKineticEnergyPowerLaw<Trajectory>;
   using Distribution2 = DistributionTimeUniform<Trajectory>;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a simulation object
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<Simulation> simulation;
   simulation = CreateSimulation<Trajectory>(argc, argv);

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
   double umag = 4.0e7 / unit_velocity_fluid;
   GeoVector u0(umag, 0.0, 0.0);
   container.Insert(u0);

// Magnetic field
   double RS = 6.957e10 / unit_length_fluid;
   double r_ref = 3.0 * RS;
   double BmagE = 5.0e-5 / unit_magnetic_fluid;
   double Bmag_ref = BmagE * Sqr((GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid) / r_ref);
   GeoVector B0(Bmag_ref, 0.0, 0.0);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

// solar rotation vector
   double w0 = M_2PI / (25.0 * 24.0 * 3600.0) / unit_frequency_fluid;
   GeoVector Omega(0.0, 0.0, w0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction for distances closer to the Sun
   double dmax_fraction = 0.1;
   container.Insert(dmax_fraction);

   simulation->AddBackground(Background(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   double init_t = 0.0;
   container.Insert(init_t);

   simulation->AddInitial(InitialTime(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial position
   GeoVector init_pos(1.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0, 0.0);
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpace(), container);

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

   simulation->AddInitial(InitialMomentum(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Inner boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_Sun = 1;
   container.Insert(max_crossings_Sun);

// Action
   std::vector<int> actions_Sun;
   actions_Sun.push_back(-1);
   actions_Sun.push_back(-1);
   container.Insert(actions_Sun);

// Origin
   container.Insert(gv_zeros);

// Radius
   double inner_boundary = 0.01 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(inner_boundary);

   simulation->AddBoundary(BoundarySphere(), container);

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
   actions_outer.push_back(0);
   container.Insert(actions_outer);

// Origin
   container.Insert(gv_zeros);

// Radius
   double outer_boundary = 80.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(outer_boundary);

   simulation->AddBoundary(BoundarySphere(), container);

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
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Max duration of the trajectory
   double maxtime = -60.0 * 60.0 * 24.0 * 365.0 / unit_time_fluid;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTime(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Parallel mean free path
   double lam0 = 0.1 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(lam0);

// Rigidity normalization factor
   double R0 = 1.0e9 / unit_rigidity_particle;
   container.Insert(R0);

// Magnetic field normalization factor
   container.Insert(BmagE);

// Power law slope for rigidity
   double pow_law_R = 0.5;
   container.Insert(pow_law_R);

// Power law slope for magnetic field
   double pow_law_B = -1.0;
   container.Insert(pow_law_B);

// Ratio of kappa_perp to kappa_para
   double kap_rat = 0.05;
   container.Insert(kap_rat);

// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(Diffusion(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1 (spectrum)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins1(100, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(EnrKin(momentum1), 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(EnrKin(momentum2), 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(1, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro variable
   double unit_distro1 = 1.0 / (Sqr(unit_length_fluid) * unit_time_fluid * M_4PI * unit_energy_particle);
   container.Insert(unit_distro1);

// Physical units of the bin variable
   GeoVector unit_val1 = {unit_energy_particle, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

// Normalization for the "hot" boundary
   double J0 = 1.0 / unit_distro1;
   container.Insert(J0);

// Characteristic energy
   double T0 = 1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0);

// Spectral power law
   double pow_law_T = -1.8;
   container.Insert(pow_law_T);

// Constant value for the "cold" condition
   double val_cold1 = 0.0;
   container.Insert(val_cold1);

   simulation->AddDistribution(Distribution1(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2 (exit time)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins2(100, 0, 0);
   container.Insert(n_bins2);
   
// Smallest value
   GeoVector minval2(0.0, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(maxtime, 0.0, 0.0);
   container.Insert(maxval2);

// Linear or logarithmic bins
   MultiIndex log_bins2(0, 0, 0);
   container.Insert(log_bins2);

// Add outlying events to the end bins
   MultiIndex bin_outside2(1, 0, 0);
   container.Insert(bin_outside2);

// Physical units of the distro variable
   double unit_distro2 = 1.0;
   container.Insert(unit_distro2);

// Physical units of the bin variable
   GeoVector unit_val2 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val2);

// Keep records
   bool keep_records2 = false;
   container.Insert(keep_records2);

// Value for the "hot" condition
   double val_hot2 = 1.0;
   container.Insert(val_hot2);

// Value for the "cold" condition
   double val_cold2 = 0.0;
   container.Insert(val_cold2);

// Value for which time to bin
   int val_time2 = 1;
   container.Insert(val_time2);

   simulation->AddDistribution(Distribution2(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "gcr_modulation_example_distro_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, "modulated_gcr_spectrum.dat", true);
   simulation->PrintDistro1D(1, 0, "time_gcr_trajec.dat", true);
   
   return 0;
};
