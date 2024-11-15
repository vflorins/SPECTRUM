#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_smooth_shock.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/boundary_momentum.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{
   BackgroundSmoothShock background;
   SpatialData spdata;
   int i,j,k;

   spdata._mask = BACKGROUND_ALL;
   DataContainer container;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a simulation object
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle type
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int specie = Specie::proton;
   simulation->SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();
   
// Origin_t0
   double t0 = 0.0;
   container.Insert(t0);

// Origin_r0
   container.Insert(gv_zeros);
   
// Velocity
   double U_up = 4.5e7 / unit_velocity_fluid;
   GeoVector u0(U_up, 0.0, 0.0);
   container.Insert(u0);

//! Reference magnetic field (persistent) in Gauss
   double B_up = 5.0e-5 / unit_magnetic_fluid;
   GeoVector B0(B_up, 0.0, 0.0);
   container.Insert(B0);

//! time reference value (persistent)
   double U_dn = 1.1e7 / unit_velocity_fluid;   ///s=4
//    double U_dn = 1.5e6 / unit_velocity_fluid;   ///s=3
//   double U_dn = 2.25e6 / unit_velocity_fluid;   ///s=2
//   double U_dn = 4.5e6 / unit_velocity_fluid;   ///s=1
   double dt_max = 1.0e4 / unit_time_fluid;
   double dmax0 = dt_max * U_dn;
   container.Insert(dmax0);

//! Shock starting position (persistent) r_shock
   container.Insert(gv_zeros);

//! Shock normal (persistent)
   GeoVector n_shock (-1.0, 0.0, 0.0);
   container.Insert(n_shock);

//! Shock velocity (persistent)
   double v_shock = 0.0;
   container.Insert(v_shock);

//! Downstream flow vector (persistent), "u0" is upstream flow vector divide by shock strength
   GeoVector u1 (U_dn, 0.0, 0.0);
   container.Insert(u1);

//! Downstream magnetic field (persistent), "B0" is upstream magnetic field
   double B_dn = 2.0e-4 / unit_magnetic_fluid;    ///s=4
//   double B_dn = 1.5e-4 / unit_magnetic_fluid;    ///s=3
//   double B_dn = 1.0e-4 / unit_magnetic_fluid;    ///s=2
//   double B_dn = 5.0e-5 / unit_magnetic_fluid;    ///s=1
//   GeoVector B1 (B_up, 0.0, 0.0);
   container.Insert(B_dn);

//! Width of shock transition region (persistent)
//   double width_shock =  0.0000001 * dt_max * U_dn;
//   double width_shock =  2.0e4 * dt_max * U_dn;    ///This gives no hump with kd=e-4 times that.
   double width_shock = 1.4e-2 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(width_shock);

   simulation->AddBackground(BackgroundSmoothShock(), container);
   
//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// initial time for interval
   double inittime = 0.0 / unit_time_fluid;
   container.Insert(inittime);

   simulation->AddInitial(InitialTimeFixed(), container);

 //----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------
   container.Clear();  
// Initial position which is actually the final position of the particle in backward simulation
   GeoVector init_pos(0.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0, 0.0);
   container.Insert(init_pos);
   
   // Initial position which is actually the final position of the particle in backward simulation
   GeoVector final_pos( 10.0 * width_shock * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0, 0.0);
   container.Insert(final_pos);
   
   int n_interval = 0;
   container.Insert(n_interval);
    simulation->AddInitial(InitialSpaceLine(), container);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------
//The initial momentum distribution is a linear biased “thick shell”, which produces random, uniformly
//distributed pitch and gyrophase angles and a linear distributed momentum magnitude between 1 MeV and 300 MeV.

   container.Clear();

// Lower bound for momentum
   double momentum1 = Mom(1.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum1);

// Upper bound for momentum
   double momentum2 = Mom(20.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum2);

// Log bias
   bool log_bias = false;
   container.Insert(log_bias);

   simulation->AddInitial(InitialMomentumThickShell(), container);

//---------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Momentum injection boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

//// here we are binning the weight=0(i.e. Action=1 i.e, Uniformcold() with weight=Value_cold) for the particles which reach the momentum inner boundary. This includes  
//// their count in the total number of particle counts. action_1_in is for the first distribution i.e. position distribution and action_2 is for the second distribution 
//// i.e. momentum distribution. For momentum distribution, we set the weight=1(i.e. Action=0 i.e, UniformHot() with weight=Value_Hot). 
   container.Clear();

// Max crossings
   int max_crossings_Sun = 1;
   container.Insert(max_crossings_Sun);

// Action
   std::vector<int> actions_Sun;
   actions_Sun.clear();
   actions_Sun.push_back(0);
   actions_Sun.push_back(0);
   container.Insert(actions_Sun);

// Radius

   double momentum_inner = Mom(1.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum_inner);
   

   simulation->AddBoundary(BoundaryMomentumInject(), container);
//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time limit
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1. Total number of crossings allowed, -1 means unlimited (persistent)
   int max_crossings_time = -1;
   container.Insert(max_crossings_time);

// Action. Actions to perform when recording events, -1 for no action (discard completely)
   std::vector<int> actions_time;
   actions_time.push_back(1);
   actions_time.push_back(1);
   container.Insert(actions_time);
   
// Max duration of the trajectory
//   double maxtime = - 5.0e6 / unit_time_fluid;          ///////t1 dt_phy=0.05*dt
   double maxtime = - 58.0 * 24.0 * 60 * 60 / unit_time_fluid;          ///////58=5e6, 30=2.5e6, 5.7=5e5
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);
   

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------
   container.Clear();

// Reference diffusion coefficient (normalized units)
   double kappa0 = 2.4e19 / unit_length_fluid / unit_velocity_fluid;     ///This value is equal to 1.4e-4 * U_up / 2.0 with dt_max = 1.0/ unit_length_fluid
   container.Insert(kappa0);

// Normalization of velocity
   double U_0 = U_dn;
   container.Insert(U_0);

// Power of flow velocity
   double power_law_U = 2.0;
   container.Insert(power_law_U);

// Ratio of perp to para diffusion
   double kap_rat = 0.0;
   container.Insert(kap_rat);
   
// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(DiffusionFlowPowerLaw(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------
   container.Clear();

// Reference diffusion coefficient (normalized units)
//   double kappa0 = width_shock * U_dn / 2.0 / dt_max;
//   double D0 = 2.4e19 / unit_length_fluid / unit_velocity_fluid; 
//   container.Insert(D0);
   
// Pass ownership of "diffusion" to simulation
//   simulation->AddDiffusion(DiffusionParaConstant(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1 (space)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins1(200, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(-1.0, 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(2.0, 0.0, 0.0);
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
   GeoVector unit_val1 = {1.0, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

//! Constant value for the "hot" condition (persistent)
   double val_hot1 = 1.0;
   container.Insert(val_hot1);
   
//! Constant value for the "cold" condition (persistent)
   double val_cold1 = 0.0;
   container.Insert(val_cold1);
   
//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time1 = 0;
   container.Insert(val_time1);

//! Which coordinate representation to use for value: 0 "native coordinates", 1 locally spherical with B || z (persisent)
   int val_coord1 = 0;
   container.Insert(val_coord1);
   
//  GeoVector momentum_1(5.5e6, 0.0, 0.0) ;
//   GeoVector momentum_1(Mom(0.9 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie), 0.0, 0.0);
//   container.Insert(momentum_1);

//  GeoVector momentum_2(8.5e8, 0.0, 0.0) ;
//   GeoVector momentum_2(Mom(5.1 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie), 0.0, 0.0);
//   container.Insert(momentum_2);


   simulation->AddDistribution(DistributionPositionUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2 (momentum)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins2(20, 0, 0);
   container.Insert(n_bins2);
   
// Smallest value
   GeoVector minval2(momentum1, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(momentum2, 0.0, 0.0);
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
   GeoVector unit_val2 = {1.0, 1.0, 1.0};
   container.Insert(unit_val2);

// Don't keep records
   bool keep_records2 = false;
   container.Insert(keep_records2);

//! Constant value for the "hot" condition (persistent)
   double val_hot2 = 1.0;
   container.Insert(val_hot2);
   
//! Constant value for the "cold" condition (persistent)
   double val_cold2 = 0.0;
   container.Insert(val_cold2);
   
//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time2 = 0;
   container.Insert(val_time2);

//! Which coordinate representation to use for value: 0 "native coordinates", 1 locally spherical with B || z (persisent)
   int val_coord2 = 0;
   container.Insert(val_coord2);
   
 // Position minimum
//  GeoVector position_1(-0.014, 0.0, 0.0);
//  container.Insert(position_1);

// Position maximum
//  GeoVector position_2(1.0, 0.0, 0.0);
//  container.Insert(position_2);

   simulation->AddDistribution(DistributionMomentumUniform(), container);
   
//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   
  
   std::string simulation_files_prefix = "back_sim/smooth_shock/6nov_uniform_distri/p1-20/distro_";

   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "pos_5e6_p1-50_x0-10xsh.dat", false); 
   simulation->PrintDistro1D(1, 0, simulation_files_prefix + "mom_5e6_p1-50_x0-10xsh.dat", false); 
   
   return 0;
};


