#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_smooth_shock.hh"
#include "src/diffusion_other.hh"
#include "src/source_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_momentum.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i;
   DataContainer container;

//--------------------------------------------------------------------------------------------------
// Create a simulation object
//--------------------------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//--------------------------------------------------------------------------------------------------
// Particle type
//--------------------------------------------------------------------------------------------------

   int specie = SPECIES_PROTON_BEAM;
   simulation->SetSpecie(specie);

//--------------------------------------------------------------------------------------------------
// Background
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   container.Insert(0.0);

// Origin
   container.Insert(gv_zeros);

// Upstream velocity
   double U_up = 4.0e7 / unit_velocity_fluid;
   GeoVector u0(U_up, 0.0, 0.0);
   container.Insert(u0);

// Upstream magnetic field
   double B_up = 5.0e-7 / unit_magnetic_fluid;
   GeoVector B0(B_up, 0.0, 0.0);
   container.Insert(B0);

// Maximum displacement
   double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   double dmax = 0.1 * one_au;
   container.Insert(dmax);

// Shock normal
   GeoVector n_shock (-1.0, 0.0, 0.0);
   container.Insert(n_shock);

// Shock velocity
   double v_shock = 0.0;
   container.Insert(v_shock);

// Compression ratio
   double s = 4.0;
   container.Insert(s);
   double U_dn = U_up / s;

// Shock width
   double w_sh = 0.01 * one_au;
   container.Insert(w_sh);

// dmax fraction
   double dmax_fraction = 0.1;
   container.Insert(dmax_fraction);

   simulation->AddBackground(BackgroundSmoothShock(), container);

//--------------------------------------------------------------------------------------------------
// Time initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Time to obtain solution
   double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;
   double t_final = 100.0 * one_day;
   container.Insert(t_final);

   simulation->AddInitial(InitialTimeFixed(), container);

//--------------------------------------------------------------------------------------------------
// Spatial initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Point to obtain solution
   double z_spectrum = 0.1 * one_au;
   GeoVector init_pos(z_spectrum, 0.0, 0.0);
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpaceFixed(), container);

//--------------------------------------------------------------------------------------------------
// Momentum initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Lower bound for momentum
   double one_MeV = SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle;
   double p0 = Mom(1.0 * one_MeV, specie);
   container.Insert(p0);

// Upper bound for momentum
   double pf = Mom(100.0 * one_MeV, specie);
   container.Insert(pf);

// Log bias
   bool log_bias = true;
   container.Insert(log_bias);

   simulation->AddInitial(InitialMomentumThickShell(), container);

//--------------------------------------------------------------------------------------------------
// Time limit
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Maximum crossings
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Actions vector
   std::vector<int> actions_time;
   actions_time.push_back(1);
   container.Insert(actions_time);

// Initial time at which to stop integrating
   container.Insert(0.0);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//--------------------------------------------------------------------------------------------------
// Diffusion model
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Reference diffusion coefficient
   double kappa_up = 1.0e20 / unit_diffusion_fluid;
   double kappa_dn = kappa_up * Sqr(U_dn / U_up);
   container.Insert(kappa_up);

// Normalization of bulk velocity
   container.Insert(U_up);

// Power of bulk velocity dependance
   double power_law_U = 2.0;
   container.Insert(power_law_U);

// Normalization of particle momentum
   double mom_0 =  0.1 * SpeciesMasses[specie] * c_code;
   container.Insert(mom_0);

// Power of particle momentum dependance
   double power_law_p = 0.0;
   container.Insert(power_law_p);

// Ratio of perpendicular to parallel diffusion
   double kap_rat = 0.0;
   container.Insert(kap_rat);

   simulation->AddDiffusion(DiffusionFlowMomentumPowerLaw(), container);

//--------------------------------------------------------------------------------------------------
// Source term
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Injection momentum
   container.Insert(p0);

// Rate at shock
   container.Insert(1.0);

   simulation->AddSource(SourceMomentumInjection(), container);

//--------------------------------------------------------------------------------------------------
// Distribution (momentum)
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins3(100, 0, 0);
   container.Insert(n_bins3);

// Smallest value
   GeoVector minval3(p0, 0.0, 0.0);
   container.Insert(minval3);

// Largest value
   GeoVector maxval3(pf, 0.0, 0.0);
   container.Insert(maxval3);

// Linear or logarithmic bins
   MultiIndex log_bins3(1, 0, 0);
   container.Insert(log_bins3);

// Add outlying events to the end bins
   MultiIndex bin_outside3(0, 0, 0);
   container.Insert(bin_outside3);

// Physical units of the distro variable
   double unit_distro3 = 1.0;
   container.Insert(unit_distro3);

// Physical units of the bin variable which is momentum here. This is for x axis.
   GeoVector unit_val3 = {unit_momentum_particle, 1.0, 1.0};
   container.Insert(unit_val3);

// Don't keep records
   bool keep_records3 = false;
   container.Insert(keep_records3);

// Constant value for the "hot" condition
   double val_hot3 = 1.0;
   container.Insert(val_hot3);

// Constant value for the "cold" condition
   double val_cold3 = 0.0;
   container.Insert(val_cold3);

// Which coordinates to use for value: 0 initial, 1 final
   int val_time3 = 0;
   container.Insert(val_time3);

// Which coordinate representation to use for value: 0 "native coordinates", 1 locally spherical with B || z
   int val_coord3 = 0;
   container.Insert(val_coord3);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//--------------------------------------------------------------------------------------------------
// Run the simulation
//--------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "main_test_diff_shock_acc_" + simulation->GetTrajectoryName() + "_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "spectrum.dat", false);

   if(MPI_Config::is_master) {
      std::cout << std::endl;
      std::cout << "DIFFUSIVE SHOCK ACCELERATION" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Trajectory type: " << simulation->GetTrajectoryName() << std::endl;
      std::cout << "z_diff = " << kappa_up / U_up / one_au << " au" << std::endl;
      std::cout << "z_spectrum = " << z_spectrum / one_au << " au" << std::endl;
      std::cout << "w_shock = " << w_sh / one_au << " au" << std::endl;
      std::cout << "t_acc = " << 4.0 * kappa_up / Sqr(U_up) / one_day << " days" << std::endl;
      std::cout << "t_final = " << t_final / one_day << " days" << std::endl;
      std::cout << "min dt_adv = " << w_sh / (U_up + (kappa_up - kappa_dn) / w_sh) / one_day << " days" << std::endl;
      std::cout << "max dt_dif = " << Sqr(w_sh) / kappa_dn << " days" << std::endl;
      std::cout << "++++++++++++++++++++" << std::endl;
      std::cout << "Distribution files outputed to " << simulation_files_prefix << std::endl;
      std::cout << std::endl;
   };

   return 0;
};
