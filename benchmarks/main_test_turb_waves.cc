#include "src/background_waves.hh"
#include "src/boundary_time.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include "src/traj_config.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{

   DataContainer container;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<TrajectoryBase> trajectory = std::make_unique<TrajectoryType>();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Connect RNG
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::shared_ptr<RNG> rng = std::make_shared<RNG>(time(NULL));
   trajectory->ConnectRNG(rng);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle type
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int specie = SPECIES_PROTON_BEAM;
   trajectory->SetSpecie(specie);

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
   double Bmag = 5.0E-6 / unit_magnetic_fluid;
   GeoVector B0(0.0, 0.0, Bmag);
   container.Insert(B0);

// Effective "mesh" resolution
   double enr = 100.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle;
   double R_L = LarmorRadius(Mom(enr, specie), Bmag, specie);
   double dmax = 0.1 * R_L;
   container.Insert(dmax);

   TurbProp turb_prop;

   //! Reference wavenumber corresponding to w_min sampled at the speed of the T/L component (6.65 AU)
   double k0 = 6.3E-14 * unit_length_fluid;

//! Fluctuation variance netween k0 and infinity
   // double var_k0_inf = 3.45E-14 / Sqr(unit_magnetic_fluid);
   double var_k0_inf = 0.2 * Sqr(Bmag);

// Longest wave
   double lambda_max = 10.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;

// Shortest wave
   double lambda_min = 0.002 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;

   turb_prop.kmax = M_2PI / lambda_min;
   turb_prop.kmin = M_2PI / lambda_max;
   turb_prop.slope = 5.0 / 3.0;

// Variance between 2*pi/Lmax and infinity
   double variance = var_k0_inf * pow(k0 / turb_prop.kmin, turb_prop.slope - 1.0);

// A modes
   int nA_modes = 0;
   turb_prop.n_waves = nA_modes;
   turb_prop.variance = 1.0 * variance;
   container.Insert(turb_prop);

// T modes
   int nT_modes = 50;
   turb_prop.n_waves = nT_modes;
   turb_prop.variance = 1.0 * variance;
   container.Insert(turb_prop);

// L modes
   int nL_modes = 0;
   turb_prop.n_waves = nL_modes;
   turb_prop.variance = 1.0 * variance;
   container.Insert(turb_prop);

// I modes
   int nI_modes = 0;
   turb_prop.n_waves = nI_modes;
   turb_prop.variance = 1.0 * variance;
   container.Insert(turb_prop);

   trajectory->AddBackground(BackgroundWaves(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   double init_t = 0.0;
   container.Insert(init_t);

   trajectory->AddInitial(InitialTimeFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Size of the initial cube
   double cube_size = 100.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(-cube_size / 2.0 * gv_ones);
   container.Insert( cube_size / 2.0 * gv_ones);

   trajectory->AddInitial(InitialSpaceBox(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double momentum = Mom(enr, specie);
   GeoVector init_mom(0.0, 0.0, momentum);

   container.Insert(init_mom);

   trajectory->AddInitial(InitialMomentumFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition (end)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings = 1;
   container.Insert(max_crossings);

// Action
   std::vector<int> actions; // empty vector because there are no distributions
   container.Insert(actions);
   
// Duration of the trajectory
   double maxtime = 1000.0 / unit_time_fluid;
   container.Insert(maxtime);

   trajectory->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   trajectory->SetStart();
   trajectory->Integrate();
   trajectory->InterpretStatus();
   
   std::string trajectory_file = "main_test_turb_waves_" + trajectory->GetName() + ".lines";
   std::cout << std::endl;
   std::cout << "TURBULENCE VIA SUPERPOSITION OF WAVES" << std::endl;
   std::cout << "++++++++++++++++++++" << std::endl;
   std::cout << "Trajectory type: " << trajectory->GetName() << std::endl;
   std::cout << "Number of modes:" << std::endl;
   std::cout << "\t Alfven/Slab  = " << nA_modes << std::endl;
   std::cout << "\t Transverse   = " << nT_modes << std::endl;
   std::cout << "\t Longitudinal = " << nL_modes << std::endl;
   std::cout << "\t Isotropic    = " << nI_modes << std::endl;
   std::cout << "++++++++++++++++++++" << std::endl;
   std::cout << "Trajectory outputed to " << trajectory_file << std::endl;
   std::cout << std::endl;

   trajectory->PrintCSV(trajectory_file,false);
   
   return 0;
};
