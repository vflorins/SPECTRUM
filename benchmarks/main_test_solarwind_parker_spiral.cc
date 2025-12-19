#include "src/background_solarwind.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/boundary_momentum.hh"
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

   Specie<default_specie> specie;

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
   double umag = 4.0e7 / Particle::unit_velocity;
   GeoVector u0(umag, 0.0, 0.0);
   container.Insert(u0);

// Magnetic field
   double RS = 6.957e10 / Particle::unit_length;
   double r_ref = 3.0 * RS;
   double BmagE = 5.0e-5 / Particle::unit_magnetic;
   double one_au = SPC_CONST_CGSM_ASTRONOMICAL_UNIT / Particle::unit_length;
   double Bmag_ref = BmagE * Sqr(one_au / r_ref);
   GeoVector B0(Bmag_ref, 0.0, 0.0);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax = one_au;
   container.Insert(dmax);

// Solar rotation vector
   double w0 = M_2PI / (25.0 * 24.0 * 3600.0) * Particle::unit_time;
   GeoVector Omega(0.0, 0.0, w0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction for distances closer to the dipole
   double dmax_fraction = 0.01;
   container.Insert(dmax_fraction);

   trajectory->AddBackground(BackgroundSolarWind(), container);

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

   GeoVector start_pos(r_ref,0.0,0.0);
   container.Insert(start_pos);

   trajectory->AddInitial(InitialSpaceFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double keV_kinetic_energy = 1.0;
   container.Insert(Particle::Mom<specie>(keV_kinetic_energy * SPC_CONST_CGSM_KILO_ELECTRON_VOLT / Particle::unit_energy));

   trajectory->AddInitial(InitialMomentumShell(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Space boundary condition 1 (source surface)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings = 1;
   container.Insert(max_crossings);

// Action
   std::vector<int> actions; // empty vector because there are no distributions
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Radius
   container.Insert(RS);

   trajectory->AddBoundary(BoundarySphereAbsorb(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Space boundary condition 2 (termination shock)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   container.Insert(max_crossings);

// Action
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Radius
   double r_out = 10.0 * one_au;
   container.Insert(r_out);

   trajectory->AddBoundary(BoundarySphereAbsorb(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   trajectory->SetStart();
   trajectory->Integrate();
   trajectory->InterpretStatus();

   std::string trajectory_file = "main_test_solarwind_parker_spiral_" + trajectory->GetName() + ".lines";
   std::cout << std::endl;
   std::cout << "PARKER SPIRAL SOLAR WIND" << std::endl;
   std::cout << "++++++++++++++++++++" << std::endl;
   std::cout << "Trajectory type: " << trajectory->GetName() << std::endl;
   std::cout << "Maximum radial distance  = " << r_out << " AU" << std::endl;
   std::cout << "Time elapsed (simulated) = " << trajectory->ElapsedTime() * Particle::unit_time << " s" << std::endl;
   std::cout << "++++++++++++++++++++" << std::endl;
   std::cout << "Trajectory outputed to " << trajectory_file << std::endl;
   std::cout << std::endl;

   trajectory->PrintCSV(trajectory_file, true);
   
   return 0;
};
