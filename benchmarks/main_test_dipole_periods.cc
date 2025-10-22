
#include "common/fields.hh"
#include "src/background_dipole.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/boundary_momentum.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
// todo when all trajectories are updated
//#include "src/trajectory.hh"
#include "src/trajectory_guiding_.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{

   DataContainer container;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Set the types
//----------------------------------------------------------------------------------------------------------------------------------------------------

//   using HConfig = HConfigDefault;

   using HConfig = HConfig<
         BuildMode::debug,
         SpecieId::proton_core,
         1,1
         ,
         MyTrajectory,
         BackgroundDefault,
         DiffusionDefault,
         DistributionDefault
      >;

//   using Fields = Fields<Elc_t, Mag_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotHatMag_t, DotAbsMag_t>;
   using Trajectory = TrajectoryGuiding<HConfig>;
   using Background = BackgroundDipole<HConfig>;

   using InitialTime = InitialTimeFixed<HConfig>;
   using InitialSpace = InitialSpaceFixed<HConfig>;
   using InitialMomentum = InitialMomentumRing<Trajectory>;
//   using Diffusion = DiffusionIsotropicConstant<Trajectory>;

   using BoundaryTime = BoundaryTimeExpire<Trajectory>;
   using Boundary1 = BoundarySphereAbsorb<Trajectory>;
   using Boundary2 = BoundaryPlanePass<Trajectory>;
   using Boundary3 = BoundaryMirror<Trajectory>;


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<Trajectory> trajectory = std::make_unique<Trajectory>();

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
   double Bmag = 0.311 / unit_magnetic_fluid;
   GeoVector B0(0.0, 0.0, Bmag);
   container.Insert(B0);

// Effective "mesh" resolution
   double RE = 6.37e8 / unit_length_fluid;
   double dmax_fraction = 0.1;
   double dmax = dmax_fraction * RE;
   container.Insert(dmax);

// Reference equatorial distance
   container.Insert(RE);

// dmax fraction for distances closer to the dipole
   container.Insert(dmax_fraction);

   trajectory->AddBackground(Background(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   double init_t = 0.0;
   container.Insert(init_t);

   trajectory->AddInitial(InitialTime(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   double L = 3.0;
   GeoVector start_pos(L*RE, 0.0, 0.0);
   container.Insert(start_pos);

   trajectory->AddInitial(InitialSpace(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double MeV_kinetic_energy = 1.0;
   container.Insert(Mom<specie>(MeV_kinetic_energy * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle));

   double theta_eq = DegToRad(30.0);
   container.Insert(theta_eq);

   trajectory->AddInitial(InitialMomentum(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition (end)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions; // empty vector because there are no distributions
   container.Insert(actions);
   
// Duration of the trajectory
   double drift_period = 3600.0 * 1.05 / MeV_kinetic_energy / L / (1.0 + 0.43 * sin(theta_eq)) / unit_time_fluid;
   double bounce_period = 2.41 * L * (1.0 - 0.43 * sin(theta_eq)) / sqrt(MeV_kinetic_energy) / unit_time_fluid;
   double maxtime = 10.0 * drift_period;
   container.Insert(maxtime);

   trajectory->AddBoundary(BoundaryTime(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Space boundary condition 1 (Earth)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_Earth = 1;
   container.Insert(max_crossings_Earth);

// Action
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Radius
   container.Insert(RE);

   trajectory->AddBoundary(Boundary1(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Space boundary condition 2 (drift)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings = -1;
   container.Insert(max_crossings);

// Action
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Normal
   GeoVector normal_drift(1.0,0.0,0.0);
   container.Insert(normal_drift);

   trajectory->AddBoundary(Boundary2(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Space boundary condition 3 (bounce)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   container.Insert(max_crossings);

// Action
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Normal
   GeoVector normal_bounce(0.0,0.0,1.0);
   container.Insert(normal_bounce);

   trajectory->AddBoundary(Boundary2(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum boundary condition (bounce)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   container.Insert(max_crossings);

// Action
   container.Insert(actions);

   trajectory->AddBoundary(Boundary3(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   trajectory->SetStart();
   trajectory->Integrate();
   trajectory->InterpretStatus();

   std::string trajectory_file = "fields_main_test_dipole_drifts_" + trajectory->GetName() + ".lines";
   std::cout << std::endl;
   std::cout << "DIPOLE FIELD DRIFT PERIODS" << std::endl;
   std::cout << "=========================================================" << std::endl;
   std::cout << "Trajectory type: " << trajectory->GetName() << std::endl;
   std::cout << "Time elapsed (simulated)     = " << trajectory->ElapsedTime() * unit_time_fluid << " s" << std::endl;
   std::cout << "drift period (theory)        = " << drift_period * unit_time_fluid << " s" << std::endl;
   std::cout << "drift period (simulation)    = " << 2.0 * trajectory->ElapsedTime() * unit_time_fluid / trajectory->Crossings(1,1) << " s" << std::endl;
   std::cout << "bounce period (theory)       = " << bounce_period * unit_time_fluid << " s" << std::endl;
   std::cout << "bounce period (simulation 1) = " << 2.0 * trajectory->ElapsedTime() * unit_time_fluid / trajectory->Crossings(1,2) << " s" << std::endl;
   std::cout << "bounce period (simulation 2) = " << 2.0 * trajectory->ElapsedTime() * unit_time_fluid / trajectory->Mirrorings() << " s" << std::endl;
   std::cout << "=========================================================" << std::endl;
   std::cout << "Trajectory outputed to " << trajectory_file << std::endl;
   std::cout << std::endl;

   trajectory->PrintCSV(trajectory_file,false);
   
   return 0;
};
