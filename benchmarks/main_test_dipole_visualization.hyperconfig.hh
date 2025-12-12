// File main_test_dipole_visualization.hyperconfig.hh
// Dipole Fieldline None

#include "common/compiletime_lists.hh"
#include "common/vectors.hh"
#include "common/fields.hh"
#include "src/hyperconfigure.hh"

using namespace Spectrum;

constexpr auto specieid_ = SpecieId::proton_core;
constexpr auto specie = Specie<specieid_>();


struct SimulationConfig1 {
   static constexpr auto specieid = specieid_;
   static constexpr auto build_mode = BuildMode::debug;
// Whether to print the last trajectory
   static constexpr auto print_last_trajectory = true;
// Whether there is a supervisor process. This does not require there to be server processes.
   static constexpr auto supervisor = false;
// These can also be set from the command line (argv) at runtime.
   static constexpr auto num_trajectories = 1;
   static constexpr auto batch_size = 1;
   static constexpr auto max_trajectories_per_worker = 1;
};

struct BackgroundConfig1{
//! Name: background
// Description: The background type used for simulation
// Options: CylindricalObstacle | DataBATL | DataCartesian | Dipole | Discontinuity | MagnetizedCylinder | MagnetizedSphere | Shock | SmoothDiscontinuity | SmoothShock | SolarWind | SphericalObstacle | Uniform | VLISMBochum | Waves
   static constexpr auto background = Config::Background::Dipole;
//! Name: derivative_method
// Description: The method used to evaluate derivatives of spatially located field quantities. Usually, the only information checked is whether or not this value is numeric.
// Options: numeric | nonnumeric | analytic | datadefined
   static constexpr auto derivative_method = BackgroundOptions::DerivativeMethod::analytic;
//! Name: num_numeric_grad_evals
// Description: The number of derivative evaluations applied in the derivative-averaging method
   static constexpr int num_numeric_grad_evals = 1;
//! Name: incr_dmax_ratio
// Description: What fraction of _dmax to use to calculate the field increment
   static constexpr double incr_dmax_ratio = 0.0001;
//! Name: dmax0
// Description: baseline simulation-wide dmax value
   static constexpr double dmax0 = 0.1;
//! Name: dmax_fraction
// Description: Fraction of dmax in a reduced region (background dependent)
   static constexpr double dmax_fraction = 0.1;
//! Name: r0
// Description: value (todo: description)
   static constexpr GeoVector r0 = {0.0, 0.0, 0.0};
//! Name: u0
// Description: u0
   static constexpr GeoVector u0 = {0.0, 0.0, 0.0};
//! Name: B0
// Description: B0
   static constexpr GeoVector B0 = {1.0, 1.0, 1.0};
//! Name: r_ref
// Description: Reference equatorial distance
   static constexpr double r_ref = 1.0;
};

struct TrajectoryConfig1{
//! Name: trajectory
// Description: The trajectory type used for simulation
// Options: Fieldline | Focused | Guiding | GuidingDiff | GuidingScatt | GuidingDiffScatt | Lorentz | Parker
   static constexpr auto trajectory = Config::Trajectory::Fieldline;
//! Name: Coordinates
// Description: The coordinates of the trajectory during the simulation.
   using Coordinates = Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
//! Name: RecordCoordinates
// Description: The coordinates (and format spec) used to record the trajectory progress.
   using RecordCoordinates = Fields<FConfig<>, Pos_t, Time_t>;
//! Name: Fields
// Description: The fields computed in the local environment of the trajectory during the simulation.
   using Fields = Fields<FConfig<specieid_>>;
//! Name: FieldlineField_t
// Description: The tracked field for a Fieldline trajectory class.
   using FieldlineField_t = Mag_t;
//! Name: time_flow
// Description: The time flow direction used for simulation
// Options: forward | backward
   static constexpr auto time_flow = TrajectoryOptions::TimeFlow::forward;
//! Name: rk_integrator
// Description: The discretization scheme used for simulation, a Runge-Kutta scheme
// Options: Euler1E | Euler1I | Midpoint_2I | Midpoint_2E | RungeKutta_4E | Kutta38_4E | RungeKuttaFehlberg_54E | CashKarp_54E | DormandPrince_54E | RungeKuttaFehlberg65E | RungeKuttaFehlberg76E | RungeKuttaFehlberg_87E | ...others (see source)
   static constexpr auto rk_integrator = RKIntegrator::DormandPrince_54E;
//! Name: record_mag_extrema
// Description: Whether to record magnetic field extrema
   static constexpr bool record_mag_extrema = false;
//! Name: record_trajectory
// Description: Whether to record segments
   static constexpr bool record_trajectory = false;
//! Name: record_trajectory_segment_presize
// Description: Default initial size of trajectory segment record
   static constexpr int record_trajectory_segment_presize = 10000;
//! Name: advance_safety_level
// Description: Trajectory advance routine safety level
// Options: low: no checks | medium: check dt only | high: check dt, number of segments, and time adaptations per step
   static constexpr auto advance_safety_level = TrajectoryOptions::SafetyLevel::low;
//! Name: max_trajectory_steps
// Description: Largest length for single trajectory
   static constexpr int max_trajectory_steps = 100000;
//! Name: max_time_adaptations
// Description: Largest number of time step adaptations for a single time step
   static constexpr int max_time_adaptations = 1;
//! Name: n_max_calls
// Description: Upper limit on the number of steps in debug mode
// Options: -1: unlimited | n ≥ 0: limited, upper bound n
   static constexpr int n_max_calls = -1;
};

struct DiffusionConfig1{
//! Name: diffusion
// Description: The diffusion type used for simulation
// Options: None | IsotropicConstant | QLTConstant | WNLTConstant | WNLTRampVLISM | ParaConstant | PerpConstant | FullConstant | FlowMomentumPowerLaw | KineticEnergyRadialDistancePowerLaw | RigidityMagneticFieldPowerLaw | StraussEtAl2013 | GuoEtAl2014 | PotgieterEtAl2015 | EmpiricalSOQLTandUNLT
   static constexpr auto diffusion = Config::Diffusion::None;
//! Name: Coordinates
// Description: The coordinates where the diffusion is computed during the simulation.
   using Coordinates = Fields<FConfig<>>;
//! Name: Fields
// Description: The fields computed in the local environment for the diffusion computation during the simulation.
   using Fields = Fields<FConfig<>>;
};



using HConfig = HyperConfigure<
      SimulationConfig1,
      BackgroundConfig1,
      TrajectoryConfig1,
      DiffusionConfig1
>;

