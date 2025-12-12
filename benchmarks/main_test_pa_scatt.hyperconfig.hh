// File main_test_pa_scatt.hyperconfig.hh
// Uniform GuidingScatt IsotropicConstant

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
   static constexpr auto background = Config::Background::Uniform;
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
};

struct TrajectoryConfig1{
//! Name: trajectory
// Description: The trajectory type used for simulation
// Options: Fieldline | Focused | Guiding | GuidingDiff | GuidingScatt | GuidingDiffScatt | Lorentz | Parker
   static constexpr auto trajectory = Config::Trajectory::GuidingScatt;
//! Name: Coordinates
// Description: The coordinates of the trajectory during the simulation.
   using Coordinates = Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
//! Name: RecordCoordinates
// Description: The coordinates (and format spec) used to record the trajectory progress.
   using RecordCoordinates = Fields<FConfig<>, Pos_t, Time_t>;
//! Name: Fields
// Description: The fields computed in the local environment of the trajectory during the simulation.
   using Fields = Fields<FConfig<specieid_>, Fluv_t, Mag_t, Ele_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>;
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
//! Name: pperp_method
// Description: Switch controlling how to calculate mu. updating mu according to the scheme does not guarantee conservation of magnetic moment, but can be used with non-adiabatic terms.
// Options: moment_cons: compute mu from magnetic moment conservation | scheme: update according to scheme
   static constexpr auto pperp_method = TrajectoryOptions::PPerpMethod::scheme;
//! Name: cfl_advection
// Description: CFL condition for advection
   static constexpr double cfl_advection = 0.5;
//! Name: drift_safety
// Description: Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
   static constexpr double drift_safety = 0.5;
//! Name: mirror_threshold
// Description: How many time steps to allow before recording a mirror event
   static constexpr int mirror_threshold = 10;
//! Name: split_scatt_fraction
// Description: Whether to split the diffusive advance into two (one before and one after the advection).
// Options: 0.0: do not split | >0.0: fraction of stochastic step to take before deterministic step
   static constexpr double split_scatt_fraction = 0.0;
//! Name: const_dmumax
// Description: Desired accuracy in pitch angle cosine or in pitch angle mu
// Options: constant_dtheta_max: dtheta_max = 2π/180 (deg to rad conversion factor) | constant_dmumax: dmumax = 0.02 (desired accuracy in pitch angle cosine)
   static constexpr auto const_dmumax = TrajectoryOptions::ConstDmumax::constant_dtheta_max;
//! Name: stochastic_method
// Description: Which stochastic method to use for scattering
// Options: Euler | Milstein | RK2
   static constexpr auto stochastic_method = TrajectoryOptions::StochasticMethod::Euler;
//! Name: cfl_pitchangle
// Description: CFL condition for pitch angle scattering
   static constexpr double cfl_pitchangle = 0.5;
};

struct DiffusionConfig1{
//! Name: diffusion
// Description: The diffusion type used for simulation
// Options: None | IsotropicConstant | QLTConstant | WNLTConstant | WNLTRampVLISM | ParaConstant | PerpConstant | FullConstant | FlowMomentumPowerLaw | KineticEnergyRadialDistancePowerLaw | RigidityMagneticFieldPowerLaw | StraussEtAl2013 | GuoEtAl2014 | PotgieterEtAl2015 | EmpiricalSOQLTandUNLT
   static constexpr auto diffusion = Config::Diffusion::IsotropicConstant;
//! Name: Coordinates
// Description: The coordinates where the diffusion is computed during the simulation.
   using Coordinates = Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>;
//! Name: Fields
// Description: The fields computed in the local environment for the diffusion computation during the simulation.
   using Fields = Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>;
//! Name: D0
// Description: diffusion coefficient (if constant)
   static constexpr double D0 = 1.0;
};



using HConfig = HyperConfigure<
      SimulationConfig1,
      BackgroundConfig1,
      TrajectoryConfig1,
      DiffusionConfig1
>;

