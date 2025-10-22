/*!
\file trajectory.config.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory class
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/


#ifndef SPECTRUM_TRAJECTORY_CONFIG_HH
#define SPECTRUM_TRAJECTORY_CONFIG_HH

#include "common/compiletime_lists.hh"
#include "common/fields.hh"

namespace Spectrum {

template <TrajectoryId trajectoryid, SpecieId specieid>
class TrajectoryCoordinates;
template <TrajectoryId trajectoryid, SpecieId specieid>
class TrajectoryFields;


template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Fieldline, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Fieldline, specieid> {
   /* Unused */
   using type = Fields<FConfig<specieid>>;
};
// TODO DONE


template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Lorentz, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::cartesian>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Lorentz, specieid> {
   using type = Fields<FConfig<specieid>, Mag_t, Elc_t>;
};
// TODO DONE


template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Focused, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Focused, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelFluv_t, DelAbsMag_t, DelMag_t, DotFluv_t>;
};
// TODO DONE



template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Guiding, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Guiding, specieid> {
   using type = Fields<FConfig<specieid>, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelAbsMag_t, DotAbsMag_t>;
};



template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Parker, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Parker, specieid> {
   using type = Fields<FConfig<specieid>, Todo_t>;
};


/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory class
\author Lucius Schoenbaum
\date 09/05/2025
*/
template<
      SpecieId specieid,
      TrajectoryId trajectoryid,
      typename TrajectoryFields,
      typename RecordCoordinates_,
      TrajectoryOptions::TimeFlow time_flow_,
      RKIntegrator rk_integrator_,
//! Whether to record magnetic field extrema
    bool record_mag_extrema_,
//! Whether to record segments
    bool record_trajectory_,
//! Default initial size of trajectory segment record
      size_t record_trajectory_segment_presize_,
//! Trajectory advance safety level: 0 means no checks, 1 means check dt only, 2 means check dt, number of segments, and time adaptations per step.
      TrajectoryOptions::SafetyLevel trajectory_adv_safety_level_,
//! Largest length for single trajectory
      int max_trajectory_steps_,
//! Largest number of time step adaptations for a single time step
      int max_time_adaptations_,
//! Upper limit on the number of steps in debug mode
//! use -1 for unlimited
      int n_max_calls_,
//! CFL condition for advection
    Ratio cfl_advection_,
//! CFL condition for diffusion
    Ratio cfl_diffusion_,
//! CFL condition for acceleration
    Ratio cfl_acceleration_,
//! CFL condition for pitch angle scattering
      Ratio cfl_pitchangle_,
//! Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
      Ratio drift_safety_,
//! How many time steps to allow before recording a mirror event
      int mirror_threshold_,
//// Switch controlling how to calculate mu. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
      TrajectoryOptions::PPerpMethod pperp_method_,
//! Flag to use gradient and curvature drifts in drift velocity calculation
      TrajectoryOptions::UseBDrifts use_B_drifts_,
//! Which stochastic method to use for PA scattering, 0 = Euler, 1 = Milstein, 2 = RK2
//! Which stochastic method to use for perpendicular diffusion, 0 = Euler, 1 = Milstein, 2 = RK2
    TrajectoryOptions::StochasticMethod stochastic_method_,
    TrajectoryOptions::StochasticMethod stochastic_method_mu_,
    TrajectoryOptions::StochasticMethod stochastic_method_perp_,
//! Whether to split the diffusive advance into two (one before and one after the advection).
// #define SPLIT_SCATT
      bool split_scatt_,
//#ifdef SPLIT_SCATT
//! Fraction of stochastic step to take before deterministic step
//const double alpha = 0.5;
//#endif
      Ratio split_scatt_fraction_,
//#if CONST_DMUMAX == 1
//! Desired accuracy in pitch angle cosine
//const double dmumax = 0.02;
//#else
//! Desired accuracy in pitch angle (deg x [deg to rad conversion factor])
//const double dthetamax = 2.0 * M_PI / 180.0;
//#endif
      TrajectoryOptions::ConstDmumax const_dmumax_,
//! Number of time steps per one orbit
    int steps_per_orbit_,
//! Which method of computation to use for divK: 0 = using direct central FD, 1 = using _spdata.grad quantities
    TrajectoryOptions::DivkMethod divk_method_,
//! Maximum allowed fraction of momentum change per step
      Ratio dlnp_max_
>
struct TrajectoryConfig {

   using Coordinates = TrajectoryCoordinates<trajectoryid, specieid>;
   using Fields = TrajectoryFields;
   using RecordCoordinates = RecordCoordinates_;

   static constexpr auto time_flow = time_flow_;
   RKIntegrator rk_integrator = rk_integrator_;
   static constexpr auto record_mag_extrema = record_mag_extrema_;
   static constexpr auto record_trajectory = record_trajectory_;
   static constexpr auto record_trajectory_segment_presize = record_trajectory_segment_presize_;
   static constexpr auto trajectory_adv_safety_level = trajectory_adv_safety_level_;
   static constexpr auto max_trajectory_steps = max_trajectory_steps_;
   static constexpr auto max_time_adaptations = max_time_adaptations_;
   static constexpr auto n_max_calls = n_max_calls_;
   static constexpr auto cfl_advection = get_fp<cfl_advection_>();
   static constexpr auto cfl_diffusion = get_fp<cfl_diffusion_>();
   static constexpr auto cfl_acceleration = get_fp<cfl_acceleration_>();
   static constexpr auto cfl_pitchangle = get_fp<cfl_pitchangle_>();
   static constexpr auto drift_safety = get_fp<drift_safety_>();
   static constexpr auto mirror_threshold = get_fp<mirror_threshold_>();
   static constexpr auto pperp_method = pperp_method_;
   static constexpr auto use_B_drifts = use_B_drifts_;
   static constexpr auto stochastic_method = stochastic_method_;
   static constexpr auto stochastic_method_mu = stochastic_method_mu_;
   static constexpr auto stochastic_method_perp = stochastic_method_perp_;
   static constexpr auto split_scatt = split_scatt_;
   static constexpr auto split_scatt_fraction = get_fp<split_scatt_fraction_>();
   static constexpr auto const_dmumax = const_dmumax_;
   static constexpr auto steps_per_orbit = steps_per_orbit_;
   static constexpr auto divk_method = divk_method_;
   static constexpr auto dlnp_max = get_fp<dlnp_max_>();
};




}


#endif