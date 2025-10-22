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


template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Lorentz, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::cartesian>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Lorentz, specieid> {
   using type = Fields<FConfig<specieid>, Mag_t, Elc_t>;
};


template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Focused, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Focused, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelFluv_t, DelAbsMag_t, DelMag_t, DotFluv_t>;
};



template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Guiding, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Guiding, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>;
};



template <SpecieId specieid>
class TrajectoryCoordinates<TrajectoryId::Parker, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>;
};
template <SpecieId specieid>
class TrajectoryFields<TrajectoryId::Parker, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t>;
};


/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory class
\author Lucius Schoenbaum
\date 09/05/2025
*/
template<
   SpecieId specieid,
   TrajectoryId trajectoryid,
   typename Fields_,
   typename RecordCoordinates_,
   TrajectoryOptions::TimeFlow time_flow_,
   RKIntegrator rk_integrator_,
   bool record_mag_extrema_,
   bool record_trajectory_,
   size_t record_trajectory_segment_presize_,
   TrajectoryOptions::SafetyLevel trajectory_adv_safety_level_,
   int max_trajectory_steps_,
   int max_time_adaptations_,
   int n_max_calls_,
   Ratio cfl_advection_,
   Ratio cfl_diffusion_,
   Ratio cfl_acceleration_,
   Ratio cfl_pitchangle_,
   Ratio drift_safety_,
   int mirror_threshold_,
   TrajectoryOptions::PPerpMethod pperp_method_,
   TrajectoryOptions::UseBDrifts use_B_drifts_,
   TrajectoryOptions::StochasticMethod stochastic_method_,
   TrajectoryOptions::StochasticMethod stochastic_method_mu_,
   TrajectoryOptions::StochasticMethod stochastic_method_perp_,
   Ratio split_scatt_fraction_,
   TrajectoryOptions::ConstDmumax const_dmumax_,
   int steps_per_orbit_,
   TrajectoryOptions::DivkMethod divk_method_,
   Ratio dlnp_max_
>
struct TrajectoryConfig {

   using Coordinates = TrajectoryCoordinates<trajectoryid, specieid>;
   using Fields = Cond<std::same_as<Fields_, Default>, TrajectoryFields<trajectoryid, specieid>, Fields_>;
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
   static constexpr auto split_scatt_fraction = get_fp<split_scatt_fraction_>();
   static constexpr auto const_dmumax = const_dmumax_;
   static constexpr auto steps_per_orbit = steps_per_orbit_;
   static constexpr auto divk_method = divk_method_;
   static constexpr auto dlnp_max = get_fp<dlnp_max_>();
};




}


#endif