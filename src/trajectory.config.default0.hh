/*!
\file trajectory.config.default.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory class
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/


#ifndef SPECTRUM_TRAJECTORY_CONFIG_DEFAULT_HH
#define SPECTRUM_TRAJECTORY_CONFIG_DEFAULT_HH

#include "trajectory.config.hh"

namespace Spectrum {




template <TrajectoryId trajectoryid, SpecieId specieid>
using TrajectoryDefault = TrajectoryConfig<
      Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::cartesian>, Pos_t, Time_t, Mom_t>,
      specieid,
      trajectoryid,
      get_timeflow<trajectoryid>(),
      RKIntegrator::DormandPrince_54E,
      /* record_mag_extrema */
      false,
      /* record_trajectory */
      false,
      /* record_trajectory_segment_presize */
      10000,
      /* trajectory_adv_safety_level */
      TrajectoryOptions::SafetyLevel::low,
      /* max_trajectory_steps */
      100000,
      /* max_time_adaptations */
      1,
      /* n_max_calls */
      1000,
      /* cfl_advection */
      std::ratio<1,10>,
      /* cfl_diffusion */
      std::ratio<1,10>,
      /* cfl_acceleration */
      std::ratio<1,10>,
      /* cfl_pitchangle */
      std::ratio<1,10>,
      /* drift_safety */
      std::ratio<1,10>,
      /* mirror_threshold */
      1,
      /* pperp_method */
      TrajectoryOptions::PPerpMethod::mag_moment_conservation,
      /* use_B_drifts */
      TrajectoryOptions::UseBDrifts::none,
      /* stochastic_method */
      TrajectoryOptions::StochasticMethod::Euler,
      /* stochastic_method_mu */
      TrajectoryOptions::StochasticMethod::Euler,
      /* stochastic_method_perp */
      TrajectoryOptions::StochasticMethod::Euler,
      /* split_scatt */
      false,
      /* split_scatt_fraction */
      std::ratio<1,10>,
      /* const_dmumax */
      TrajectoryOptions::ConstDmumax::constant_dmu_max,
      /* steps_per_orbit */
      1,
      /* divk_method */
      TrajectoryOptions::DivkMethod::direct,
      /* dlnp_max */
      std::ratio<1,10>
>;



template <SpecieId specieid>
using TrajectoryGuidingDefault = TrajectoryConfig<
      Fields<FConfig<>>,
      specieid,
      TrajectoryId::Guiding,
      get_timeflow<TrajectoryId::Guiding>(),
      RKIntegrator::DormandPrince_54E,
      /* record_mag_extrema */
      false,
      /* record_trajectory */
      false,
      /* record_trajectory_segment_presize */
      10000,
      /* trajectory_adv_safety_level */
      TrajectoryOptions::SafetyLevel::low,
      /* max_trajectory_steps */
      100000,
      /* max_time_adaptations */
      1,
      /* n_max_calls */
      1000,
      /* cfl_advection */
      std::ratio<1,10>,
      /* cfl_diffusion */
      std::ratio<1,10>,
      /* cfl_acceleration */
      std::ratio<1,10>,
      /* cfl_pitchangle */
      std::ratio<1,10>,
      /* drift_safety */
      std::ratio<1,10>,
      /* mirror_threshold */
      1,
      /* pperp_method */
      TrajectoryOptions::PPerpMethod::mag_moment_conservation,
      /* use_B_drifts */
      TrajectoryOptions::UseBDrifts::none,
      /* stochastic_method */
      TrajectoryOptions::StochasticMethod::Euler,
      /* stochastic_method_mu */
      TrajectoryOptions::StochasticMethod::Euler,
      /* stochastic_method_perp */
      TrajectoryOptions::StochasticMethod::Euler,
      /* split_scatt */
      false,
      /* split_scatt_fraction */
      std::ratio<1,10>,
      /* const_dmumax */
      TrajectoryOptions::ConstDmumax::constant_dmu_max,
      /* steps_per_orbit */
      1,
      /* divk_method */
      TrajectoryOptions::DivkMethod::direct,
      /* dlnp_max */
      std::ratio<1,10>
>;





}

#endif
