//
// Created by Lucius Schoenbaum on 8/19/25.
//

#ifndef SPECTRUM_TRAJECTORY_PARAMS_HH
#define SPECTRUM_TRAJECTORY_PARAMS_HH

#include "common/extremadata.hh"

/*!
\author Lucius Schoenbaum
\date 08/16/2025
The flow direction (temporal) for particle trace simulations.
 */
enum class TimeFlow {
   forward,
   backward
};


/*!
\author Lucius Schoenbaum
\date 08/19/2025
The safety level during the simulation, that sets a basket of simulation parameters.
 It can be set via integers 0, 1, 2 or via labels none, weak, strong.
 */
enum TrajectorySafetyLevel {
   none = 0,
   low = 1,
   high = 2
};

template <TrajectorySafetyLevel level>
struct TrajectorySafety;
template<>
struct TrajectorySafety<none> {};
template<>
struct TrajectorySafety<low> {
//   ! Largest length for single trajectory
//   static constexpr int max_trajectory_steps = 10000000;
   // todo - step carefully through source for this case
};
template<>
struct TrajectorySafety<high> {
//! Largest length for single trajectory
   static constexpr int max_trajectory_steps = 10000000;
//! Largest number of time step adaptations for a single time step
   static constexpr int max_time_step_adaptations = 100;
//! Time step adaptations for a single time step
   int time_step_adaptations;
};

template <bool record_trajectory>
struct TrajectoryRecords;
template <>
struct TrajectoryRecords<true> {
//! Number of trajectory segments (transient)
   int n_segs;
};
template <>
struct TrajectoryRecords<false> {};



// todo if the entire class 'params' is declared constexpr, then the
//  fields in TrajectoryParams do not need to be declared constexpr. <-- check

struct TrajectoryParams {

   TimeFlow timeflow = TimeFlow::forward;

   bool record_trajectory = false;

   bool record_bmag_extrema = true;

   TrajectorySafetyLevel traj_adv_safety_level = TrajectorySafetyLevel::low;

};


#endif
