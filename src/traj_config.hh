/*!
\file traj_config.hh
\brief Trajectory configurator
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJ_CONFIG_HH
#define SPECTRUM_TRAJ_CONFIG_HH

#include "config.h"

#if TRAJ_TYPE == TRAJ_FIELDLINE
#include "trajectory_fieldline.hh"
#elif TRAJ_TYPE == TRAJ_LORENTZ
#include "trajectory_lorentz.hh"
#elif TRAJ_TYPE == TRAJ_GUIDING
#include "trajectory_guiding.hh"
#elif TRAJ_TYPE == TRAJ_GUIDING_SCATT
#include "trajectory_guiding_scatt.hh"
#elif TRAJ_TYPE == TRAJ_GUIDING_DIFF
#include "trajectory_guiding_diff.hh"
#elif TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT
#include "trajectory_guiding_diff_scatt.hh"
#elif TRAJ_TYPE == TRAJ_FOCUSED
#include "trajectory_focused.hh"
#elif TRAJ_TYPE == TRAJ_PARKER
#include "trajectory_parker.hh"
#elif TRAJ_TYPE == TRAJ_PARKER_SOURCE
#include "trajectory_parker_source.hh"
#else
#error Unsupported Trajectory type
#endif

#endif
