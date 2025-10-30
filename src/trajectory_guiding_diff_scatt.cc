/*!
\file trajectory_guiding_diff_scatt.cc
\brief Declares a class for trajectory based on guiding center equations with perpendicular diffusion and pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_guiding_diff_scatt.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingDiffScatt methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
TrajectoryGuidingDiffScatt::TrajectoryGuidingDiffScatt(void)
                          : TrajectoryGuidingDiff(traj_name_guidingdiffscatt, 0, STATE_NONE, defsize_guidingdiffscatt),
                            TrajectoryGuidingScatt(traj_name_guidingdiffscatt, 0, STATE_NONE, defsize_guidingdiffscatt),
                            TrajectoryGuiding(traj_name_guidingdiffscatt, 0, STATE_NONE, defsize_guidingdiffscatt)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
*/
bool TrajectoryGuidingDiffScatt::IsSimmulationReady(void) const
{
   return TrajectoryGuidingDiff::IsSimmulationReady();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
void TrajectoryGuidingDiffScatt::DiffusionCoeff(void)
try {
   TrajectoryGuidingDiff::DiffusionCoeff();
   TrajectoryGuidingScatt::DiffusionCoeff();
}

catch(ExFieldError& exception) {
//   PrintError(__FILE__, __LINE__, "Error in field increment evaluation", true);
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
void TrajectoryGuidingDiffScatt::PhysicalStep(void)
{
   TrajectoryGuidingDiff::PhysicalStep();
   TrajectoryGuidingScatt::PhysicalStep();
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 04/29/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryGuidingDiffScatt::Advance(void)
{
// Retrieve latest point of the trajectory
   Load();

// Compute drift and diffusion coefficients for physical step computation
   DiffusionCoeff();
   DriftCoeff();

   PhysicalStep();
   dt = fmin(dt_physical, dt_adaptive);
   TimeBoundaryBefore();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// First half of stochastic pitch angle contribution and advection term
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Perform first half of PA scattering
#if STOCHASTIC_METHOD_MU == 0
   EulerPitchAngleScatt(0);
#elif STOCHASTIC_METHOD_MU == 1
   MilsteinPitchAngleScatt(0);
#elif STOCHASTIC_METHOD_MU == 2
   RK2PitchAngleScatt(0);
#endif

// Store position and momentum locally
   StoreLocal();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute the RK slopes.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The commomn fields and "dmax" have been computed at the end of Advance() or in SetStart() before the first step.
// Compute the slopes. The first two components for momentum are always zero for GC (the perpendicular momentum is determined from conservation of magnetic moment).
   Slopes(slope_pos[0], slope_mom[0], slope_amp[0], slope_wgt[0]);

// Stochastic RK slopes
#if STOCHASTIC_METHOD_PERP == 0
   EulerPerpDiffSlopes();
#elif STOCHASTIC_METHOD_PERP == 1
   MilsteinPerpDiffSlopes();
#elif STOCHASTIC_METHOD_PERP == 2
   if (RK2PerpDiffSlopes()) return true;
#endif

// If trajectory terminated (or is invalid) while computing slopes, exit advance function with true (step was taken)
   if (RKSlopes()) return true;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Advance the trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If adaptive method error is unacceptable, exit advance function with false (step was not taken)
   if (RKStep()) return false;

// Stochastic displacement
   _pos += dr_perp;

#ifdef SPLIT_SCATT

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Second half of stochastic pitch angle contribution and advection term
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If an exit spatial boundary was crossed, the fields may no longer be available, so the second half of the scattering cannot be performed. In that case the function should return immediately and the current position will be saved.
   if (SpaceTerminateCheck()) return true;

// Compute diffusion coefficients
   CommonFields();
   TrajectoryGuidingScatt::DiffusionCoeff();

// Perform second half of PA scattering
#if STOCHASTIC_METHOD_MU == 0
   EulerPitchAngleScatt(1);
#elif STOCHASTIC_METHOD_MU == 1
   MilsteinPitchAngleScatt(1);
#elif STOCHASTIC_METHOD_MU == 2
   RK2PitchAngleScatt(1);
#endif

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------
   HandleBoundaries();

// If trajectory is not finished (in particular, spatial boundary not crossed), the fields can be computed and momentum corrected
   if (BITS_LOWERED(_status, TRAJ_FINISH)) {
      CommonFields();
      MomentumCorrection();
   };

// Add the new point to the trajectory.
   Store();

   return true;
};

};
