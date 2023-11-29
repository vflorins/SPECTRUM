/*!
\file trajectory_guiding_scatt.cc
\brief Defines a class for trajectory based on guiding center equations with pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_guiding_scatt.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
TrajectoryGuidingScatt::TrajectoryGuidingScatt(void)
                      : TrajectoryGuiding(traj_name_guidingscatt, 0, STATE_NONE, defsize_guidingscatt)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
\param[in] presize_in Whether to pre-allocate memory for trajectory arrays
*/
TrajectoryGuidingScatt::TrajectoryGuidingScatt(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                      : TrajectoryGuiding(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
*/
bool TrajectoryGuidingScatt::IsSimmulationReady(void) const
{
   if(!TrajectoryBase::IsSimmulationReady()) return false;

// A diffusion object is required
   if(diffusion == nullptr) return false;
   else if(BITS_LOWERED(diffusion->GetStatus(), STATE_SETUP_COMPLETE)) return false;
   return true;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
void TrajectoryGuidingScatt::DiffusionCoeff(void)
try {
   Dmumu = diffusion->GetComponent(2, _t, _pos, ConvertMomentum(), _spdata);

// Compute the derivative in mu
   Vmu = diffusion->GetMuDerivative();
}

catch(ExFieldError& exception) {
//   PrintError(__FILE__, __LINE__, "Error in field increment evaluation", true);
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/10/2022
\param[in] second True if this is the second step of a split scheme
*/
void TrajectoryGuidingScatt::EulerPitchAngleScatt(bool second)
{
   double mu_new, dt_local;
   GeoVector mom_conv = ConvertMomentum();

#ifdef SPLIT_SCATT
   dt_local = (second ? 1.0 - alpha : alpha) * dt;
#else
   dt_local = dt;
#endif

// TODO The loop ensures that "mu" does not become larger than 1, but a simple reflection might be sufficient
   do {
      mu_new = mom_conv[1] + Vmu * dt_local + sqrt(2.0 * Dmumu * dt_local) * rng->GetNormal();

#ifdef GEO_DEBUG
      if(fabs(mu_new) > 1.0) Nabsmugt1 += 1;
#endif
   } while(fabs(mu_new) > 1.0);

   _mom[0] = mom_conv[0] * sqrt(1 - Sqr(mu_new));
   _mom[2] = mom_conv[0] * mu_new;
   _vel = Vel(_mom, specie);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/10/2022
\param[in] second True if this is the second step of a split scheme
*/
void TrajectoryGuidingScatt::MilsteinPitchAngleScatt(bool second)
{
   double mu_new, dmu, mu = _mom[2] / _mom.Norm();
   double dt_local, b, b1, dW, Dmumu_new;
   GeoVector mom_conv = ConvertMomentum();

#ifdef SPLIT_SCATT
   dt_local = (second ? 1.0 - alpha : alpha) * dt;
#else
   dt_local = dt;
#endif

// Get diffusion slope and derivative of diffusion slope
   b = sqrt(2.0 * Dmumu);
   dmu = small * (mu + small < 1.0 ? 1.0 : -1.0);
   mom_conv[1] += dmu;
   Dmumu_new = diffusion->GetComponent(2, _t, _pos, mom_conv, _spdata);
   b1 = (sqrt(2.0 * Dmumu_new) - b) / dmu;

// TODO The loop ensures that "mu" does not become larger than 1, but a simple reflection might be sufficient
   do {
      dW = sqrt(dt_local) * rng->GetNormal();
      mu_new = mu + Vmu * dt_local + b * dW + 0.5 * b * b1 * (Sqr(dW) - dt_local);

#ifdef GEO_DEBUG
      if(fabs(mu_new) > 1.0) Nabsmugt1 += 1;
#endif
   } while(fabs(mu_new) > 1.0);

   _mom[0] = mom_conv[0] * sqrt(1 - Sqr(mu_new));
   _mom[2] = mom_conv[0] * mu_new;
   _vel = Vel(_mom, specie);
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2022
\param[in] second True if this is the second step of a split scheme

Computes RK stochastic step with 2 evaluations of the drift term, 3 evaluations of the variance term, and 1 derivative of the variance term
*/
void TrajectoryGuidingScatt::RK2PitchAngleScatt(bool second)
{
   double slope_Vmu[2], slope_Dmumu[3];
   double gambar = 1.0 / sqrt(3.0);
   double dmu, mu_new, mu = _mom[2] / _mom.Norm();
   double Dmumu_new, dfit, dW, dt_local;
   GeoVector mom_conv = ConvertMomentum();

#ifdef SPLIT_SCATT
   dt_local = (second ? 1.0 - alpha : alpha) * dt;
#else
   dt_local = dt;
#endif

// Vmu and Dmumu at current position should already be computed
   slope_Vmu[0] = Vmu;
   slope_Dmumu[0] = sqrt(2.0 * Dmumu);

// TODO The loop ensures that "mu" does not become larger than 1, but a simple reflection might be sufficient
   do {

// Generate stochastic factor
      dW = sqrt(dt_local) * rng->GetNormal();

// Compute second slope for Vmu
      mu_new = mu + slope_Dmumu[0] * dW;
      mu_new += slope_Vmu[0] * dt_local;
      if(fabs(mu_new) > 1.0) {
#ifdef GEO_DEBUG
         Nabsmugt1 += 1;
#endif
         continue;
      };

      mom_conv[1] = mu_new;
      Dmumu = diffusion->GetComponent(2, _t + dt_local, _pos, mom_conv, _spdata);

      dmu = small * (mom_conv[1] + small < 1.0 ? 1.0 : -1.0);
      mom_conv[1] += dmu;

      Dmumu_new = diffusion->GetComponent(2, _t + dt_local, _pos, mom_conv, _spdata);
      slope_Vmu[1] = (Dmumu_new - Dmumu) / dmu;

// Compute second slope for Dmumu
      mu_new = mu + gambar * slope_Dmumu[0] * dW;
      mu_new += slope_Vmu[0] * dt_local;
      if(fabs(mu_new) > 1.0) {
#ifdef GEO_DEBUG
         Nabsmugt1 += 1;
#endif
         continue;
      };

      mom_conv[1] = mu_new;
      Dmumu = diffusion->GetComponent(2, _t + dt_local, _pos, mom_conv, _spdata);
      slope_Dmumu[1] = sqrt(2.0 * Dmumu);

// Compute third slope for Dmumu
      mu_new = mu - slope_Dmumu[0] / (3.0 * gambar) * dW;
      mu_new += slope_Vmu[0] * dt_local;
      if(fabs(mu_new) > 1.0) {
#ifdef GEO_DEBUG
         Nabsmugt1 += 1;
#endif
         continue;
      };

      mom_conv[1] = mu_new;
      Dmumu = diffusion->GetComponent(2, _t + dt_local, _pos, mom_conv, _spdata);
      slope_Dmumu[2] = sqrt(2.0 * Dmumu);

// Compute additional fit term
      dmu = small * (mu + small < 1.0 ? 1.0 : -1.0);
      mom_conv[1] = mu + dmu;
      Dmumu_new = diffusion->GetComponent(2, _t, _pos, mom_conv, _spdata);
      dfit = (sqrt(2.0 * Dmumu_new) - slope_Dmumu[0]) / dmu;

// Add all the stuff
      mu_new = mu + (0.5 * slope_Dmumu[0] + (slope_Dmumu[1] + 3.0 * Sqr(gambar) * slope_Dmumu[2]) / (2.0 + 6.0 * Sqr(gambar)) ) * dW;
      mu_new += 0.5 * dt_local * (slope_Vmu[0] + slope_Vmu[1]);
      mu_new += 0.5 * slope_Dmumu[0] * dfit * (Sqr(dW) - dt_local);

#ifdef GEO_DEBUG
      if(fabs(mu_new) > 1.0) Nabsmugt1 += 1;
#endif
   } while(fabs(mu_new) > 1.0);

   _mom[0] = mom_conv[0] * sqrt(1 - Sqr(mu_new));
   _mom[2] = mom_conv[0] * mu_new;
   _vel = Vel(_mom, specie);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
void TrajectoryGuidingScatt::PhysicalStep(void)
{
#if CONST_DMUMAX == 0
   double dmumax = sqrt(1 - fabs(_mom[2] / _mom.Norm())) * dthetamax + 0.5 * fabs(_mom[2] / _mom.Norm()) * Sqr(dthetamax);
#endif
   dt_physical = fmin(dt_physical, cfl_pa_gs * Sqr(dmumax) / Dmumu);
   dt_physical = fmin(dt_physical, cfl_pa_gs * dmumax / fabs(Vmu));
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 04/29/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryGuidingScatt::Advance(void)
{
// Retrieve latest point of the trajectory
   Load();

// Compute drift and diffusion coefficients (including PA advection term) for physical step computation
   DiffusionCoeff();
   DriftCoeff();

   TrajectoryGuiding::PhysicalStep();
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
   Slopes(slope_pos[0], slope_mom[0]);

// If trajectory terminated (or is invalid) while computing slopes, exit advance function with true (step was taken)
   if(RKSlopes()) return true;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Advance the trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If adaptive method error is unacceptable, exit advance function with false (step was not taken)
   if(RKStep()) return false;

#ifdef SPLIT_SCATT

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Second half of stochastic pitch angle contribution and advection term
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If an exit spatial boundary was crossed, the fields may no longer be available, so the second half of the scattering cannot be performed. In that case the function should return immediately and the current position will be saved.
   if(SpaceTerminateCheck()) return true;

// Compute diffusion coefficients (including PA advection term)
   CommonFields();
   DiffusionCoeff();

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
   if(BITS_LOWERED(_status, TRAJ_FINISH)) {
      CommonFields();
      MomentumCorrection();
   };

// Add the new point to the trajectory.
   Store();

   return true;
};

};
