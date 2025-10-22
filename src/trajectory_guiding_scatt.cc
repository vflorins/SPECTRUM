/*!
\file trajectory_guiding_scatt.cc
\brief Defines a class for trajectory based on guiding center equations with pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso Guzman

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
template <typename Background, typename Diffusion>
TrajectoryGuidingScatt<Background, Diffusion>::TrajectoryGuidingScatt(void)
      : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Background, typename Diffusion>
TrajectoryGuidingScatt<Background, Diffusion>::TrajectoryGuidingScatt(const std::string& name_in, uint16_t status_in)
      : TrajectoryBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
*/
template <typename Background, typename Diffusion>
bool TrajectoryGuidingScatt<Background, Diffusion>::IsSimulationReady(void) const
{
   if (!TrajectoryBase::IsSimulationReady()) return false;

// A diffusion object is required
   if (diffusion == nullptr) return false;
   else if (BITS_LOWERED(diffusion->GetStatus(), STATE_SETUP_COMPLETE)) return false;
   return true;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
template <typename Background, typename Diffusion>
void TrajectoryGuidingScatt<Background, Diffusion>::DiffusionCoeff(void)
try {

   auto dcoords = DiffusionCoordinates::Convert(_coords);
   DiffusionFields dfields = DiffusionCoordinates::Get(_fields);
   CommonFields_Diffusion(dcoords, dfields);
   diffusion->Stage(dcoords, dfields);
   Dmumu = diffusion->Get(Component::mu);

// Compute the derivative in mu
   Vmu = diffusion->GetMuDerivative(Component::mu);
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
template <typename Background, typename Diffusion>
void TrajectoryGuidingScatt<Background, Diffusion>::EulerPitchAngleScatt(bool second)
{
   double mu_new, dt_local;
   /*
    *
    * Reminder:
    * received Mom() is in coords: perp, para
    * converted Mom() is in coords: p, mu
    *
    */
   auto dcoords = DiffusionCoordinates::Convert(_coords);

   if constexpr (split_scatt) {
      auto alpha = split_scatt_fraction;
      dt_local = (second ? 1.0 - alpha : alpha) * dt;
   }
   else {
      dt_local = dt;
   }

// TODO The loop ensures that "mu" does not become larger than 1, but a simple reflection might be sufficient
   do {
// Vmu is always added because the contribution from this term is outside of the DriftCoeff() function and the "drift_vel" GeoVector
      mu_new = dcoords.MomMu() + Vmu * dt_local + sqrt(2.0 * Dmumu * dt_local) * rng->GetNormal();

      if constexpr (HConfig::build_mode == BuildMode::debug) {
         if (fabs(mu_new) > 1.0) Nabsmugt1 += 1;
      }
   } while (fabs(mu_new) > 1.0);

   _coords = TrajectoryCoordinates::Convert(dcoords);
   // todo review: it performs these
//   _coords.MomPerp() = mom_conv[0] * sqrt(1 - Sqr(mu_new));
//   _coords.MomPara() = mom_conv[0] * mu_new;
//   _coords.Vel() = Vel<specie>(_coords.Mom());
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/10/2022
\param[in] second True if this is the second step of a split scheme
*/
template <typename Background, typename Diffusion>
void TrajectoryGuidingScatt<Background, Diffusion>::MilsteinPitchAngleScatt(bool second)
{
// Calculate Dmumu_new
   auto dcoords = DiffusionCoordinates::Convert(_coords);

   double dmu = sp_small * (dcoords.MomMu() + sp_small < 1.0 ? 1.0 : -1.0);
   dcoords.MomMu() += dmu;
   DiffusionFields dfields = DiffusionCoordinates::Get(_fields);
   CommonFields_Diffusion(dcoords, dfields);
   diffusion->Stage(dcoords, dfields);
   // todo: fields/coords ---- this code is invalid as-is. mom_conv is created and *then* modified, _coords is left as-is - review and fix
   double Dmumu_new = diffusion->Get(Component::mu);

   double mu_new;
   double dt_local, b, b1, dW;

   if constexpr (split_scatt) {
      auto alpha = split_scatt_fraction;
      dt_local = (second ? 1.0 - alpha : alpha) * dt;
   }
   else {
      dt_local = dt;
   }

// Get diffusion slope and derivative of diffusion slope
   b = sqrt(2.0 * Dmumu);
   b1 = (sqrt(2.0 * Dmumu_new) - b) / dmu;

// TODO The loop ensures that "mu" does not become larger than 1, but a simple reflection might be sufficient
   do {
      dW = sqrt(dt_local) * rng->GetNormal();
// Vmu is always added because the contribution from this term is outside of the DriftCoeff() function and the "drift_vel" GeoVector
      mu_new = dcoords.MomMu() + Vmu * dt_local + b * dW + 0.5 * b * b1 * (Sqr(dW) - dt_local);

      if constexpr (HConfig::build_mode == BuildMode::debug) {
         if (fabs(mu_new) > 1.0) Nabsmugt1 += 1;
      }
   } while (fabs(mu_new) > 1.0);

   _coords = TrajectoryCoordinates::Convert(dcoords);
   // todo review: it performs these
//   _coords.MomPerp() = mom_conv[0] * sqrt(1 - Sqr(mu_new));
//   _coords.MomPara() = mom_conv[0] * mu_new;
//   _coords.Vel() = Vel<specie>(_coords.Mom());
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2022
\param[in] second True if this is the second step of a split scheme

Computes RK stochastic step with 2 evaluations of the drift term, 3 evaluations of the variance term, and 1 derivative of the variance term
*/
template <typename Background, typename Diffusion>
void TrajectoryGuidingScatt<Background, Diffusion>::RK2PitchAngleScatt(bool second)
{
// todo review after spdata/fields
   double slope_Vmu[2], slope_Dmumu[3];
   double gambar = 1.0 / sqrt(3.0);
   double dmu, mu_new;
   auto dcoords = DiffusionCoordinates::Convert(_coords);
   auto dcoords_orig = dcoords;
   auto dfields = DiffusionCoordinates();

//   double mu = _coords.MomPara() / _coords.AbsMom();
//   GeoVector mom_conv = ConvertMomentum();
   double Dmumu_new, dfit, dW, dt_local;

   if constexpr (split_scatt) {
      auto alpha = split_scatt_fraction;
      dt_local = (second ? 1.0 - alpha : alpha) * dt;
   }
   else {
      dt_local = dt;
   }

// Vmu and Dmumu at current position should already be computed
   slope_Vmu[0] = Vmu;
   slope_Dmumu[0] = sqrt(2.0 * Dmumu);

// TODO The loop ensures that "mu" does not become larger than 1, but a simple reflection might be sufficient
   do {

// Generate stochastic factor
      dW = sqrt(dt_local) * rng->GetNormal();

// Compute second slope for Vmu
      mu_new = dcoords.MomMu() + slope_Dmumu[0] * dW;
      mu_new += slope_Vmu[0] * dt_local;
      if (fabs(mu_new) > 1.0) {
         if constexpr (HConfig::build_mode == BuildMode::debug) {
            Nabsmugt1 += 1;
         }
         continue;
      };

      dcoords.MomMu() = mu_new;
      CommonFields(dcoords, dfields);
      diffusion->Stage(dcoords, dfields);
      Dmumu = diffusion->Get(Component::mu);

      dmu = sp_small * (dcoords.MomMu() + sp_small < 1.0 ? 1.0 : -1.0);
      dcoords.MomMu() += dmu;
      dcoords.Time() += dt_local; // todo review _______

      CommonFields(dcoords, dfields);
      diffusion->Stage(dcoords, dfields);
      Dmumu_new = diffusion->Get(Component::mu);
      slope_Vmu[1] = (Dmumu_new - Dmumu) / dmu;

// Compute second slope for Dmumu
      dcoords = dcoords_orig;
      mu_new = dcoords.MomMu() + gambar * slope_Dmumu[0] * dW;
      mu_new += slope_Vmu[0] * dt_local;
      if (fabs(mu_new) > 1.0) {
         if constexpr (HConfig::build_mode == BuildMode::debug) {
            Nabsmugt1 += 1;
         }
         continue;
      };

      dcoords.MomMu() = mu_new;
      dcoords.Time() += dt_local;
      CommonFields(dcoords, dfields);
      diffusion->Stage(dcoords, dfields);
      Dmumu = diffusion->Get(Component::mu);
      slope_Dmumu[1] = sqrt(2.0 * Dmumu);

// Compute third slope for Dmumu
      dcoords = dcoords_orig;
      mu_new = dcoords.MomMu() - slope_Dmumu[0] / (3.0 * gambar) * dW;
      mu_new += slope_Vmu[0] * dt_local;
      if (fabs(mu_new) > 1.0) {
         if constexpr (HConfig::build_mode == BuildMode::debug)
            Nabsmugt1 += 1;
         continue;
      };

      dcoords.MomMu() = mu_new;
      dcoords.Time() += dt_local;
      CommonFields(dcoords, dfields);
      diffusion->Stage(dcoords, dfields);
      Dmumu = diffusion->Get(Component::mu);
      slope_Dmumu[2] = sqrt(2.0 * Dmumu);

// Compute additional fit term
      dcoords = dcoords_orig;
      dmu = sp_small * (dcoords.MomMu() + sp_small < 1.0 ? 1.0 : -1.0);
      dcoords.MomMu() += dmu;
      CommonFields(dcoords, dfields);
      diffusion->Stage(dcoords, dfields);
      Dmumu_new = diffusion->Get(Component::mu);
      dfit = (sqrt(2.0 * Dmumu_new) - slope_Dmumu[0]) / dmu;

// Add all the stuff
      mu_new = dcoords_orig.MomMu() + (0.5 * slope_Dmumu[0] + (slope_Dmumu[1] + 3.0 * Sqr(gambar) * slope_Dmumu[2]) / (2.0 + 6.0 * Sqr(gambar)) ) * dW;
      mu_new += 0.5 * dt_local * (slope_Vmu[0] + slope_Vmu[1]);
      mu_new += 0.5 * slope_Dmumu[0] * dfit * (Sqr(dW) - dt_local);

      if constexpr (HConfig::build_mode == BuildMode::debug) {
         if (fabs(mu_new) > 1.0) Nabsmugt1 += 1;
      }
   } while (fabs(mu_new) > 1.0);

   dcoords.MomMu() = mu_new;
   _coords = TrajectoryCoordinates::Convert(dcoords);
   // todo review: it performs these
//   _coords.MomPerp() = mom_conv[0] * sqrt(1 - Sqr(mu_new));
//   _coords.MomPara() = mom_conv[0] * mu_new;
//   _coords.Vel() = Vel<specie>(_coords.Mom());
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
template <typename Background, typename Diffusion>
void TrajectoryGuidingScatt<Background, Diffusion>::PhysicalStep(void)
{
   constexpr double cfl_pa = cfl_pitchangle;
   double dmumax;
   if constexpr (const_dmumax == TrajectoryOptions::ConstDmumax::constant_dtheta_max) {
      double dthetamax = 2.0*M_PI/180.0;
      dmumax = sqrt(1 - fabs(_coords.MomPara() / _coords.AbsMom())) * dthetamax + 0.5 * fabs(_coords.MomPara() / _coords.AbsMom()) * Sqr(dthetamax);
   }
   else {
      dmumax = 0.02;
   }
   dt_physical = fmin(dt_physical, cfl_pa * Sqr(dmumax) / Dmumu);
   dt_physical = fmin(dt_physical, cfl_pa * dmumax / fabs(Vmu));
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 04/29/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Background, typename Diffusion>
bool TrajectoryGuidingScatt<Background, Diffusion>::Advance(void)
{

// Compute drift and diffusion coefficients (including PA advection term) for physical step computation
   DiffusionCoeff();
   DriftCoeff();

   TrajectoryGuiding::PhysicalStep();
   PhysicalStep();
   dt = fmin(dt_physical, dt_adaptive);
   TimeBoundaryProximityCheck();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// First half of stochastic pitch angle contribution and advection term
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Perform first half of PA scattering
   if constexpr (stochastic_method_mu == TrajectoryOptions::StochasticMethod::Euler) {
      EulerPitchAngleScatt(0);
   }
   else if constexpr (stochastic_method_mu == TrajectoryOptions::StochasticMethod::Milstein) {
      MilsteinPitchAngleScatt(0);
   }
   else if constexpr (stochastic_method_mu == TrajectoryOptions::StochasticMethod::RK2) {
      RK2PitchAngleScatt(0);
   }

// Store position and momentum locally
   StoreLocal();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute the RK slopes.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The commomn fields and "dmax" have been computed at the end of Advance() or in SetStart() before the first step.
// Compute the slopes. The first two components for momentum are always zero for GC (the perpendicular momentum is determined from conservation of magnetic moment).
   Slopes(slope_pos[0], slope_mom[0]);

// If trajectory terminated (or is invalid) while computing slopes, exit advance function with true (step was taken)
   if (RKSlopes()) return true;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Advance the trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If adaptive method error is unacceptable, exit advance function with false (step was not taken)
   if (RKStep()) return false;

   // todo replace bool split_scatt with check whether fraction < 1
   if constexpr (split_scatt) {
//----------------------------------------------------------------------------------------------------------------------------------------------------
// Second half of stochastic pitch angle contribution and advection term
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If an exit spatial boundary was crossed, the fields may no longer be available, so the second half of the scattering cannot be performed. In that case the function should return immediately and the current position will be saved.
      if (SpaceTerminateCheck()) return true;

// Compute diffusion coefficients (including PA advection term)
      CommonFields(_coords, _fields);
      DiffusionCoeff();

// Perform second half of PA scattering
      if constexpr (stochastic_method_mu == TrajectoryOptions::StochasticMethod::Euler) {
         EulerPitchAngleScatt(1);
      }
      else if constexpr (stochastic_method_mu == TrajectoryOptions::StochasticMethod::Milstein) {
         MilsteinPitchAngleScatt(1);
      }
      else if constexpr (stochastic_method_mu == TrajectoryOptions::StochasticMethod::RK2) {
         RK2PitchAngleScatt(1);
      }
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------
   HandleBoundaries();

// If trajectory is not finished (in particular, spatial boundary not crossed), the fields can be computed and momentum corrected
   if (BITS_LOWERED(_status, TRAJ_FINISH)) {
      CommonFields(_coords, _fields);
      MomentumCorrection();
   };

// Add the new point to the trajectory.
   records.Store(_coords);

   return true;
};

};
