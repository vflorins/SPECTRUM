/*!
\file trajectory_parker.cc
\brief Defines a class for trajectory based on the Parker Transport Equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_parker.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryParker methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Background, typename Diffusion>
TrajectoryParker<Background, Diffusion>::TrajectoryParker(void)
      : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Background, typename Diffusion>
TrajectoryParker<Background, Diffusion>::TrajectoryParker(const std::string& name_in, uint16_t status_in)
      : TrajectoryBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Background, typename Diffusion>
bool TrajectoryParker<Background, Diffusion>::IsSimulationReady(void) const
{
   if (!TrajectoryBase::IsSimulationReady()) return false;

// A diffusion object is required
   if (diffusion == nullptr) return false;
   else if (BITS_LOWERED(diffusion->GetStatus(), STATE_SETUP_COMPLETE)) return false;
   return true;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::FieldAlignedFrame(void)
{
   fa_basis[2] = _fields.HatMag();
   fa_basis[0] = GetSecondUnitVec(_fields.HatMag());
   fa_basis[1] = fa_basis[2] ^ fa_basis[0];
};

/*!
\author Juan G Alonso Guzman
\date 03/11/2024
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::DiffusionCoeff(void)
try {
   int i,j;
   // todo review impl for parker coords -> diffusion coords (any)
   /*
    *
    * REMINDER TO SELF: the dcoords are always (t, x, y, z, p, mu) !
    *
    */
   auto dcoords = DiffusionCoordinates::Convert(_coords);
   DiffusionFields dfields = DiffusionCoordinates::Get(_fields);
   CommonFields<DiffusionCoordinates, DiffusionFields, DiffusionFieldsRemainder>(dcoords, dfields);
   diffusion->Stage(dcoords, dfields);

   if constexpr (HConfig::divk_method == TrajectoryOptions::DivkMethod::direct) {
// Compute using Diffusion Fields type, this only compute fields needed by diffusion
      double Kperp_forw, Kperp_back, Kpara_forw, Kpara_back, Kappa_forw, Kappa_back;
      double delta = fmin(LarmorRadius<specie>(_coords.AbsMom(), _fields.AbsMag()), _dmax);

// Compute perpendicular and parallel diffusion coefficients and diffusion tensor.
      Kperp = diffusion->Get(Component::perp);
      Kpara = diffusion->Get(Component::para);

// Loop over dimensions to find derivatives of Kappa.
      divK = gv_zeros;
// TODO: check if CommonFields returns a STATE_INVALID flag and revert to forward/backward (1st order) FD.
      for (j = 0; j < 3; j++) {
// Forward evaluation
         dcoords.Pos()[j] += delta;
         //\\ The position has changed, so the CommonFields method must be called for all dfields.
         CommonFields(dcoords, dfields);
         //> Stage diffusion again.
         diffusion->Stage(dcoords, dfields);
         Kperp_forw = diffusion->Get(Component::perp);
         Kpara_forw = diffusion->Get(Component::para);
         for (i = 0; i < 3; i++) {
            Kappa_forw = Kperp_forw * (i == j ? 1.0 : 0.0) + (Kpara_forw - Kperp_forw) * dfields.HatMag()[j] * dfields.HatMag()[i];
         };
// Backward evaluation
         dcoords.Pos()[j] -= 2.0 * delta;
         CommonFields(dcoords, dfields);
         diffusion->Stage(dcoords, dfields);
         Kperp_back = diffusion->Get(Component::perp);
         Kpara_back = diffusion->Get(Component::para);
         for (i = 0; i < 3; i++) {
            Kappa_back = Kperp_back * (i == j ? 1.0 : 0.0) + (Kpara_back - Kperp_back) * dfields.HatMag()[j] * dfields.HatMag()[i];
            divK[i] += 0.5 * (Kappa_forw - Kappa_back) / delta;
         };
         dcoords.Pos()[j] += delta;
      };
   }
   else {
      GeoVector gradKpara, gradKperp;

// Compute Kperp and grad(Kperp)
      Kperp = diffusion->Get(Component::perp);
      auto ddata = background->GetDerivativeData();
      gradKperp[0] = diffusion->GetDirectionalDerivative(Component::perp, 0, ddata);
      gradKperp[1] = diffusion->GetDirectionalDerivative(Component::perp, 1, ddata);
      gradKperp[2] = diffusion->GetDirectionalDerivative(Component::perp, 2, ddata);

// Compute Kpara and grad(Kpara)
      Kpara = diffusion->Get(Component::para);
      ddata = background->GetDerivativeData();
      gradKpara[0] = diffusion->GetDirectionalDerivative(Component::para, 0, ddata);
      gradKpara[1] = diffusion->GetDirectionalDerivative(Component::para, 1, ddata);
      gradKpara[2] = diffusion->GetDirectionalDerivative(Component::para, 2, ddata);

// Assemble diffusion tensor
      GeoVector HatMag = _fields.HatMag();
      GeoMatrix DyadHatMag = Dyadic(HatMag);
      divK = gradKperp + DyadHatMag * (gradKpara - gradKperp) + (Kpara - Kperp) * (_fields.DivHatMag() * HatMag + HatMag * _fields.DelHatMag());
   }
}

catch(ExFieldError& exception) {
   // PrintError(__FILE__, __LINE__, "Error in field increment evaluation");
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Juan G Alonso Guzman
\date 03/11/2024
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::EulerDiffSlopes(void)
{
   double dWx, dWy, dWz;
   FieldAlignedFrame();

// Generate stochastic factors
   dWx = sqrt(dt) * rng->GetNormal();
   dWy = sqrt(dt) * rng->GetNormal();
   dWz = sqrt(dt) * rng->GetNormal();

// Recompute Kperp and Kpara at the beginning of the step
   auto dcoords = DiffusionCoordinates::Convert(_coords);
   DiffusionFields dfields = DiffusionCoordinates::Get(_fields);
   CommonFields<DiffusionCoordinates, DiffusionFields, DiffusionFieldsRemainder>(dcoords, dfields);
   diffusion->Stage(dcoords, dfields);
   Kperp = diffusion->Get(Component::perp);
   Kpara = diffusion->Get(Component::para);

// Compute random displacement
   dr_perp[0] = sqrt(2.0 * Kperp) * dWx;
   dr_perp[1] = sqrt(2.0 * Kperp) * dWy;
   dr_perp[2] = sqrt(2.0 * Kpara) * dWz;

// Convert to global coordinates
   dr_perp.ChangeFromBasis(fa_basis);
};

//TODO: Implement Milstein and RK2 schemes

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\note "Bvec", "Bmag", and "bhat" must be available, so "CommonFields()" must be called prior to this function.
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::DriftCoeff(void)
{
   if constexpr (HConfig::use_B_drifts == TrajectoryOptions::UseBDrifts::gradient_curvature) {
      // Compute |B|*curl(b/|B|)
      drift_vel = (_fields.CurlMag() - 2.0 * (_fields.DelAbsMag() ^ _fields.HatMag())) / _fields.AbsMag();
// Scale by pvc/3q|B| = r_L*v/3
      drift_vel *= LarmorRadius<specie>(_coords.Mom()[0], _fields.AbsMag()) * _coords.Vel()[0] / 3.0;
// Scale magnitude to an upper limit of v/2 if necessary.
      if (drift_vel.Norm() > 0.5 * _coords.Vel()[0]) {
         drift_vel.Normalize();
         drift_vel *= 0.5 * _coords.Vel()[0];
      };
// Add bulk flow velocity
      drift_vel += _fields.Vel();
   }
   else {
      drift_vel = _fields.Vel();
   }
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::PhysicalStep(void)
{
   constexpr double cfl_adv = HConfig::cfl_advection;
   constexpr double cfl_dif = HConfig::cfl_diffusion;
   constexpr double cfl_acc = HConfig::cfl_acceleration;
   constexpr double dlnpmax = HConfig::dlnp_max;
   if constexpr (HConfig::time_flow == TimeFlow::forward) {
      dt_physical = cfl_adv * _dmax / (drift_vel + divK).Norm();
   }
   else {
      dt_physical = cfl_adv * _dmax / (drift_vel - divK).Norm();
   }
   dt_physical = fmin(dt_physical, cfl_dif * Sqr(_dmax) / fmax(Kperp, Kpara));
   dt_physical = fmin(dt_physical, cfl_acc * 3.0 * dlnpmax / fabs(_fields.DivVel()));
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename Background, typename Diffusion>
void TrajectoryParker<Background, Diffusion>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   DriftCoeff();
   DiffusionCoeff();

// Position slopes
   slope_pos_istage = drift_vel;


   if constexpr (HConfig::time_flow == TimeFlow::forward) {
      slope_pos_istage += divK;
   }
   else {
      slope_pos_istage -= divK;
   }

// Momentum slopes
   slope_mom_istage[0] = -_coords.Mom()[0] * _fields.DivVel() / 3.0;
   slope_mom_istage[1] = 0.0;
   slope_mom_istage[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Background, typename Diffusion>
bool TrajectoryParker<Background, Diffusion>::Advance(void)
{
// store locally
   StoreLocal();

// The commomn fields and "dmax" have been computed at the end of Advance() or in SetStart() before the first step.
// Compute the slopes. The first two components for momentum are always zero for GC (the perpendicular momentum is determined from conservation of magnetic moment).
   Slopes(slope_pos[0], slope_mom[0]);

   PhysicalStep();
   dt = fmin(dt_physical, dt_adaptive);
   TimeBoundaryProximityCheck();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute the RK slopes.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Stochastic RK slopes
// TODO: Add other (higher order) options for the stochastic step (e.g. Milstein or RK2)
   if constexpr (HConfig::stochastic_method == TrajectoryOptions::StochasticMethod::Euler) {
      EulerDiffSlopes();
   }
   else if constexpr (HConfig::stochastic_method == TrajectoryOptions::StochasticMethod::Milstein) {
      dr_perp = gv_zeros;
   }
   else if constexpr (HConfig::stochastic_method == TrajectoryOptions::StochasticMethod::RK2) {
      dr_perp = gv_zeros;
   }

// If trajectory terminated (or is invalid) while computing slopes, exit advance function with true (step was taken)
   if (RKSlopes()) return true;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Advance trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If adaptive method error is unacceptable, exit advance function with false (step was not taken)
   if (RKStep()) return false;

// Stochastic displacement
   _coords.Pos() += dr_perp;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------
   HandleBoundaries();

// If trajectory is not finished (in particular, spatial boundary not crossed), the fields can be computed
   if (BITS_LOWERED(_status, TRAJ_FINISH)) CommonFields(_coords, _fields);

// Add the new point to the trajectory.
   records.Store(_coords);

   return true;
};

};
