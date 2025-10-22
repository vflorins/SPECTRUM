/*!
\file trajectory_focused.cc
\brief Defines a class for trajectory based on the focused transport equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_focused.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFocused methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
template <typename Background, typename Diffusion>
TrajectoryFocused<Background, Diffusion>::TrajectoryFocused(void)
                 : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Background, typename Diffusion>
TrajectoryFocused<Background, Diffusion>::TrajectoryFocused(const std::string& name_in, uint16_t status_in)
                 : TrajectoryBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
template <typename Background, typename Diffusion>
void TrajectoryFocused<Background, Diffusion>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Magnetic moment is conserved (in the absence of scattering)
   mag_mom = MagneticMoment<specie>(_coords.Mom()[0] * sqrt(1.0 - Sqr(_coords.Mom()[1])), _fields.AbsMag());
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
*/
template <typename Background, typename Diffusion>
void TrajectoryFocused<Background, Diffusion>::DriftCoeff(void)
{
   auto U = _fields.Vel();
   auto bhat = _fields.HatMag();
   if constexpr (HConfig::use_B_drifts) {
// Compute st2 and ct2
      ct2 = Sqr(_coords.Mom()[1]);
      st2 = 1.0 - ct2;
// Compute gradient drift
      auto B = _fields.Mag();
      auto absB = _fields.AbsMag();
      drift_vel = 0.5 * st2 * (bhat ^ _fields.DelAbsMag()) / absB;
// Add curvature drift. Note that (Bvec * grad)Bvec = [Bvec]^T * [gradBvec].
      drift_vel += ( 0.5 * st2 * bhat * (B * _fields.CurlMag())
                     + ct2 * (bhat ^ (B * _fields.DelMag())) ) / Sqr(absB);
// Scale by pvc/qB
      drift_vel *= LarmorRadius<specie>(_coords.Mom()[0], _fields.AbsMag()) * _coords.Vel()[0];
// Add bulk flow and parallel velocities
      drift_vel += U + _coords.Vel()[0] * _coords.Mom()[1] * bhat;
   }
   else {
      drift_vel = U + _coords.Vel()[0] * _coords.Mom()[1] * bhat;
   }
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
*/
template <typename Background, typename Diffusion>
void TrajectoryFocused<Background, Diffusion>::PhysicalStep(void)
{
   constexpr double cfl_adv = HConfig::cfl_advection;
   constexpr double drift_safety = HConfig::drift_safety;
// If the pitch angle is at 90 degrees we only have the perpendicular component of "drift_vel", which may be too small, but can increase by a large (relative) factor during the integration step. For this reason a small fraction of the total speed is added to the characteristic speed.
   dt_physical = cfl_adv * _dmax / (drift_vel.Norm() + drift_safety * _coords.Vel()[0]);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename Background, typename Diffusion>
void TrajectoryFocused<Background, Diffusion>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   GeoMatrix bhatbhat;
   GeoVector cdUvecdt;
   double bhat_cdUvecdt, bhatbhat_gradUvec;

   DriftCoeff();
   slope_pos_istage = drift_vel;

   if constexpr (!HConfig::use_B_drifts) {
      st2 = 1.0 - Sqr(_coords.Mom()[1]);
   }
// Compute bb : grad U
   GeoVector U = _fields.Vel();
   GeoVector bhat = _fields.HatMag();
   GeoMatrix gradU = _fields.DelVel();
   GeoVector DotU = _fields.DotVel();
   bhatbhat.Dyadic(bhat);
   bhatbhat_gradUvec = bhatbhat % gradU;
// Compute 2.0 * b * (convective)dU/dt / v. Note that (Uvec * grad)Uvec = [Uvec]^T * [gradUvec].
   GeoVector tmp = U * gradU;
   cdUvecdt = DotU + tmp;
   bhat_cdUvecdt = 2.0 * bhat * cdUvecdt / _coords.Vel()[0];

   slope_mom_istage[0] = 0.5 * _coords.Mom()[0] * ( (3.0 * st2 - 2.0) * bhatbhat_gradUvec
                                         - st2 * _fields.DivFluv() - _coords.Mom()[1] * bhat_cdUvecdt);

   if constexpr (HConfig::pperp_method == TrajectoryOptions::PPerpMethod::mag_moment_conservation) {
      slope_mom_istage[1] = 0.0;
   }
   else if constexpr (HConfig::pperp_method == TrajectoryOptions::PPerpMethod::scheme) {
      slope_mom_istage[1] = 0.5 * st2 * ( _coords.Vel()[0] * _fields.DivHatMag() - bhat_cdUvecdt
                                          + _coords.Mom()[1] * (_fields.DivFluv() - 3.0 * bhatbhat_gradUvec) );
   }
   
   slope_mom_istage[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Background, typename Diffusion>
bool TrajectoryFocused<Background, Diffusion>::Advance(void)
{
   return RKAdvance();
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
template <typename Background, typename Diffusion>
inline void TrajectoryFocused<Background, Diffusion>::MomentumCorrection(void)
{
   if constexpr (HConfig::pperp_method == TrajectoryOptions::PPerpMethod::mag_moment_conservation) {
// Adjust perp component to conserve magnetic moment
      _coords.Mom()[1] = sqrt(1.0 - Sqr(PerpMomentum<specie>(mag_mom, _fields.AbsMag()) / _coords.Mom()[0]));
   }

// Check to enforce |mu| <= 1.0
   if (_coords.Mom()[1] > 1.0) _coords.Mom()[1] = 1.0;
   else if (_coords.Mom()[1] < -1.0) _coords.Mom()[1] = -1.0;
};

};
