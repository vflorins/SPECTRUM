/*!
\file trajectory_guiding.cc
\brief Defines a class for trajectory based on relativistic guiding center equations
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_guiding.hh"
#include "common/print_warn.hh"

// debug
#include <iostream>
using std::cout; using std::endl;

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuiding methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
template <typename HConfig>
TrajectoryGuiding<HConfig>::TrajectoryGuiding(void)
      : TrajectoryBase(traj_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
TrajectoryGuiding<HConfig>::TrajectoryGuiding(const std::string_view& name_in, status_t status_in)
      : TrajectoryBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
void TrajectoryGuiding<HConfig>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Magnetic moment is conserved (in the absence of scattering)
   mag_mom = MagneticMoment<specie>(_coords.MomPerp(), _fields.AbsMag());
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 10/12/2025
\note Magnetic field data must be available, so "CommonFields()" must be called prior to this function.

Ref: Tao, X., Chan, A. A., and Brizard, A. J., Hamiltonian theory of adiabatic motion of relativistic charged particles, Phys. Plasmas, v. 14, p. 09107 (2007).
*/
template <typename HConfig>
void TrajectoryGuiding<HConfig>::ModifiedFields(void)
try {
auto MP = _coords.MomPerp();
auto AM = _fields.AbsMag();

// Modified fields
   auto rL = LarmorRadius<specie>(MP, AM);
   auto rR = LarmorRadius<specie>(_coords.MomPara(), _fields.AbsMag());
   Evec_star = _fields.Ele();
   Evec_star = Evec_star - rR * _fields.AbsMag() * _fields.DotHatMag() / c_code;
   Evec_star = Evec_star - rL * _coords.VelPerp() * _fields.DelAbsMag() / (2.0 * c_code);
   Bvec_star = _fields.Mag();
   Bvec_star = Bvec_star + rR * _fields.AbsMag() * _fields.CurlHatMag();
}

catch(ExFieldError& exception) {
//   PrintError(__FILE__, __LINE__, exception.what());
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Vladimir Florinski
\date 02/11/2022
*/
template <typename HConfig>
void TrajectoryGuiding<HConfig>::DriftCoeff(void)
{
   ModifiedFields();

   drift_vel = (_coords.VelPara() * Bvec_star + c_code * (Evec_star ^ _fields.HatMag())) / (Bvec_star * _fields.HatMag());
   cout << "VelPara, VelPerp " << _coords.VelPara() << ", " << _coords.VelPerp()<< endl;
};

/*!
\author Vladimir Florinski
\date 02/10/2022
*/
template <typename HConfig>
void TrajectoryGuiding<HConfig>::PhysicalStep(void)
{
// If the pitch angle is at 90 degrees we only have the perpendicular component of "drift_vel", which may be too small, but can increase by a large (relative) factor during the integration step. For this reason a small fraction of the total speed is added to the characteristic speed.
   constexpr double cfl_adv = Config::cfl_advection;
   constexpr double drift_safety = Config::drift_safety;
   dt_physical = cfl_adv * _dmax/(drift_vel.Norm() + drift_safety * _coords.Vel().Norm());
//   dt_physical = cfl_adv * _dmax / (drift_vel.Norm() + drift_safety * _coords.Vel().Norm());
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 02/11/2022
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename HConfig>
void TrajectoryGuiding<HConfig>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   DriftCoeff();
   slope_pos_istage = drift_vel;

   if constexpr (Config::pperp_method == TrajectoryOptions::PPerpMethod::mag_moment_conservation) {
      slope_mom_istage[0] = 0.0;
   }
   else if constexpr (Config::pperp_method == TrajectoryOptions::PPerpMethod::scheme) {
      slope_mom_istage[0] = 0.5 * _coords.MomPerp() / _fields.AbsMag() * (_fields.DotAbsMag() + drift_vel * _fields.DelAbsMag());
   }

   slope_mom_istage[1] = 0.0;
   slope_mom_istage[2] = specie.q * (Evec_star * Bvec_star) / (Bvec_star * _fields.HatMag());
};

/*!
\author Juan G Alonso Guzman
\date 04/19/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename HConfig>
bool TrajectoryGuiding<HConfig>::Advance(void)
{
   return RKAdvance();
};

/*!
\author Juan G Alonso Guzman
\date 04/19/2022
*/
template <typename HConfig>
inline void TrajectoryGuiding<HConfig>::MomentumCorrection(void)
{
   if constexpr (Config::pperp_method == TrajectoryOptions::PPerpMethod::mag_moment_conservation) {
// Adjust perp component to conserve magnetic moment
      _coords.MomPerp('w') = PerpMomentum<specie>(mag_mom, _fields.AbsMag());
   }
};

};
