/*!
\file trajectory_guiding.cc
\brief Defines a class for trajectory based on relativistic guiding center equations
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_guiding.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuiding methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
TrajectoryGuiding::TrajectoryGuiding(void)
                 : TrajectoryBase(traj_name_guiding, 0, STATE_NONE, defsize_guiding)
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
TrajectoryGuiding::TrajectoryGuiding(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                 : TrajectoryBase(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
void TrajectoryGuiding::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Redefine mask
   _spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB | BACKGROUND_dBdt;
   spdata0._mask = BACKGROUND_ALL | BACKGROUND_gradB | BACKGROUND_dBdt;

// Magnetic moment is conserved (in the absence of scattering)
   mag_mom = MagneticMoment(traj_mom[0][0], _spdata.Bmag, specie);
};

/*!
\author Vladimir Florinski
\date 03/24/2022
\note "Bvec", "Bmag", and "bhat" must be available, so "CommonFields()" must be called prior to this function.

Ref: Tao, X., Chan, A. A., and Brizard, A. J., Hamiltonian theory of adiabatic motion of relativistic charged particles, Phys. Plasmas, v. 14, p. 09107 (2007).
*/
void TrajectoryGuiding::ModifiedFields(void)
try {
   double rL, rR;

// Compute derived quantities from vector gradient of Bvec
   gradBmag = _spdata.gradBmag();
   curlbhat = (_spdata.curlB() - (gradBmag ^ _spdata.bhat)) / _spdata.Bmag;
   dBmagdt =  _spdata.dBvecdt * _spdata.bhat;
   dbhatdt = (_spdata.dBvecdt - (dBmagdt * _spdata.bhat)) / _spdata.Bmag;

// Modified fields
   rL = LarmorRadius(_mom[0], _spdata.Bmag, specie);
   rR = LarmorRadius(_mom[2], _spdata.Bmag, specie);
   Evec_star = _spdata.Evec - rR * _spdata.Bmag * dbhatdt / c_code - rL * _vel[0] * gradBmag / (2.0 * c_code);
   Bvec_star = _spdata.Bvec + rR * _spdata.Bmag * curlbhat;
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
void TrajectoryGuiding::DriftCoeff(void)
{
   ModifiedFields();
   drift_vel = (_vel[2] * Bvec_star + c_code * (Evec_star ^ _spdata.bhat)) / (Bvec_star * _spdata.bhat);
};

/*!
\author Vladimir Florinski
\date 02/10/2022
*/
void TrajectoryGuiding::PhysicalStep(void)
{
// If the pitch angle is at 90 degrees we only have the perpendicular component of "drift_vel", which may be too small, but can increase by a large (relative) factor during the integration step. For this reason a small fraction of the total speed is added to the characteristic speed.
   dt_physical = cfl_adv_tg * _spdata.dmax / (drift_vel.Norm() + drift_safety_tg * _vel.Norm());
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 02/11/2022
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
void TrajectoryGuiding::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   DriftCoeff();
   slope_pos_istage = drift_vel;

#if PPERP_METHOD == 1
   slope_mom_istage[0] = _mom[0] / 2.0 / _spdata.Bmag * (dBmagdt + drift_vel * gradBmag);
#else
   slope_mom_istage[0] = 0.0;
#endif

   slope_mom_istage[1] = 0.0;
   slope_mom_istage[2] = q * (Evec_star * Bvec_star) / (Bvec_star * _spdata.bhat);
};

/*!
\author Juan G Alonso Guzman
\date 04/19/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryGuiding::Advance(void)
{
   return RKAdvance();
};

/*!
\author Juan G Alonso Guzman
\date 04/19/2022
*/
inline void TrajectoryGuiding::MomentumCorrection(void)
{
// Adjust perp component to conserve magnetic moment
#if PPERP_METHOD == 0
   _mom[0] = PerpMomentum(mag_mom, _spdata.Bmag, specie);
#endif
};

};
