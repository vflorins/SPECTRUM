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
TrajectoryFocused::TrajectoryFocused(void)
                 : TrajectoryBase(traj_name_focused, 0, STATE_NONE, defsize_focused)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
\param[in] presize_in Whether to pre-allocate memory for trajectory arrays
*/
TrajectoryFocused::TrajectoryFocused(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                 : TrajectoryBase(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
void TrajectoryFocused::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Redefine mask
   _spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB | BACKGROUND_gradU | BACKGROUND_dUdt;
   spdata0._mask = BACKGROUND_ALL | BACKGROUND_gradB | BACKGROUND_gradU | BACKGROUND_dUdt;

// Magnetic moment is conserved (in the absence of scattering)
   mag_mom = MagneticMoment(traj_mom[0][0] * sqrt(1.0 - Sqr(traj_mom[0][1])), _spdata.Bmag, specie);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
*/
void TrajectoryFocused::DriftCoeff(void)
{
#ifdef TRAJ_FOCUSED_USE_B_DRIFTS
// Compute st2 and ct2
   ct2 = Sqr(_mom[1]);
   st2 = 1.0 - ct2;
// Compute gradient drift
   drift_vel = 0.5 * st2 * (_spdata.bhat ^ _spdata.gradBmag()) / _spdata.Bmag;
// Add curvature drift. Note that (Bvec * grad)Bvec = [gradBvec] * [Bvec].
   drift_vel += ( 0.5 * st2 * _spdata.bhat * (_spdata.Bvec * _spdata.curlB())
                + ct2 * (_spdata.bhat ^ (_spdata.gradBvec * _spdata.Bvec)) ) / Sqr(_spdata.Bmag);
// Scale by pvc/qB
   drift_vel *= LarmorRadius(_mom[0], _spdata.Bmag, specie) * _vel[0];
// Add bulk flow and parallel velocities
   drift_vel += _spdata.Uvec + _vel[0] * _mom[1] * _spdata.bhat;
#else
   drift_vel = _spdata.Uvec + _vel[0] * _mom[1] * _spdata.bhat;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
*/
void TrajectoryFocused::PhysicalStep(void)
{
// If the pitch angle is at 90 degrees we only have the perpendicular component of "drift_vel", which may be too small, but can increase by a large (relative) factor during the integration step. For this reason a small fraction of the total speed is added to the characteristic speed.
   dt_physical = cfl_adv_tf * _spdata.dmax / (drift_vel.Norm() + drift_safety_tf * _vel[0]);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
void TrajectoryFocused::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   GeoMatrix bhatbhat;
   GeoVector cdUvecdt;
   double divbhat, bhat_cdUvecdt, bhatbhat_gradUvec, divUvec;

   DriftCoeff();
   slope_pos_istage = drift_vel;

#ifndef TRAJ_FOCUSED_USE_B_DRIFTS
   st2 = 1.0 - Sqr(_mom[1]);
#endif
// Compute div b
   divbhat = (_spdata.divB() - _spdata.bhat * _spdata.gradBmag()) / _spdata.Bmag;
// Compute bb : grad U
   bhatbhat.Dyadic(_spdata.bhat);
   bhatbhat_gradUvec = bhatbhat % _spdata.gradUvec;
// Compute div U
   divUvec = _spdata.divU();
// Compute 2.0 * b * (convective)dU/dt / v. Note that (Uvec * grad)Uvec = [gradUvec] * [Uvec].
   cdUvecdt = _spdata.dUvecdt + (_spdata.gradUvec * _spdata.Uvec);
   bhat_cdUvecdt = 2.0 * _spdata.bhat * cdUvecdt / _vel[0];

   slope_mom_istage[0] = 0.5 * _mom[0] * ( (3.0 * st2 - 2.0) * bhatbhat_gradUvec 
                                         - st2 * divUvec - _mom[1] * bhat_cdUvecdt);

#if PPERP_METHOD == 1
   slope_mom_istage[1] = 0.5 * st2 * ( _vel[0] * divbhat - bhat_cdUvecdt 
                                     + _mom[1] * (divUvec - 3.0 * bhatbhat_gradUvec) );
#else
   slope_mom_istage[1] = 0.0;
#endif
   
   slope_mom_istage[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryFocused::Advance(void)
{
   return RKAdvance();
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
inline void TrajectoryFocused::MomentumCorrection(void)
{
// Adjust perp component to conserve magnetic moment
#if PPERP_METHOD == 0
   _mom[1] = sqrt(1.0 - Sqr(PerpMomentum(mag_mom, _spdata.Bmag, specie) / _mom[0]));
#endif

//Check to enforce |mu| <= 1.0
   if(_mom[1] > 1.0) _mom[1] = 1.0;
   else if(_mom[1] < -1.0) _mom[1] = -1.0;
};

};
