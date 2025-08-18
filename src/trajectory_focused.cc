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
template <typename Fields>
TrajectoryFocused<Fields>::TrajectoryFocused(void)
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
template <typename Fields>
TrajectoryFocused<Fields>::TrajectoryFocused(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                 : TrajectoryBase(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
template <typename Fields>
void TrajectoryFocused<Fields>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Magnetic moment is conserved (in the absence of scattering)
   mag_mom = MagneticMoment(traj_mom[0][0] * sqrt(1.0 - Sqr(traj_mom[0][1])), _fields.AbsMag(), specie);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
*/
template <typename Fields>
void TrajectoryFocused<Fields>::DriftCoeff(void)
{
// TODO: improve Field types for these algebraic procedures to avoid explicit casting
   auto U = static_cast<GeoVector>(_fields.Vel());
   auto bhat = static_cast<GeoVector>(_fields.HatMag());
#ifdef TRAJ_FOCUSED_USE_B_DRIFTS
// Compute st2 and ct2
   ct2 = Sqr(_mom[1]);
   st2 = 1.0 - ct2;
// Compute gradient drift
   auto B = static_cast<GeoVector>(_fields.Mag());
   auto absB = static_cast<double>(_fields.AbsMag());
   drift_vel = 0.5 * st2 * (bhat ^ static_cast<GeoVector>(_fields.DelAbsMag())) / absB;
// Add curvature drift. Note that (Bvec * grad)Bvec = [Bvec]^T * [gradBvec].
   drift_vel += ( 0.5 * st2 * bhat * (B * static_cast<GeoVector>(_fields.curlB()))
                + ct2 * (bhat ^ (B * static_cast<GeoMatrix>(_fields.DelMag()))) ) / Sqr(absB);
// Scale by pvc/qB
   drift_vel *= LarmorRadius(_mom[0], _fields.AbsMag(), specie) * _vel[0];
// Add bulk flow and parallel velocities
   drift_vel += U + _vel[0] * _mom[1] * bhat;
#else
   drift_vel = U + _vel[0] * _mom[1] * bhat;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
*/
template <typename Fields>
void TrajectoryFocused<Fields>::PhysicalStep(void)
{
// If the pitch angle is at 90 degrees we only have the perpendicular component of "drift_vel", which may be too small, but can increase by a large (relative) factor during the integration step. For this reason a small fraction of the total speed is added to the characteristic speed.
   dt_physical = cfl_adv_tf * _dmax / (drift_vel.Norm() + drift_safety_tf * _vel[0]);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename Fields>
void TrajectoryFocused<Fields>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   GeoMatrix bhatbhat;
   GeoVector cdUvecdt;
   double bhat_cdUvecdt, bhatbhat_gradUvec;

   DriftCoeff();
   slope_pos_istage = drift_vel;

#ifndef TRAJ_FOCUSED_USE_B_DRIFTS
   st2 = 1.0 - Sqr(_mom[1]);
#endif
// Compute bb : grad U
   GeoVector U = _fields.Vel();
   GeoVector bhat = static_cast<GeoVector>(_fields.HatMag());
   GeoMatrix gradU = _fields.DelVel();
   GeoVector DdtU = _fields.DdtVel();
   bhatbhat.Dyadic(bhat);
   bhatbhat_gradUvec = bhatbhat % gradU;
// Compute 2.0 * b * (convective)dU/dt / v. Note that (Uvec * grad)Uvec = [Uvec]^T * [gradUvec].
// TODO: improve Field types for these algebraic procedures
   GeoVector tmp = U * gradU;
   cdUvecdt = DdtU + tmp;
   bhat_cdUvecdt = 2.0 * bhat * cdUvecdt / _vel[0];

   slope_mom_istage[0] = 0.5 * _mom[0] * ( (3.0 * st2 - 2.0) * bhatbhat_gradUvec 
                                         - st2 * _fields.divU() - _mom[1] * bhat_cdUvecdt);

#if PPERP_METHOD == 1
   slope_mom_istage[1] = 0.5 * st2 * ( _vel[0] * _fields.divbhat() - bhat_cdUvecdt
                                     + _mom[1] * (_fields.divU() - 3.0 * bhatbhat_gradUvec) );
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
template <typename Fields>
bool TrajectoryFocused<Fields>::Advance(void)
{
   return RKAdvance();
};

/*!
\author Juan G Alonso Guzman
\date 08/07/2023
*/
template <typename Fields>
inline void TrajectoryFocused<Fields>::MomentumCorrection(void)
{
// Adjust perp component to conserve magnetic moment
#if PPERP_METHOD == 0
   _mom[1] = sqrt(1.0 - Sqr(PerpMomentum(mag_mom, _fields.AbsMag(), specie) / _mom[0]));
#endif

//Check to enforce |mu| <= 1.0
   if (_mom[1] > 1.0) _mom[1] = 1.0;
   else if (_mom[1] < -1.0) _mom[1] = -1.0;
};

};
