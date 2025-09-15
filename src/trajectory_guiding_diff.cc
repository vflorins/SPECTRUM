/*!
\file trajectory_guiding_diff.cc
\brief Defines a class for trajectory based on guiding center equations with perpendicular diffusion
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "trajectory_guiding_diff.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingDiff methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
*/
template <typename Fields>
TrajectoryGuidingDiff<Fields>::TrajectoryGuidingDiff(void)
                     : TrajectoryGuidingBase(traj_name_guidingdiff, 0, STATE_NONE, defsize_guidingdiff)
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
template <typename Fields>
TrajectoryGuidingDiff<Fields>::TrajectoryGuidingDiff(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                     : TrajectoryGuidingBase(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
*/
template <typename Fields>
bool TrajectoryGuidingDiff<Fields>::IsSimulationReady(void) const
{
   if(!TrajectoryBase::IsSimulationReady()) return false;

// A diffusion object is required
   if(diffusion == nullptr) return false;
   else if(BITS_LOWERED(diffusion->GetStatus(), STATE_SETUP_COMPLETE)) return false;
   return true;
};

/*!
\author Vladimir Florinski
\date 04/22/2022
*/
template <typename Fields>
void TrajectoryGuidingDiff<Fields>::FieldAlignedFrame(void)
{
   auto bhat = _fields.HatMag();
   fa_basis[2] = bhat;
   fa_basis[0] = GetSecondUnitVec(bhat);
   fa_basis[1] = fa_basis[2] ^ fa_basis[0];
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/16/2025
*/
template <typename Fields>
void TrajectoryGuidingDiff<Fields>::DiffusionCoeff(void)
{
   try {
      GeoVector gradDperp;
      FieldAlignedFrame();
      Dperp = diffusion->GetComponent(0, _t, _pos, ConvertMomentum(), _fields);

// Compute gradient of Dperp
      auto ddata = background->GetDerivativeData();
      gradDperp[0] = diffusion->GetDirectionalDerivative(0, ddata);
      gradDperp[1] = diffusion->GetDirectionalDerivative(1, ddata);
      gradDperp[2] = diffusion->GetDirectionalDerivative(2, ddata);

// Find components of Vperp in field aligned frame
      Vperp[0] = gradDperp * fa_basis[0];
      Vperp[1] = gradDperp * fa_basis[1];
      Vperp[2] = 0.0;

// Convert "Vperp" to global coordinates
      Vperp.ChangeFromBasis(fa_basis);
   }

   catch (ExFieldError &exception) {
      // PrintError(__FILE__, __LINE__, "Error in field increment evaluation");
      RAISE_BITS(_status, TRAJ_DISCARD);
      throw;
   };
}

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/16/2022
*/
template <typename Fields>
void TrajectoryGuidingDiff<Fields>::EulerPerpDiffSlopes(void)
{
   double dWx, dWy;
   GeoVector mom_conv = ConvertMomentum();

// Generate stochastic factors
   dWx = sqrt(dt) * rng->GetNormal();
   dWy = sqrt(dt) * rng->GetNormal();

// Recompute Dperp at the beginning of the step
   Dperp = diffusion->GetComponent(0, _t, _pos, mom_conv, _fields);

// Compute random displacement
   dr_perp[0] = sqrt(2.0 * Dperp) * dWx;
   dr_perp[1] = sqrt(2.0 * Dperp) * dWy;
   dr_perp[2] = 0.0;

// Convert to global coordinates
   dr_perp.ChangeFromBasis(fa_basis);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/05/2023
*/
template <typename Fields>
void TrajectoryGuidingDiff<Fields>::MilsteinPerpDiffSlopes(void)
{
   double dWx, dWy, Vxy, random, Dperp_new;
   double dbdx, dbdy, dx, dy, slope_Dperp;
   GeoVector mom_conv = ConvertMomentum(), pos_new, xhat, yhat;
   Fields fields_new;

// Generate stochastic factors
   dWx = sqrt(dt) * rng->GetNormal();
   dWy = sqrt(dt) * rng->GetNormal();

// Recompute Dperp at the beginning of the step
   Dperp = diffusion->GetComponent(0, _t, _pos, mom_conv, _fields);
   slope_Dperp = sqrt(2.0 * Dperp);

// Compute random displacement
   dr_perp[0] = slope_Dperp * dWx;
   dr_perp[1] = slope_Dperp * dWy;
   dr_perp[2] = 0.0;

// Calculate derivatives of sqrt(2.0 * Dperp) = b_11 = b_22. Note that b_12 = b_21 = 0.
   xhat = gv_nx;
   xhat.ChangeFromBasis(fa_basis);
   dx = background->GetSafeIncr(xhat);
   pos_new = _pos + dx * xhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   dbdx = (sqrt(2.0 * Dperp_new) - slope_Dperp) / dx;

   yhat = gv_ny;
   yhat.ChangeFromBasis(fa_basis);
   dy = background->GetSafeIncr(yhat);
   pos_new = _pos + dy * yhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   dbdy = (sqrt(2.0 * Dperp_new) - slope_Dperp) / dy;

// Add extra Miltein terms to Euler step
   random = rng->GetUniform();
   Vxy = (random < 0.5 ? dt : -dt);

   dr_perp[0] += 0.5 * slope_Dperp * (dbdx * (Sqr(dWx) - dt) + dbdy * (dWx * dWy - Vxy));
   dr_perp[1] += 0.5 * slope_Dperp * (dbdy * (Sqr(dWy) - dt) + dbdx * (dWx * dWy + Vxy));

// Convert to global coordinates
   dr_perp.ChangeFromBasis(fa_basis);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/16/2022
*/
template <typename Fields>
bool TrajectoryGuidingDiff<Fields>::RK2PerpDiffSlopes(void)
{
   double dWx, dWy, Rx, Ry, Vxy, random;
   double dbdx, dbdy, d2bdx2, d2bdy2, d2bdxdy, dx, dy;
   double slope_Dperp[3], Dperp_newx, Dperp_newy, Dperp_new;
   double gambar = 1.0 / sqrt(3.0);
   GeoVector mom_conv = ConvertMomentum(), dpos, pos_new, xhat, yhat;
   Fields fields_new;

// Generate stochastic factors
   dWx = sqrt(dt) * rng->GetNormal();
   dWy = sqrt(dt) * rng->GetNormal();

// Compute first contribution to stochastic slope
   Dperp = diffusion->GetComponent(0, _t, _pos, mom_conv, _fields);
   slope_Dperp[0] = sqrt(2.0 * Dperp);

// Compute displacement for slope calculation
   dpos[0] = slope_Dperp[0] * dWx;
   dpos[1] = slope_Dperp[0] * dWy;
   dpos[2] = 0.0;
   dpos.ChangeFromBasis(fa_basis);

// Compute second contribution to stochastic slope
   _pos += gambar * dpos;
   if(SpaceTerminateCheck()) return true;

   CommonFields(_t + dt, _pos, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t + dt, _pos, mom_conv, fields_new);
   slope_Dperp[1] = sqrt(2.0 * Dperp_new);
   _pos = local_pos;

// Compute third contribution to stochastic slope
   _pos -= dpos / (3.0 * gambar);
   if(SpaceTerminateCheck()) return true;

   CommonFields(_t + dt, _pos, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t + dt, _pos, mom_conv, fields_new);
   slope_Dperp[2] = sqrt(2.0 * Dperp_new);
   _pos = local_pos;

// Compute additional fit term (this is the most computationally intensive part)

// first and second x derivative
   xhat = gv_nx;
   xhat.ChangeFromBasis(fa_basis);
   dx = background->GetSafeIncr(xhat);
   pos_new = _pos + dx * xhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   dx *= 0.5;
   pos_new = _pos + dx * xhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_newx = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   dbdx = (sqrt(2.0 * Dperp_newx) - slope_Dperp[0]) / dx;
   d2bdx2 = (sqrt(2.0 * Dperp_new) - 2.0 * sqrt(2.0 * Dperp_newx) + slope_Dperp[0]) / Sqr(dx);

// first and second y derivative
   yhat = gv_ny;
   yhat.ChangeFromBasis(fa_basis);
   dy = background->GetSafeIncr(yhat);
   pos_new = _pos + dy * yhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   dy *= 0.5;
   pos_new = _pos + dy * yhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_newy = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   dbdy = (sqrt(2.0 * Dperp_newy) - slope_Dperp[0]) / dy;
   d2bdy2 = (sqrt(2.0 * Dperp_new) - 2.0 * sqrt(2.0 * Dperp_newy) + slope_Dperp[0]) / Sqr(dy);

// mixed second derivative
   pos_new = _pos + dx * xhat + dy * yhat;
   CommonFields(_t, pos_new, mom_conv, fields_new);
   Dperp_new = diffusion->GetComponent(0, _t, pos_new, mom_conv, fields_new);
   d2bdxdy = (sqrt(2.0 * Dperp_new) - sqrt(2.0 * Dperp_newx) - sqrt(2.0 * Dperp_newy) + slope_Dperp[0]) / (dx * dy);

// combination
   random = rng->GetUniform();
   Vxy = (random < 0.5 ? dt : -dt);
   Rx = 0.5 * slope_Dperp[0] * (dbdx * (Sqr(dWx) - dt) + dbdy * (dWx * dWy - Vxy))
         + dt * Sqr(slope_Dperp[0]) * (d2bdy2 * dWx - d2bdxdy * dWy) / 6.0;
   Ry = 0.5 * slope_Dperp[0] * (dbdy * (Sqr(dWy) - dt) + dbdx * (dWx * dWy + Vxy))
         + dt * Sqr(slope_Dperp[0]) * (d2bdx2 * dWy - d2bdxdy * dWx) / 6.0;


// Add all the stuff
   dr_perp[0] = 0.5 * (slope_Dperp[0] + (slope_Dperp[1] + 3.0 * Sqr(gambar) * slope_Dperp[2]) / (1.0 + 3.0 * Sqr(gambar))) * dWx + Rx;
   dr_perp[1] = 0.5 * (slope_Dperp[0] + (slope_Dperp[1] + 3.0 * Sqr(gambar) * slope_Dperp[2]) / (1.0 + 3.0 * Sqr(gambar))) * dWy + Ry;
   // dr_perp[0] = sqrt(2.0 * Dperp) * dWx;
   // dr_perp[1] = sqrt(2.0 * Dperp) * dWy;
   dr_perp[2] = 0.0;

// Convert to global coordinates
   dr_perp.ChangeFromBasis(fa_basis);

// Reset common fields to actual current location (is this needed???)
   CommonFields();
   return false;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 03/12/2024
\param[out] position slopes
\param[out] momentum slopes
*/
template <typename Fields>
void TrajectoryGuidingDiff<Fields>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   TrajectoryGuidingBase::Slopes(slope_pos_istage, slope_mom_istage);
   TrajectoryGuidingDiff::DiffusionCoeff();
#if TRAJ_TIME_FLOW == TRAJ_TIME_FLOW_FORWARD
   slope_pos_istage += Vperp;
#else
   slope_pos_istage -= Vperp;
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 03/12/2024
*/
template <typename Fields>
void TrajectoryGuidingDiff<Fields>::PhysicalStep(void)
{
#if TRAJ_TIME_FLOW == TRAJ_TIME_FLOW_FORWARD
   dt_physical = cfl_adv_tg * _dmax / ((drift_vel + Vperp).Norm() + drift_safety_tg * _vel.Norm());
#else
   dt_physical = cfl_adv_tg * _dmax / ((drift_vel - Vperp).Norm() + drift_safety_tg * _vel.Norm());
#endif
   dt_physical = fmin(dt_physical, cfl_dif_gd * Sqr(_dmax) / Dperp);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/29/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Fields>
bool TrajectoryGuidingDiff<Fields>::Advance(void)
{
// Retrieve latest point of the trajectory and store locally
   Load();
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
#if STOCHASTIC_METHOD_PERP == 0
   EulerPerpDiffSlopes();
#elif STOCHASTIC_METHOD_PERP == 1
   MilsteinPerpDiffSlopes();
#elif STOCHASTIC_METHOD_PERP == 2
   if(RK2PerpDiffSlopes()) return true;
#else
   dr_perp = gv_zeros;
#endif

// If trajectory terminated (or is invalid) while computing slopes, exit advance function with true (step was taken)
   if(RKSlopes()) return true;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Advance trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If adaptive method error is unacceptable, exit advance function with false (step was not taken)
   if(RKStep()) return false;

// Stochastic displacement
   _pos += dr_perp;

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
