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
template <typename Fields>
TrajectoryParker<Fields>::TrajectoryParker(void)
                        : TrajectoryBase(traj_name_parker, 0, STATE_NONE, defsize_parker)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
\param[in] presize_in Whether to pre-allocate memory for trajectory arrays
*/
template <typename Fields>
TrajectoryParker<Fields>::TrajectoryParker(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                        : TrajectoryBase(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Fields>
bool TrajectoryParker<Fields>::IsSimmulationReady(void) const
{
   if (!TrajectoryBase::IsSimmulationReady()) return false;

// A diffusion object is required
   if (diffusion == nullptr) return false;
   else if (BITS_LOWERED(diffusion->GetStatus(), STATE_SETUP_COMPLETE)) return false;
   return true;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Fields>
void TrajectoryParker<Fields>::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Fields>
void TrajectoryParker<Fields>::FieldAlignedFrame(void)
{
   fa_basis[2] = _fields.HatMag();
   fa_basis[0] = GetSecondUnitVec(_fields.HatMag());
   fa_basis[1] = fa_basis[2] ^ fa_basis[0];
};

/*!
\author Juan G Alonso Guzman
\date 03/11/2024
*/
template <typename Fields>
void TrajectoryParker<Fields>::DiffusionCoeff(void)
try {
   int i,j;

#if TRAJ_PARKER_DIVK_METHOD == 0
   GeoVector pos_tmp;
   Fields fields_forw, fields_back;
   double Kperp_forw, Kperp_back, Kpara_forw, Kpara_back, Kappa_forw, Kappa_back;
   double delta = fmin(LarmorRadius(_mom[0], _fields.AbsMag(), specie), _dmax);

// Compute perpendicular and parallel diffusion coefficients and diffusion tensor.
   Kperp = diffusion->GetComponent(0, _t, _pos, _mom, _fields);
   Kpara = diffusion->GetComponent(1, _t, _pos, _mom, _fields);

// Loop over dimensions to find derivatives of Kappa.
   divK = gv_zeros;
// TODO: check if CommonFields returns a STATE_INVALID flag and revert to forward/backward (1st order) FD.
   for (j = 0; j < 3; j++) {
// Forward evaluation
      pos_tmp = _pos + delta * cart_unit_vec[j];
      CommonFields(_t, pos_tmp, _mom, fields_forw);
      Kperp_forw = diffusion->GetComponent(0, _t, pos_tmp, _mom, fields_forw);
      Kpara_forw = diffusion->GetComponent(1, _t, pos_tmp, _mom, fields_forw);
// Backward evaluation
      pos_tmp[j] -= 2.0 * delta;
      CommonFields(_t, pos_tmp, _mom, fields_back);
      Kperp_back = diffusion->GetComponent(0, _t, pos_tmp, _mom, fields_back);
      Kpara_back = diffusion->GetComponent(1, _t, pos_tmp, _mom, fields_back);
      for (i = 0; i < 3; i++) {
         Kappa_forw = Kperp_forw * (i == j ? 1.0 : 0.0) + (Kpara_forw - Kperp_forw) * fields_forw.HatMag()[j] * fields_forw.HatMag()[i];
         Kappa_back = Kperp_back * (i == j ? 1.0 : 0.0) + (Kpara_back - Kperp_back) * fields_back.HatMag()[j] * fields_back.HatMag()[i];
         divK[i] += 0.5 * (Kappa_forw - Kappa_back) / delta;
      };
   };
#else
   GeoVector gradKpara, gradKperp;
   GeoMatrix bhatbhat;

// Compute Kperp and grad(Kperp)
   Kperp = diffusion->GetComponent(0, _t, _pos, _mom, _fields);
   auto ddata = background->GetDerivativeData();
   gradKperp[0] = diffusion->GetDirectionalDerivative(0, ddata);
   gradKperp[1] = diffusion->GetDirectionalDerivative(1, ddata);
   gradKperp[2] = diffusion->GetDirectionalDerivative(2, ddata);

// Compute Kpara and grad(Kpara)
   Kpara = diffusion->GetComponent(1, _t, _pos, _mom, _fields);
   ddata = background->GetDerivativeData();
   gradKpara[0] = diffusion->GetDirectionalDerivative(0, ddata);
   gradKpara[1] = diffusion->GetDirectionalDerivative(1, ddata);
   gradKpara[2] = diffusion->GetDirectionalDerivative(2, ddata);

// Assemble diffusion tensor
   GeoVector bhat = _fields.HatMag();
   bhatbhat.Dyadic(bhat);
   divK = gradKperp + bhatbhat * (gradKpara - gradKperp)
        + (Kpara - Kperp) * (_fields.divbhat() * bhat + bhat * _fields.DelAbsMag());
#endif
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
template <typename Fields>
void TrajectoryParker<Fields>::EulerDiffSlopes(void)
{
   double dWx, dWy, dWz;
   FieldAlignedFrame();

// Generate stochastic factors
   dWx = sqrt(dt) * rng->GetNormal();
   dWy = sqrt(dt) * rng->GetNormal();
   dWz = sqrt(dt) * rng->GetNormal();

// Recompute Kperp and Kpara at the beginning of the step
   Kperp = diffusion->GetComponent(0, _t, _pos, _mom, _fields);
   Kpara = diffusion->GetComponent(1, _t, _pos, _mom, _fields);

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
template <typename Fields>
void TrajectoryParker<Fields>::DriftCoeff(void)
{
#ifdef TRAJ_PARKER_USE_B_DRIFTS
// Compute |B|*curl(b/|B|)
   drift_vel = (_fields.curlB() - 2.0 * (_fields.DelAbsMag() ^ _fields.HatMag())) / _fields.AbsMag();
// Scale by pvc/3q|B| = r_L*v/3
   drift_vel *= LarmorRadius(_mom[0], _fields.AbsMag(), specie) * _vel[0] / 3.0;
// Scale magnitude to an upper limit of v/2 if necessary.
   if (drift_vel.Norm() > 0.5 * _vel[0]) {
      drift_vel.Normalize();
      drift_vel *= 0.5 * _vel[0];
   };
// Add bulk flow velocity
   drift_vel += _fields.Vel();
#else
   drift_vel = _fields.Vel();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
template <typename Fields>
void TrajectoryParker<Fields>::PhysicalStep(void)
{
#if TRAJ_TIME_FLOW == TRAJ_TIME_FLOW_FORWARD
   dt_physical = cfl_adv_tp * _dmax / (drift_vel + divK).Norm();
#else
   dt_physical = cfl_adv_tp * _dmax / (drift_vel - divK).Norm();
#endif
   dt_physical = fmin(dt_physical, cfl_dif_tp * Sqr(_dmax) / fmax(Kperp, Kpara));
   dt_physical = fmin(dt_physical, cfl_acc_tp * 3.0 * dlnpmax / fabs(_fields.divU()));
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
template <typename Fields>
void TrajectoryParker<Fields>::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   DriftCoeff();
   DiffusionCoeff();

// Position slopes
   slope_pos_istage = drift_vel;
#if TRAJ_TIME_FLOW == TRAJ_TIME_FLOW_FORWARD
   slope_pos_istage += divK;
#else
   slope_pos_istage -= divK;
#endif

// Momentum slopes
   slope_mom_istage[0] = -_mom[0] * _fields.divU() / 3.0;
   slope_mom_istage[1] = 0.0;
   slope_mom_istage[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Fields>
bool TrajectoryParker<Fields>::Advance(void)
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
// TODO: Add other (higher order) options for the stochastic step (e.g. Milstein or RK2)
#if TRAJ_PARKER_STOCHASTIC_METHOD_DIFF == 0
   EulerDiffSlopes();
#else
   dr_perp = gv_zeros;
#endif

// If trajectory terminated (or is invalid) while computing slopes, exit advance function with true (step was taken)
   if (RKSlopes()) return true;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Advance trajectory
//----------------------------------------------------------------------------------------------------------------------------------------------------

// If adaptive method error is unacceptable, exit advance function with false (step was not taken)
   if (RKStep()) return false;

// Stochastic displacement
   _pos += dr_perp;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------
   HandleBoundaries();

// If trajectory is not finished (in particular, spatial boundary not crossed), the fields can be computed
   if (BITS_LOWERED(_status, TRAJ_FINISH)) CommonFields();

// Add the new point to the trajectory.
   Store();

   return true;
};

};
