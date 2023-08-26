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
TrajectoryParker::TrajectoryParker(void)
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
TrajectoryParker::TrajectoryParker(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
                : TrajectoryBase(name_in, specie_in, status_in, presize_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
bool TrajectoryParker::IsSimmulationReady(void) const
{
   if(!TrajectoryBase::IsSimmulationReady()) return false;

// A diffusion object is required
   if(diffusion == nullptr) return false;
   else if(BITS_LOWERED(diffusion->GetStatus(), STATE_SETUP_COMPLETE)) return false;
   return true;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
void TrajectoryParker::SetStart(void)
{
// Call the base version of this function.
   TrajectoryBase::SetStart();

// Redefine mask
   _spdata._mask = BACKGROUND_U | BACKGROUND_B | BACKGROUND_gradU | BACKGROUND_gradB;
   spdata0._mask = BACKGROUND_U | BACKGROUND_B | BACKGROUND_gradU | BACKGROUND_gradB;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
void TrajectoryParker::FieldAlignedFrame(void)
{
   fa_basis[2] = _spdata.bhat;
   fa_basis[0] = GetSecondUnitVec(_spdata.bhat);
   fa_basis[1] = fa_basis[2] ^ fa_basis[0];
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
void TrajectoryParker::DiffusionCoeff(void)
try {
   GeoVector gradDperp, gradDpara;
   FieldAlignedFrame();

// Compute Dperp and grad(Dperp)
   Dperp = diffusion->GetComponent(0, _t, _pos, _mom, _spdata);
   gradDperp[0] = diffusion->GetDirectionalDerivative(0);
   gradDperp[1] = diffusion->GetDirectionalDerivative(1);
   gradDperp[2] = diffusion->GetDirectionalDerivative(2);

// Compute Dpara and grad(Dpara)
   Dpara = diffusion->GetComponent(1, _t, _pos, _mom, _spdata);
   gradDpara[0] = diffusion->GetDirectionalDerivative(0);
   gradDpara[1] = diffusion->GetDirectionalDerivative(1);
   gradDpara[2] = diffusion->GetDirectionalDerivative(2);

// Find components of divK in field aligned frame 
   divK[0] = fa_basis[0] * gradDperp;
   divK[1] = fa_basis[1] * gradDperp;
   divK[2] = fa_basis[2] * gradDpara;

// Convert "divK" to global coordinates
   divK.ChangeFromBasis(fa_basis);
}

catch(ExFieldError& exception) {
   // PrintError(__FILE__, __LINE__, "Error in field increment evaluation");
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
void TrajectoryParker::EulerDiffSlopes(void)
{
   double dWx, dWy, dWz;

// Generate stochastic factors
   dWx = sqrt(dt) * rng->GetNormal();
   dWy = sqrt(dt) * rng->GetNormal();
   dWz = sqrt(dt) * rng->GetNormal();

// Recompute Dperp and Dpara at the beginning of the step
   Dperp = diffusion->GetComponent(0, _t, _pos, _mom, _spdata);
   Dpara = diffusion->GetComponent(1, _t, _pos, _mom, _spdata);

// Compute random displacement
   dr_perp[0] = sqrt(2.0 * Dperp) * dWx;
   dr_perp[1] = sqrt(2.0 * Dperp) * dWy;
   dr_perp[2] = sqrt(2.0 * Dpara) * dWz;

// Convert to global coordinates
   dr_perp.ChangeFromBasis(fa_basis);
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\note "Bvec", "Bmag", and "bhat" must be available, so "CommonFields()" must be called prior to this function.
*/
void TrajectoryParker::DriftCoeff(void)
{
   // Compute grad(|B|) from grad(Bvec)
   gradBmag = _spdata.gradBvec * _spdata.bhat;
   // Compute curl(bhat/|B|) (times |B|^2)
   drift_vel = (_spdata.curlB() - 2.0 * (gradBmag ^ _spdata.bhat));
   // Scale by pvc/3q|B|^3 (1/|B|^2 from previous calculation)
   drift_vel *= _mom[0] * _vel[0] * c_code / (3.0 * charge[specie] * Cube(_spdata.Bmag));
   // Add bulk flow velocity
   drift_vel += _spdata.Uvec;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
*/
void TrajectoryParker::PhysicalStep(void)
{
// In a background with U = 0 and B = 0, then drift_vel = 0. For this reason a small fraction of the total speed is added to the characteristic speed.
   dt_physical = cfl_adv_tp * _spdata.dmax / (drift_vel.Norm() + drift_safety_tp * _vel[0]);
   dt_physical = fmin(dt_physical, cfl_dif_tp * Sqr(_spdata.dmax) / sqrt(Sqr(Dperp)+Sqr(Dpara)));
   dt_physical = fmin(dt_physical, cfl_dif_tp * _spdata.dmax / divK.Norm());
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\param[out] slope_pos_istage RK slope for position
\param[out] slope_mom_istage RK slope for momentum
*/
void TrajectoryParker::Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage)
{
   DriftCoeff();
   DiffusionCoeff();
   slope_pos_istage = drift_vel;
#if TIME_FLOW == 0
   slope_pos_istage -= divK;
#else
   slope_pos_istage += divK;
#endif
   slope_mom_istage[0] = -_mom[0] * _spdata.divU() / 3.0;
   slope_mom_istage[1] = 0.0;
   slope_mom_istage[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 06/07/2023
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
bool TrajectoryParker::Advance(void)
{
// Retrieve latest point of the trajectory and store locally
   Load();
   StoreLocal();

// The commomn fields and "dmax" have been computed at the end of Advance() or in SetStart() before the first step.
// Compute the slopes. The first two components for momentum are always zero for GC (the perpendicular momentum is determined from conservation of magnetic moment).
   Slopes(slope_pos[0], slope_mom[0]);

   PhysicalStep();
   dt = fmin(dt_physical, dt_adaptive);
   TimeBoundaryBefore();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute the RK slopes.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Stochastic RK slopes
// TODO: Add other (higher order) options for the stochastic step (e.g. Milstein or RK2)
#if STOCHASTIC_METHOD_DIFF == 0
   EulerDiffSlopes();
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

// If trajectory is not finished (in particular, spatial boundary not crossed), the fields can be computed
   if(BITS_LOWERED(_status, TRAJ_FINISH)) {
      CommonFields();
   };

// Add the new point to the trajectory.
   Store();

   return true;
};

};
