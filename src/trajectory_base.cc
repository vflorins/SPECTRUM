/*!
\file trajectory_base.cc
\brief Implements a base class to calculate a trajectory
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <common/print_warn.hh>
#include <src/trajectory_base.hh>

namespace Spectrum {

#ifdef GEO_DEBUG
//! Upper limit on the number of steps in debug mode, use -1 for unlimited
constexpr unsigned int n_max_calls = -1;
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryBase protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/24/2020
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
TrajectoryBase<Trajectory, Fields, params>::TrajectoryBase(void)
              : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in    Readable name of the class
\param[in] status_in  Initial status
\param[in] specie_in  Particle's specie
\param[in] presize_in Initial lengths of the containers
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
TrajectoryBase<Trajectory, Fields, params>::TrajectoryBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in)
              : Params(name_in, specie_in, status_in)
{
   PreSize(presize_in);
};

/*!
\author Vladimir Florinski
\date 04/15/2022
\param[in] init_cap Initial array capacity
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::PreSize(int init_cap)
{
   traj_t.clear();
   traj_pos.clear();
   traj_mom.clear();
   presize = init_cap;
   traj_t.reserve(presize);
   traj_pos.reserve(presize);
   traj_mom.reserve(presize);
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::ResetAllBoundaries(void)
{
   unsigned int bnd;

   for (bnd = 0; bnd < bcond_t.size(); bnd++) {
      bcond_t[bnd]->SetScale(_dmax / c_code);
      bcond_t[bnd]->ResetBoundary(_t, _pos, _mom, _fields);
   };
   for (bnd = 0; bnd < bcond_s.size(); bnd++) {
      bcond_s[bnd]->SetScale(_dmax);
      bcond_s[bnd]->ResetBoundary(_t, _pos, _mom, _fields);
   };
   for (bnd = 0; bnd < bcond_m.size(); bnd++) {
      bcond_m[bnd]->SetScale(_mom.Norm());
      bcond_m[bnd]->ResetBoundary(_t, _pos, _mom, _fields);
   };
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::ComputeAllBoundaries(void)
{
   unsigned int bnd;

   for (bnd = 0; bnd < bcond_t.size(); bnd++) bcond_t[bnd]->ComputeBoundary(_t, _pos, _mom, _fields);
   for (bnd = 0; bnd < bcond_s.size(); bnd++) bcond_s[bnd]->ComputeBoundary(_t, _pos, _mom, _fields);
   for (bnd = 0; bnd < bcond_m.size(); bnd++) bcond_m[bnd]->ComputeBoundary(_t, _pos, _mom, _fields);
};

/*!
\author Vladimir Florinski
\date 02/06/2021
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::UpdateAllBoundaries(void)
{
   unsigned int bnd;

   for (bnd = 0; bnd < bcond_t.size(); bnd++) bcond_t[bnd]->RecordBoundary();
   for (bnd = 0; bnd < bcond_s.size(); bnd++) bcond_s[bnd]->RecordBoundary();
   for (bnd = 0; bnd < bcond_m.size(); bnd++) bcond_m[bnd]->RecordBoundary();
};

/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in]  t_in   Time point (use a negative value for trajectory end)
\param[out] pt     Nearest index on the small side
\param[out] weight Weight of the nearest point
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::GetIdx(double t_in, int& pt, double& weight) const
{
// For a negative "t_in", return an index of "-1" and weight of 0. The calling program must check for this to avoid a memory access error.
   if (t_in < 0.0) {
      pt = -1;
      weight = 0.0;
   }

// For a very large "t_in", return the next to last point, so the interpolator will use only the last point with a weight of 1. If the trajectory has only the starting point in it, the index returned will be "-1", and the calling program must check for this to avoid a memory access error.
   else if (t_in >= traj_t.back()) {
      pt = traj_t.size() - 2;
      weight = 0.0;
   }

// For additional safety, use the last argument to cap the output of "LocateInArray()".
   else {
      pt = LocateInArray(0, traj_t.size() - 1, traj_t.data(), t_in, true);
      weight = (traj_t[pt + 1] - t_in) / (traj_t[pt + 1] - traj_t[pt]);
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/07/2023
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::ReverseMomentum(void)
{
   _mom *= -1.0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024

This function should be called near the _beginning_ of the "Advance()" routine, after a call to "PhysicalStep()". Its only purpose is to adjust the time step to prevent an overshoot.
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::TimeBoundaryProximityCheck(void)
{
   unsigned int bnd;
   double delta, delta_next;
   if constexpr (params.timeflow == TimeFlow::forward) {
      delta_next = -sp_large * _dmax / c_code;
   }
   else {
      delta_next = sp_large * _dmax / c_code;
   }

// All boundaries have been evaluated at the end of the previous time step. For the first step this is done in "SetStart()".
   for (bnd = 0; bnd < bcond_t.size(); bnd++) {
      delta = bcond_t[bnd]->GetDelta();

// This gives the smallest delta in magnitude. Note that if two boundaries share a time stamp, only the first one will be processed. The first check is done to skip the event boundaries for which the crossing has already happened.
      if constexpr (params.timeflow == TimeFlow::forward) {
         if ((delta <= 0.0) && (delta > delta_next)) delta_next = delta;
      }
      else {
         if ((delta >= 0.0) && (delta < delta_next)) delta_next = delta;
      }
   };

// Check whether any boundaries _may_ be crossed and adjust the time step. For adaptive stepping the actual crossing may not occur until later.
   if constexpr (params.timeflow == TimeFlow::forward) {
      if (dt >= -delta_next) dt = fmax(-(1.0 + sp_little) * delta_next, sp_small * _dmax / c_code);
   }
   else {
      if (dt >=  delta_next) dt = fmax( (1.0 + sp_little) * delta_next, sp_small * _dmax / c_code);
   }
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\return True if a spatial boundary was crossed, or False otherwise

This function should be called each time the position is updated (e.g., inside the RK loop). The purpose is to catch the situation where the trajectory leaves the domain and the fields become unavailable making any further integration impossible. 
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
bool TrajectoryBase<Trajectory, Fields, params>::SpaceTerminateCheck(void)
try {
   uint16_t bnd_status;
   unsigned int bnd = 0, distro;

// Only check whether at least one absorbing boundary was crossed.
   bactive_s = -1;
   while ((bactive_s == -1) && (bnd < bcond_s.size())) {
      bcond_s[bnd]->ComputeBoundary(_t, _pos, _mom, _fields);
      bnd_status = bcond_s[bnd]->GetStatus();

// Terminal boundary crossed, so the trajectory cannot continue
      if (BITS_RAISED(bnd_status, BOUNDARY_CROSSED) && BITS_RAISED(bnd_status, BOUNDARY_TERMINAL)) {
         RAISE_BITS(_status, TRAJ_SPATIAL_CROSSED);
         RAISE_BITS(_status, TRAJ_FINISH);
         bactive_s = bnd;
      }
      else bnd++;
   };

// Handle distribution events
   if (bactive_s >= 0) {
      for (distro = 0; distro < distributions.size(); distro++) {
         action = bcond_s[bactive_s]->GetAction(distro);
         // todo modify signature of ProcessTrajectory
         if (action >= 0) distributions[distro]->ProcessTrajectory(traj_t[0], traj_pos[0], traj_mom[0], fields0, _t, _pos, _mom, _fields, magedata, action);
      };
   };

// If an exit spatial boundary was crossed, the fields may no longer be available, so the full RK step cannot be completed. In that case the function should return immediately and the last recorded position and momentum will be saved as if the step has completed. A check for momentum boundary is not needed; if one was crossed it will be recorded at the end of the step.
   if (BITS_RAISED(_status, TRAJ_SPATIAL_CROSSED) && BITS_RAISED(_status, TRAJ_FINISH)) {

#ifdef GEO_DEBUG
     PrintMessage(__FILE__, __LINE__, "Advance: The trajectory will terminate inside the RK loop", true);
#endif

      Store();
      return true;
   }
   else return false;
}

catch (ExUninitialized& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}

catch (ExBoundaryError& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/14/2025
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::CommonFields(void)
try {
// Compute fields and reset derivative data
   background->GetFields(_t, _pos, ConvertMomentum(), _fields);
// Set field-dependent dmax, for use by trajectories while advancing
   _dmax = background->GetDmax();
}

catch (ExUninitialized& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}

catch (ExCoordinates& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}

#if SERVER_TYPE != SERVER_SELF
catch (ExServerError& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}
#endif

catch (ExFieldError& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 02/21/2025
\param[in]  t_in   Time at which to compute fields
\param[in]  pos_in Position at which to compute fields
\param[in]  mom_in Momentum (p,mu,phi) coordinates
\param[out] spdata Spatial data at t_in and pos_in for output
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::CommonFields(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, Fields& fields)
try {
   background->GetFields(t_in, pos_in, mom_in, fields);
}

catch (ExUninitialized& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}

catch (ExCoordinates& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}

#if SERVER_TYPE != SERVER_SELF
catch (ExServerError& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}
#endif

catch (ExFieldError& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/08/2024
\return True if the domain was exited while computing the RK slopes, or False otherwise

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
bool TrajectoryBase<Trajectory, Fields, params>::RKSlopes(void)
{
   unsigned int istage, islope;

   for (istage = 1; istage < RK_Table.stages; istage++) {

// Advance to the current stage.
      for (islope = 0; islope < istage; islope++) {
         if constexpr (params.timeflow == TimeFlow::forward) {
            _t += RK_Table.a[istage] * dt;
            _pos += dt * RK_Table.b[istage][islope] * slope_pos[islope];
            _mom += dt * RK_Table.b[istage][islope] * slope_mom[islope];
         }
         else {
            _t -= RK_Table.a[istage] * dt;
            _pos -= dt * RK_Table.b[istage][islope] * slope_pos[islope];
            _mom -= dt * RK_Table.b[istage][islope] * slope_mom[islope];
         }
      };

// If an exit spatial boundary was crossed, the fields may no longer be available, so the full RK step cannot be completed. In that case the function should return immediately and the last recorded position and momentum will be saved as if the step has completed. A check for momentum boundary is not needed; if one was crossed it will be recorded at the end of the step.
      if (SpaceTerminateCheck()) return true;

// Obtain the fields at the new position. We can now compute p_perp and velocity even when using MM conservation.
      CommonFields();

// Compute/Recompute relevant momentum components based on transport
      MomentumCorrection();

// Find velocity and acceleration.
      _vel = Vel(_mom, specie);
      Slopes(slope_pos[istage], slope_mom[istage]);

// The slopes have been computed, so we can reset "_t", "_pos", and "_mom" to their values at the beginning of the step.
      LoadLocal();
   };

   return false;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024
\return True if adaptive error is unacceptable (> 1), or False otherwise

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
bool TrajectoryBase<Trajectory, Fields, params>::RKStep(void)
{
   unsigned int islope;
   double error = 1.0;
   GeoVector pos_lo;

   if constexpr (params.timeflow == TimeFlow::forward) {
      _t += dt;
   }
   else {
      _t -= dt;
   }
// For adaptive schemes "pos_lo" is computed with a lower order version (we only use position to test for accuracy).
   if constexpr (RK_Table.adaptive) pos_lo = _pos;
   for (islope = 0; islope < RK_Table.stages; islope++) {
      if constexpr (params.timeflow == TimeFlow::forward) {
         _pos += dt * RK_Table.v[islope] * slope_pos[islope];
         _mom += dt * RK_Table.v[islope] * slope_mom[islope];
         if (RK_Table.adaptive) pos_lo += dt * RK_Table.w[islope] * slope_pos[islope];
      }
      else {
         _pos -= dt * RK_Table.v[islope] * slope_pos[islope];
         _mom -= dt * RK_Table.v[islope] * slope_mom[islope];
         if (RK_Table.adaptive) pos_lo -= dt * RK_Table.w[islope] * slope_pos[islope];
      }
   };
   _vel = Vel(_mom, specie);

// Estimate the error in the adaptive RK method using position and compute the recommended time step.
   if constexpr (RK_Table.adaptive) {
      error = sqrt((_pos - pos_lo).Norm2() / Sqr(rk_tol_abs + rk_tol_rel * (_pos.Norm() + pos_lo.Norm())));
      dt_adaptive = dt * rk_adjust * pow(error, -1.0 / RK_Table.order);
      dt_adaptive = fmin(dt * rk_safety, dt_adaptive);
      dt_adaptive = fmax(dt / rk_safety, dt_adaptive);

// Don't make this step if the error is unacceptable. The FINISH flag must be cleared.
      if (error > 1.0 + sp_tiny) {
         LOWER_BITS(_status, TRAJ_FINISH);
         return true;
      };
   };

   return false;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 04/01/2024
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::HandleBoundaries(void)
{
   int bnd;
   unsigned int bnd_status, distro;

// Check _all_ boundary crossings. More than one may be crossed at any time.
   ComputeAllBoundaries();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle momentum boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------

   for (bnd = 0; bnd < bcond_m.size(); bnd++) {
      bnd_status = bcond_m[bnd]->GetStatus();
      if (BITS_RAISED(bnd_status, BOUNDARY_CROSSED)) {
         bactive_m = bnd;

// Mirroring event
         if (BITS_RAISED(bnd_status, BOUNDARY_REFLECT)) n_mirr++;

// Absorbing boundary. In the event of multiple momentum boundary crossing this will record the last one in the list (but this should not happen under normal operation).
         else if (BITS_RAISED(bnd_status, BOUNDARY_TERMINAL)) {
            RAISE_BITS(_status, TRAJ_MOMENTUM_CROSSED);
            RAISE_BITS(_status, TRAJ_FINISH);
         }

// TODO Event boundary can be added in the future
         else {};
      };
   };

// Handle distribution events
   if (bactive_m >= 0) {
      for (distro = 0; distro < distributions.size(); distro++) {
         action = bcond_m[bactive_m]->GetAction(distro);
         if (action >= 0) distributions[distro]->ProcessTrajectory(traj_t[0], traj_pos[0], traj_mom[0], fields0, _t, _pos, _mom, _fields, magedata, action);
      };
      bactive_m = -1;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle spatial boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------

   for (bnd = 0; bnd < bcond_s.size(); bnd++) {
      bnd_status = bcond_s[bnd]->GetStatus();
      if (BITS_RAISED(bnd_status, BOUNDARY_CROSSED)) {
         bactive_s = bnd;

// Reflection
         if (BITS_RAISED(bnd_status, BOUNDARY_REFLECT)) {

#ifdef GEO_DEBUG
            std::cerr << "Position and momentum before reflection: " << _pos << " " << _mom << std::endl;
#endif

// This is a very crude way to do a reflection. In the future one could improve on it by computing the precise boundary crossing location and reflecting the trajectory along the field line from there. However, this requires a lot of extra code.
            _pos -= 2.0 * bcond_s[bnd]->GetDelta() * bcond_s[bnd]->GetNormal();
            ReverseMomentum();
            n_refl++;

#ifdef GEO_DEBUG
            std::cerr << "Position and momentum after reflection: " << _pos << " " << _mom << std::endl;
#endif

// Recompute the boundary since the position and momentum have changed.
            bcond_s[bnd]->ComputeBoundary(_t, _pos, _mom, _fields);
// Manually decrement crossings_left, because recomputing the boundary changes _delta again so that _delta * _old_delta > 0.0
            bcond_s[bnd]->DecrCrossingsLeft();
// FIXME: This _should_ in principle count as mirroring, but doing so requires an evaluation of the mirroring boundary (TODO). Manually increment n_mirr for now.
            n_mirr++;
// NOTE: A solid reflective boundary does not exist in reality, but is an approximation for a region of stronger field outside of the simulation domain.
         }

// Absorbing boundary. Again, multiple boundaries are not handled properly at this time. (TODO)
         else if (BITS_RAISED(bnd_status, BOUNDARY_TERMINAL)) {
            RAISE_BITS(_status, TRAJ_SPATIAL_CROSSED);
            RAISE_BITS(_status, TRAJ_FINISH);
         }

// TODO Event boundary can be added in the future
         else {};
      };
   };

// Handle distribution events
   if (bactive_s >= 0) {
      for (distro = 0; distro < distributions.size(); distro++) {
         action = bcond_s[bactive_s]->GetAction(distro);
         if (action >= 0) distributions[distro]->ProcessTrajectory(traj_t[0], traj_pos[0], traj_mom[0], fields0,  _t, _pos, _mom, _fields, magedata, action);
      };
      bactive_s = -1;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handle temporal boundaries
//----------------------------------------------------------------------------------------------------------------------------------------------------

   for (bnd = 0; bnd < bcond_t.size(); bnd++) {
      bnd_status = bcond_t[bnd]->GetStatus();
      if (BITS_RAISED(bnd_status, BOUNDARY_CROSSED)) {
         bactive_t = bnd;
         if (BITS_RAISED(bnd_status, BOUNDARY_TERMINAL)) {
            RAISE_BITS(_status, TRAJ_TIME_CROSSED);
            RAISE_BITS(_status, TRAJ_FINISH);
         }

// TODO Event boundary can be added in the future
         else {};
      };
   };

// Handle distribution events
   if (bactive_t >= 0) {
      for (distro = 0; distro < distributions.size(); distro++) {
         action = bcond_t[bactive_t]->GetAction(distro);
         if (action >= 0) distributions[distro]->ProcessTrajectory(traj_t[0], traj_pos[0], traj_mom[0], fields0, _t, _pos, _mom, _fields, magedata, action);
      };
      bactive_t = -1;
   };

// Update all boundary objects, recording the crossings.
   UpdateAllBoundaries();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 04/19/2022
\return True if a step was taken

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
bool TrajectoryBase<Trajectory, Fields, params>::RKAdvance(void)
{
// Retrieve latest point of the trajectory and store locally
   Load();
   StoreLocal();

// The common fields and "dmax" have been computed at the end of Advance() or in SetStart() before the first step.
// Compute the slopes. The first two components for momentum are always zero for GC (the perpendicular momentum is determined from conservation of magnetic moment).
   Slopes(slope_pos[0], slope_mom[0]);

// Figure out the physical time step and take into account the adaptive step recommendation. Check if the time step is too to step over the end of the run and adjust if necessary.
   PhysicalStep();
   dt = fmin(dt_physical, dt_adaptive);
   TimeBoundaryProximityCheck();

// Compute the RK slopes. If a trajectory terminated (or is invalid) while computing slopes, exit the function.
   if (RKSlopes()) return true;

// Advance the trajectory. If the adaptive method error is unacceptable, exit the function.
   if (RKStep()) return false;

// Handle boundaries
   HandleBoundaries();

// If trajectory is not finished (in particular, spatial boundary not crossed), the fields can be computed and momentum corrected
   if (BITS_LOWERED(_status, TRAJ_FINISH)) {
      CommonFields();
      MomentumCorrection();
   };

// Add the new point to the trajectory.
   Store();

   return true;
};

/*!
\author Vladimir Florinski
\date 09/30/2022
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::MomentumCorrection(void)
{
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
bool TrajectoryBase<Trajectory, Fields, params>::IsSimmulationReady(void) const
{
// Particle specie must be known
   if ((specie < 0) || (specie >= MAX_PARTICLE_SPECIES)) return false;

// A background object is required
   if (!background) return false;
   else if (BITS_LOWERED(background->GetStatus(), STATE_SETUP_COMPLETE)) return false;

// Time initial condition is required
   if (!icond_t) return false;
   else if (BITS_LOWERED(icond_t->GetStatus(), STATE_SETUP_COMPLETE)) return false;

// Space initial condition is required
   if (!icond_s) return false;
   else if (BITS_LOWERED(icond_s->GetStatus(), STATE_SETUP_COMPLETE)) return false;

// Momentum initial condition is required
   if (!icond_m) return false;
   else if (BITS_LOWERED(icond_m->GetStatus(), STATE_SETUP_COMPLETE)) return false;

// A minimum of one time boundary is required
   if (bcond_t.size() == 0) return false;
   else if (BITS_LOWERED(bcond_t[0]->GetStatus(), STATE_SETUP_COMPLETE)) return false;

// A distribution need not be present if we are simulating a single trajectory
   return true;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryBase public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] specie_in Index of the particle species defined in physics.hh
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::SetSpecie(unsigned int specie_in)
{
   Params::SetSpecie(specie_in);
// The factor multiplying "SpeciesCharges[]" is applied in order to marry particle and fluid scales. See "LarmorRadius()" and "CyclotronFrequency()" functions in physics.hh for reference.
   q = charge_mass_particle * SpeciesCharges[specie];
   if (background != nullptr) background->SetSpecie(specie);
   if (diffusion != nullptr) diffusion->SetSpecie(specie);

   for (auto& bnd : bcond_t) bnd->SetSpecie(specie);
   for (auto& bnd : bcond_s) bnd->SetSpecie(specie);
   for (auto& bnd : bcond_m) bnd->SetSpecie(specie);

   if (icond_t != nullptr) icond_t->SetSpecie(specie);
   if (icond_s != nullptr) icond_s->SetSpecie(specie);
   if (icond_m != nullptr) icond_m->SetSpecie(specie);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] background_in Background object for type recognition
\param[in] container_in  Data container for initializating the background object
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::AddBackground(const BackgroundBase& background_in, const DataContainer& container_in)
{
   background = background_in.Clone();
   background->SetSpecie(specie);
   background->ConnectRNG(rng);
   background->SetupObject(container_in);
   
   if (IsSimmulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] diffusion_in Diffusion object for type recognitions
\param[in] container_in Data container for initializating the diffusion object
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::AddDiffusion(const DiffusionBase& diffusion_in, const DataContainer& container_in)
{
   diffusion = diffusion_in.Clone();
   diffusion->SetSpecie(specie);
   diffusion->SetupObject(container_in);

   if (IsSimmulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] boundary_in  Boundary object for type recognition
\param[in] container_in Data container for initializating the boundary object
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in)
{
// Time boundary
   if (BITS_RAISED(boundary_in.GetStatus(), BOUNDARY_TIME)) {
      bcond_t.push_back(boundary_in.Clone());
      bcond_t.back()->SetSpecie(specie);
      bcond_t.back()->SetupObject(container_in);
   }

// Spatial boundary
   else if (BITS_RAISED(boundary_in.GetStatus(), BOUNDARY_SPACE)) {
      bcond_s.push_back(boundary_in.Clone());
      bcond_s.back()->SetSpecie(specie);
      bcond_s.back()->SetupObject(container_in);
   }

// Momentum boundary
   else if (BITS_RAISED(boundary_in.GetStatus(), BOUNDARY_MOMENTUM)) {
      bcond_m.push_back(boundary_in.Clone());
      bcond_m.back()->SetSpecie(specie);
      bcond_m.back()->SetupObject(container_in);
   }

   if (IsSimmulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] initial_in   Initial object for type recognition
\param[in] container_in Data container for initializating the initial object
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::AddInitial(const InitialBase& initial_in, const DataContainer& container_in)
{
// Time condition
   if (BITS_RAISED(initial_in.GetStatus(), INITIAL_TIME)) {
      icond_t = initial_in.Clone();
      icond_t->SetSpecie(specie);
      icond_t->ConnectRNG(rng);
      icond_t->SetupObject(container_in);
   }

// Spatial condition
   else if (BITS_RAISED(initial_in.GetStatus(), INITIAL_SPACE)) {
      icond_s = initial_in.Clone();
      icond_s->SetSpecie(specie);
      icond_s->ConnectRNG(rng);
      icond_s->SetupObject(container_in);
   }

// Momentum condition
   else if (BITS_RAISED(initial_in.GetStatus(), INITIAL_MOMENTUM)) {
      icond_m = initial_in.Clone();
      icond_m->SetSpecie(specie);
      icond_m->ConnectRNG(rng);
      icond_m->SetupObject(container_in);
   };

   if (IsSimmulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 06/22/2023
\return Minimum |B| along trajectory
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
double TrajectoryBase<Trajectory, Fields, params>::GetBmagMin(void) const
{
   if constexpr (params.record_bmag_extrema)
      return magedata.Bmag_min;
};

/*!
\author Juan G Alonso Guzman
\date 06/22/2023
\return Maximum |B| along trajectory
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
double TrajectoryBase<Trajectory, Fields, params>::GetBmagMax(void) const
{
   if constexpr (params.record_bmag_extrema)
      return magedata.Bmag_max;
};

/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Position
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
GeoVector TrajectoryBase<Trajectory, Fields, params>::GetPosition(double t_in) const
{
   int pt;
   double weight;

   if constexpr (params.record_trajectory) {
      GetIdx(t_in, pt, weight);
      if (pt < 0) return traj_pos[0];
      else return weight * traj_pos[pt] + (1.0 - weight) * traj_pos[pt + 1];
   }
   else {
      std::cerr << "Cannot get position with respect to time because it is not being recorded." << std::endl;
      return gv_zeros;
   }
};

/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Velocity
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
GeoVector TrajectoryBase<Trajectory, Fields, params>::GetVelocity(double t_in) const
{
   int pt;
   double weight, mom1, vel1, mom2, vel2;

   if constexpr (params.record_trajectory) {
      GetIdx(t_in, pt, weight);
      if (pt < 0) {
         mom1 = traj_mom[0].Norm();
         vel1 = Vel(mom1, specie);
         return (vel1 / mom1) * traj_mom[0];
      }

// Linear interpolation between nearest frames
      else {
         mom1 = traj_mom[pt].Norm();
         vel1 = Vel(mom1, specie);
         mom2 = traj_mom[pt + 1].Norm();
         vel2 = Vel(mom2, specie);
         return weight * (vel1 / mom1) * traj_mom[pt] + (1.0 - weight) * (vel2 / mom2) * traj_mom[pt + 1];
      };
   }
   else {
      std::cerr << "Cannot get velocity with respect to time because it is not being recorded." << std::endl;
      return gv_zeros;
   }
};

/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Kinetic energy
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
double TrajectoryBase<Trajectory, Fields, params>::GetEnergy(double t_in) const
{
   int pt;
   double weight;

   if constexpr (params.record_trajectory) {
      GetIdx(t_in, pt, weight);
      if (pt < 0) return EnrKin(traj_mom[0].Norm(), specie);
      else return weight * EnrKin(traj_mom[pt].Norm(), specie) + (1.0 - weight) * EnrKin(traj_mom[pt + 1].Norm(), specie);
   }
   else {
      std::cerr << "Cannot get energy with respect to time because it is not being recorded." << std::endl;
      return 0.0;
   }
};

/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Integral along trajectory
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
double TrajectoryBase<Trajectory, Fields, params>::GetDistance(double t_in) const
{
   int pt, ipt;
   double weight, length = 0.0;
   GeoVector pos_final;

   if constexpr (params.record_trajectory) {
      GetIdx(t_in, pt, weight);
      if (pt >= 0) {
         for (ipt = 0; ipt < pt; ipt++) length += (traj_pos[ipt + 1] - traj_pos[ipt]).Norm();
         pos_final = weight * traj_pos[pt] + (1.0 - weight) * traj_pos[pt + 1];
         length += (pos_final - traj_pos[pt]).Norm();
      };
      return length;
   }
   else {
      std::cerr << "Cannot get distance along trajectory because it is not being recorded." << std::endl;
      return 0.0;
   }
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024

To start a new trajectory its objects must be set to their initial state. This function determines the initial position and momentum from the respective distributions, calculates the fields, initializes the boundaries at the initial poasition, and resets the counters. A time step evaluation is not performed because it is done in"Advance()" at the beginning of each step.
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::SetStart(void)
try {

// Get the starting time from the initial time distribution
   _t = icond_t->GetTimeSample();
// Get the starting position from the initial space distribution.
   _pos = icond_s->GetPosSample();
// Get a momentum sample along an arbitrary axis (bhat is unknown at this step). Only the momentum magnitude is needed for the first call to CommonFields().
   _mom = icond_m->GetMomSample(gv_ones);

// Obtain the fields for that position (this initializes _dmax)
   CommonFields();

// Record the initial spatial data for distribution purposes.
   fields0 = _fields;
   if constexpr (params.record_bmag_extrema) {
      magedata.Bmag_min = _fields.AbsMag();
      magedata.Bmag_max = _fields.AbsMag();
      magedata.Bmag_min_initial = _fields.AbsMag();
      magedata.Bmag_max_initial = _fields.AbsMag();
   }

// Get the starting momentum from the distribution along the correct axis (bhat is now determined).
   if constexpr (Fields::HatMag_found()) {
      _mom = icond_m->GetMomSample(_fields.HatMag());
   }
   else {
// TODO
      ;
   }
   _vel = Vel(_mom, specie);

// Adaptive step must be large at first so that "dt" starts with a physical step.
   dt_adaptive = sp_large * _dmax / c_code;

// Re-initialize the trajectory arrays
   if constexpr (params.record_trajectory) {
      PreSize(presize);
   }
   else {
      PreSize(1);
   }

// The first element of traj_* arrays is necessary even if trajectories are not being recorded because it is used by the distributions in "ProcessTrajectory"
   traj_t.push_back(_t);
   traj_pos.push_back(_pos);
   traj_mom.push_back(_mom);
   
// Lower all flags
   LOWER_BITS(_status, STATE_INVALID);
   LOWER_BITS(_status, TRAJ_FINISH);
   LOWER_BITS(_status, TRAJ_TIME_CROSSED);
   LOWER_BITS(_status, TRAJ_SPATIAL_CROSSED);
   LOWER_BITS(_status, TRAJ_MOMENTUM_CROSSED);
   LOWER_BITS(_status, TRAJ_DISCARD);

// Reset reflection counters
   n_refl = 0;
   n_mirr = 0;

// Reset all boundary objects. From now on these objects will track and record boundary crossings automatically.
   bactive_t = bactive_s = bactive_m = -1;
   ResetAllBoundaries();
}

catch (ExUninitialized& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
}

catch (ExFieldError& exception) {
   RAISE_BITS(_status, TRAJ_DISCARD);
   throw;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 12/17/2020
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::Integrate(void)
{
   bool was_advanced;
   if constexpr (params.traj_adv_safety_level == TrajectorySafetyLevel::high) {
      traj_safety.time_step_adaptations = 0;
   }

// Time loop is very simple - a single call to "Advance()" followed by a global boundary update. It is the responsibility of Advance() to record the distribution on boundary crossing events.
   while (BITS_LOWERED(_status, TRAJ_FINISH) && BITS_LOWERED(_status, TRAJ_DISCARD)) {

// Attempt to advance trajectory by one segment
      was_advanced = Advance();

// Update |B| extrema
      if constexpr (params.record_bmag_extrema)
         if (was_advanced) {
            magedata.Bmag_min = fmin(magedata.Bmag_min, _fields.AbsMag());
            magedata.Bmag_max = fmax(magedata.Bmag_max, _fields.AbsMag());
         }

      if constexpr (params.traj_adv_safety_level > 1) {
// Too many steps were taken - terminate
         if (Segments() > traj_safety.max_trajectory_steps) {
            RAISE_BITS(_status, TRAJ_DISCARD);
            throw ExMaxStepsReached();
         };

// Too many time adaptations were performed - terminate
         if (was_advanced) traj_safety.time_step_adaptations = 0;
         else {
            traj_safety.time_step_adaptations++;
            if (traj_safety.time_step_adaptations > traj_safety.max_time_adaptations) {
               RAISE_BITS(_status, TRAJ_DISCARD);
               throw ExMaxTimeAdaptsReached();
            };
         };
      };

      if constexpr (params.traj_adv_safety_level > 0) {
// Time step is too small - terminate
         if (dt < sp_tiny * _dmax / c_code) {
            RAISE_BITS(_status, TRAJ_DISCARD);
            throw ExTimeStepTooSmall();
         }
         else if (!std::isnormal(dt)) {
            RAISE_BITS(_status, TRAJ_DISCARD);
            throw ExTimeStepNan();
         };
      };
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 02/17/2023
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::StopBackground(void)
{
   background->StopServerFront();
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 04/22/2022
\param[in] output Which boundary array to use (time, space, or momentum)
\param[in] bnd    Which boundary condition to use
\return int number of crossings
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
int TrajectoryBase<Trajectory, Fields, params>::Crossings(unsigned int output, unsigned int bnd) const
{
   if (bnd < 0) return 0;

   if ((output == 0) && (bnd < bcond_t.size())) return bcond_t[bnd]->CrossingsMade();
   else if ((output == 1) && (bnd < bcond_s.size())) return bcond_s[bnd]->CrossingsMade();
   else if ((output == 2) && (bnd < bcond_m.size())) return bcond_m[bnd]->CrossingsMade();
   return 0;
};

/*!
\author Vladimir Florinski
\date 07/13/2020
\param[in] traj_name  File name
\param[in] phys_units Use physical units for output
\param[in] output     Which coordinates to print
\param[in] stride     Distance between points in the output (optional). If stride = 0, output based on dt_out.
\param[in] dt_out     Time increment at which to output quantities when stride = 0
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::PrintTrajectory(const std::string traj_name, bool phys_units, unsigned int output,
                                     unsigned int stride, double dt_out) const
{
   unsigned int pt, iter_out = 0, max_out = 1000000;
   double mom_mag, vm_ratio, t_out = 0.0, engkin_t;
   GeoVector pos_t, vel_t;
   std::ofstream trajfile;

   if constexpr (params.record_trajectory) {
      trajfile.open(traj_name.c_str());

// Generate multiple column output
      trajfile << std::setprecision(12);

      if (stride) {
         for (pt = 0; pt < traj_t.size(); pt += stride) {
//FIXME: This computation of momentum magnitude is not guaranteed to work for focused transport. It is only approximately correct when magnitude (_mom[0]) >> pitch angle cosine (_mom[1]).
            mom_mag = traj_mom[pt].Norm();
            vm_ratio = Vel(mom_mag, specie) / mom_mag;

            if (output & 0x01) trajfile << std::setw(20) << traj_t[pt] * (phys_units ? unit_time_fluid : 1.0);
            if (output & 0x02) trajfile << std::setw(20) << traj_pos[pt][0] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x04) trajfile << std::setw(20) << traj_pos[pt][1] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x08) trajfile << std::setw(20) << traj_pos[pt][2] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x10) trajfile << std::setw(20) << vm_ratio * traj_mom[pt][0] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x20) trajfile << std::setw(20) << vm_ratio * traj_mom[pt][1] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x40) trajfile << std::setw(20) << vm_ratio * traj_mom[pt][2] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x80) trajfile << std::setw(20) << EnrKin(mom_mag, specie) * (phys_units ? unit_energy_particle : 1.0);
            trajfile << std::endl;
         };
      }
      else {
         while (t_out < traj_t.back() && iter_out < max_out) {
            pos_t = GetPosition(t_out);
            vel_t = GetVelocity(t_out);
            engkin_t = GetEnergy(t_out);

            if (output & 0x01) trajfile << std::setw(20) << t_out * (phys_units ? unit_time_fluid : 1.0);
            if (output & 0x02) trajfile << std::setw(20) << pos_t[0] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x04) trajfile << std::setw(20) << pos_t[1] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x08) trajfile << std::setw(20) << pos_t[2] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x10) trajfile << std::setw(20) << vel_t[0] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x20) trajfile << std::setw(20) << vel_t[1] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x40) trajfile << std::setw(20) << vel_t[2] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x80) trajfile << std::setw(20) << engkin_t * (phys_units ? unit_energy_particle : 1.0);
            trajfile << std::endl;

            t_out += dt_out;
         };
      };

      trajfile.close();
   }
   else {
      std::cerr << "Cannot print trajectory because it is not being recorded." << std::endl;
      return;
   }
};

/*!
\author Vladimir Florinski
\date 10/15/2020
\param[in] traj_name  File name
\param[in] phys_units Use physical units for output
\param[in] stride     Distance between points in the output (optional)
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride) const
{
   unsigned int pt;
   std::ofstream trajfile;

   if constexpr (params.record_trajectory) {
      trajfile.open(traj_name.c_str());

// Generate CSV output
      trajfile << std::setprecision(12);
      for (pt = 0; pt < traj_t.size(); pt += stride) {
         trajfile << std::setw(20) << traj_pos[pt][0] * (phys_units ? unit_length_fluid : 1.0);
         trajfile << ",";
         trajfile << std::setw(20) << traj_pos[pt][1] * (phys_units ? unit_length_fluid : 1.0);
         trajfile << ",";
         trajfile << std::setw(20) << traj_pos[pt][2] * (phys_units ? unit_length_fluid : 1.0);
         trajfile << std::endl;
      };

      trajfile.close();
   }
   else {
      std::cerr << "Cannot print trajectory because it is not being recorded." << std::endl;
      return;
   }
};

/*!
\author Vladimir Florinski
\date 12/27/2021
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::InterpretStatus(void) const
{
   std::cerr << "Trajectory status: ";
   if (BITS_RAISED(_status, TRAJ_DISCARD)) std::cerr << "discarded\n";

// These three states correspond to an absorbing boundary leading to a termination. The status is preserved and can be checked after a trajectory completion.
   else if (BITS_RAISED(_status, TRAJ_FINISH)) {
      if (BITS_RAISED(_status, TRAJ_TIME_CROSSED)) std::cerr << "time expired\n";
      else if (BITS_RAISED(_status, TRAJ_MOMENTUM_CROSSED)) std::cerr << "momentum boundary crossed\n";
      else if (BITS_RAISED(_status, TRAJ_SPATIAL_CROSSED)) std::cerr << "spatial boundary crossed\n";
   }
   else std::cerr << "in progress\n";
};

/*!
\author Vladimir Florinski
\date 02/22/2023
*/
template <typename Trajectory, typename Fields, TrajectoryParams params>
void TrajectoryBase<Trajectory, Fields, params>::PrintInfo(void) const
{
   int obj;
   std::cerr << std::endl;
   std::cerr << "Printing trajectory object information\n";
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Type\n";
   std::cerr << "   " << class_name << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Distributions\n";
   for (obj = 0; obj < distributions.size(); obj++) {
      std::cerr << "   " << distributions[obj]->GetName() << std::endl;
   };
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Background\n";
   if (background) std::cerr << "   " << background->GetName() << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Diffusion\n";
   if (diffusion) std::cerr << "   " << diffusion->GetName() << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Boundaries\n";
   for (obj = 0; obj < bcond_t.size(); obj++) {
      std::cerr << "   " << bcond_t[obj]->GetName() << std::endl;
   };
   for (obj = 0; obj < bcond_s.size(); obj++) {
      std::cerr << "   " << bcond_s[obj]->GetName() << std::endl;
   };
   for (obj = 0; obj < bcond_m.size(); obj++) {
      std::cerr << "   " << bcond_m[obj]->GetName() << std::endl;
   };
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Initials\n";
   if (icond_s) std::cerr << "   " << icond_s->GetName() << std::endl;
   if (icond_m) std::cerr << "   " << icond_m->GetName() << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
};

};
