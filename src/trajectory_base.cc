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

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryBase protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/24/2020
*/
template <typename Background, typename Diffusion>
TrajectoryBase<Background, Diffusion>::TrajectoryBase(void)
              : Params("", STATE_NONE),
              records(HConfig::record_trajectory_segment_presize)
{
};

/*!
\author Vladimir Florinski
\date 01/28/2022
\param[in] name_in    Readable name of the class
\param[in] status_in  Initial status
*/
template <typename Background, typename Diffusion>
TrajectoryBase<Background, Diffusion>::TrajectoryBase(const std::string& name_in, uint16_t status_in)
              : Params(name_in, status_in),
                records(HConfig::record_trajectory_segment_presize)
{
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::ResetAllBoundaries(void)
{
   unsigned int bnd;

   for (bnd = 0; bnd < bcond_t.size(); bnd++) {
      bcond_t[bnd]->SetScale(_dmax / c_code);
      bcond_t[bnd]->ResetBoundary(_coords, _fields);
   };
   for (bnd = 0; bnd < bcond_s.size(); bnd++) {
      bcond_s[bnd]->SetScale(_dmax);
      bcond_s[bnd]->ResetBoundary(_coords, _fields);
   };
   for (bnd = 0; bnd < bcond_m.size(); bnd++) {
      bcond_m[bnd]->SetScale(_coords.Mom().Norm());
      bcond_m[bnd]->ResetBoundary(_coords, _fields);
   };
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::ComputeAllBoundaries(void)
{
   unsigned int bnd;

   for (bnd = 0; bnd < bcond_t.size(); bnd++) bcond_t[bnd]->ComputeBoundary(_coords, _fields);
   for (bnd = 0; bnd < bcond_s.size(); bnd++) bcond_s[bnd]->ComputeBoundary(_coords, _fields);
   for (bnd = 0; bnd < bcond_m.size(); bnd++) bcond_m[bnd]->ComputeBoundary(_coords, _fields);
};

/*!
\author Vladimir Florinski
\date 02/06/2021
*/
template <typename Background>
void TrajectoryBase<Background>::UpdateAllBoundaries(void)
{
   unsigned int bnd;

   for (bnd = 0; bnd < bcond_t.size(); bnd++) bcond_t[bnd]->RecordBoundary();
   for (bnd = 0; bnd < bcond_s.size(); bnd++) bcond_s[bnd]->RecordBoundary();
   for (bnd = 0; bnd < bcond_m.size(); bnd++) bcond_m[bnd]->RecordBoundary();
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/11/2025
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::ReverseMomentum(void)
{
   _coords.Mom() *= -1.0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024

This function should be called near the _beginning_ of the "Advance()" routine, after a call to "PhysicalStep()". Its only purpose is to adjust the time step to prevent an overshoot.
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::TimeBoundaryProximityCheck(void)
{
   unsigned int bnd;
   double delta, delta_next;
   if constexpr (HConfig::time_flow == TimeFlow::forward) {
      delta_next = -sp_large * _dmax / c_code;
   }
   else {
      delta_next = sp_large * _dmax / c_code;
   }

// All boundaries have been evaluated at the end of the previous time step. For the first step this is done in "SetStart()".
   for (bnd = 0; bnd < bcond_t.size(); bnd++) {
      delta = bcond_t[bnd]->GetDelta();

// This gives the smallest delta in magnitude. Note that if two boundaries share a time stamp, only the first one will be processed. The first check is done to skip the event boundaries for which the crossing has already happened.
      if constexpr (HConfig::time_flow == TimeFlow::forward) {
         if ((delta <= 0.0) && (delta > delta_next)) delta_next = delta;
      }
      else {
         if ((delta >= 0.0) && (delta < delta_next)) delta_next = delta;
      }
   };

// Check whether any boundaries _may_ be crossed and adjust the time step. For adaptive stepping the actual crossing may not occur until later.
   if constexpr (HConfig::time_flow == TimeFlow::forward) {
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
template <typename Background, typename Diffusion>
bool TrajectoryBase<Background, Diffusion>::SpaceTerminateCheck(void)
try {
   uint16_t bnd_status;
   int bnd = 0, distro;

// Only check whether at least one absorbing boundary was crossed.
   bactive_s = -1;
   while ((bactive_s == -1) && (bnd < bcond_s.size())) {
      bcond_s[bnd]->ComputeBoundary(_coords, _fields);
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
         if (action >= 0) distributions[distro]->ProcessTrajectory(coords0, fields0, records.GetMagInitial(), _coords, _fields, records.GetMagExtrema(), action);
      };
   };

// If an exit spatial boundary was crossed, the fields may no longer be available, so the full RK step cannot be completed. In that case the function should return immediately and the last recorded position and momentum will be saved as if the step has completed. A check for momentum boundary is not needed; if one was crossed it will be recorded at the end of the step.
   if (BITS_RAISED(_status, TRAJ_SPATIAL_CROSSED) && BITS_RAISED(_status, TRAJ_FINISH)) {

      if constexpr (HConfig::build_mode == BuildMode::debug) {
         PrintMessage(__FILE__, __LINE__, "Advance: The trajectory will terminate inside the RK loop", true);
      }

      records.Store(_coords);
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
\brief Compute fields for a non-canonical coordinate
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 02/21/2025
\param[in]  coords   Coordinates (typically Time, Position, Momentum) at which to compute fields
\param[out] fields Fields requested to be populated by the background
*/
template <typename Background, typename Diffusion>
template <typename Coordinates, typename Fields, typename RequestedFields>
void TrajectoryBase<Background, Diffusion>::CommonFields(Coordinates& coords, Fields& fields)
try {
// Compute fields and reset derivative data
   background->template GetFields<Coordinates, Fields, RequestedFields>(coords, fields);
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
\date 10/08/2024
\return True if the domain was exited while computing the RK slopes, or False otherwise

If the state at return contains the TRAJ_TERMINATE flag, the calling program must stop this trajectory. If the state at the end contains the TRAJ_DISCARD flag, the calling program must reject this trajectory (and possibly repeat the trial with a different random number).
*/
template <typename Background, typename Diffusion>
bool TrajectoryBase<Background, Diffusion>::RKSlopes(void)
{
// When the function exits we are finished with the input coords, and can free the memory.
   unsigned int istage, islope;

   for (istage = 1; istage < ButcherTable::data.rk_stages; istage++) {

// Advance to the current stage.
      for (islope = 0; islope < istage; islope++) {
         if constexpr (HConfig::time_flow == TimeFlow::forward) {
            _coords.Time() += butcher_table.a[istage] * dt;
            _coords.Pos() += dt * butcher_table.b[istage][islope] * slope_pos[islope];
            _coords.Mom() += dt * butcher_table.b[istage][islope] * slope_mom[islope];
         }
         else {
            _coords.Time() -= butcher_table.a[istage] * dt;
            _coords.Pos() -= dt * butcher_table.b[istage][islope] * slope_pos[islope];
            _coords.Mom() -= dt * butcher_table.b[istage][islope] * slope_mom[islope];
         }
      };

// If an exit spatial boundary was crossed, the fields may no longer be available, so the full RK step cannot be completed. In that case the function should return immediately and the last recorded position and momentum will be saved as if the step has completed. A check for momentum boundary is not needed; if one was crossed it will be recorded at the end of the step.
      if (SpaceTerminateCheck()) return true;

// Obtain the fields at the new position. We can now compute p_perp and velocity even when using MM conservation.
      CommonFields(_coords, _fields);

// Compute/Recompute relevant momentum components based on transport
      MomentumCorrection();

// Find velocity and acceleration.
      if constexpr (TrajectoryCoordinates::Vel_found())
         _coords.Vel() = Vel<specie>(_coords.Mom());
      Slopes(slope_pos[istage], slope_mom[istage]);

// The slopes have been computed, so we can reset _coords to their values at the beginning of the step.
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
template <typename Background, typename Diffusion>
bool TrajectoryBase<Background, Diffusion>::RKStep(void)
{
   using BT = ButcherTable;
   unsigned int islope;
   double error = 1.0;
   GeoVector pos_lo;

   if constexpr (HConfig::time_flow == TimeFlow::forward) {
      _coords.Time() += dt;
   }
   else {
      _coords.Time() -= dt;
   }
// For adaptive schemes "pos_lo" is computed with a lower order version (we only use position to test for accuracy).
   if constexpr (BT::data.adaptive)
      pos_lo = _coords.Pos();
   for (islope = 0; islope < BT::data.rk_stages; islope++) {

      if constexpr (HConfig::time_flow == TimeFlow::forward) {
         _coords.Pos() += dt * butcher_table.v[islope] * slope_pos[islope];
         _coords.Mom() += dt * butcher_table.v[islope] * slope_mom[islope];
         if constexpr (BT::data.adaptive)
            pos_lo += dt * butcher_table.w[islope] * slope_pos[islope];
      }
      else {
         _coords.Pos() -= dt * butcher_table.v[islope] * slope_pos[islope];
         _coords.Mom() -= dt * butcher_table.v[islope] * slope_mom[islope];
         if constexpr (BT::data.adaptive)
            pos_lo -= dt * butcher_table.w[islope] * slope_pos[islope];
      }

   };
// Update velocity fields
   _coords.Vel() = Vel<specie>(_coords.Mom());


// Estimate the error in the adaptive RK method using position and compute the recommended time step.
   if (BT::data.adaptive) {
      error = sqrt((_coords.Pos() - pos_lo).Norm2() / Sqr(BT::rk_tol_abs + BT::rk_tol_rel * (_coords.Pos().Norm() + pos_lo.Norm())));
      dt_adaptive = dt * BT::rk_adjust * pow(error, -1.0 / BT::data.order);
      dt_adaptive = fmin(dt * BT::rk_safety, dt_adaptive);
      dt_adaptive = fmax(dt / BT::rk_safety, dt_adaptive);

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
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::HandleBoundaries(void)
{
   int bnd, bnd_status, distro;

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
         if (action >= 0) distributions[distro]->ProcessTrajectory(coords0, fields0, records.GetMagInitial(), _coords, _fields, records.GetMagExtrema(), action);
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

            if constexpr (HConfig::build_mode == BuildMode::debug) {
               std::cerr << "Position and momentum before reflection: " << _coords.Pos() << " " << _coords.Mom() << std::endl;
            }

// This is a very crude way to do a reflection. In the future one could improve on it by computing the precise boundary crossing location and reflecting the trajectory along the field line from there. However, this requires a lot of extra code.
            _coords.Pos() -= 2.0 * bcond_s[bnd]->GetDelta() * bcond_s[bnd]->GetNormal();
            ReverseMomentum();
            n_refl++;

            if constexpr (HConfig::build_mode == BuildMode::debug) {
               std::cerr << "Position and momentum after reflection: " << _coords.Pos() << " " << _coords.Mom() << std::endl;
            }

// Recompute the boundary since the position and momentum have changed.
            bcond_s[bnd]->ComputeBoundary(_coords, _fields);
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
         if (action >= 0) distributions[distro]->ProcessTrajectory(coords0, fields0, records.GetMagInitial(), _coords, _fields, records.GetMagExtrema(), action);
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
         if (action >= 0) distributions[distro]->ProcessTrajectory(coords0, fields0, records.GetMagInitial(), _coords, _fields, records.GetMagExtrema(), action);
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
template <typename Background, typename Diffusion>
bool TrajectoryBase<Background, Diffusion>::RKAdvance(void)
{

// Retrieve latest point of the trajectory and store it locally.
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
      CommonFields(_coords, _fields);
      MomentumCorrection();
   };

// A new point is obtained, update records.
   records.Store(_coords);

   return true;
};

/*!
\author Vladimir Florinski
\date 09/30/2022
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::MomentumCorrection(void)
{
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename Background, typename Diffusion>
bool TrajectoryBase<Background, Diffusion>::IsSimulationReady(void) const
{
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
\date 05/27/2022
\param[in] background_in Background object for type recognition
\param[in] container_in  Data container for initializating the background object
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::AddBackground(const DataContainer& container_in)
{
   background = std::make_unique(Background());
   background->ConnectRNG(rng);
   background->SetupObject(container_in);
   
   if (IsSimulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] diffusion_in Diffusion object for type recognitions
\param[in] container_in Data container for initializating the diffusion object
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::AddDiffusion(const DiffusionBase& diffusion_in, const DataContainer& container_in)
{
   diffusion = diffusion_in.Clone();
   diffusion->SetupObject(container_in);

   if (IsSimulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] boundary_in  Boundary object for type recognition
\param[in] container_in Data container for initializating the boundary object
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in)
{
// Time boundary
   if (BITS_RAISED(boundary_in.GetStatus(), BOUNDARY_TIME)) {
      bcond_t.push_back(boundary_in.Clone());
      bcond_t.back()->SetupObject(container_in);
   }

// Spatial boundary
   else if (BITS_RAISED(boundary_in.GetStatus(), BOUNDARY_SPACE)) {
      bcond_s.push_back(boundary_in.Clone());
      bcond_s.back()->SetupObject(container_in);
   }

// Momentum boundary
   else if (BITS_RAISED(boundary_in.GetStatus(), BOUNDARY_MOMENTUM)) {
      bcond_m.push_back(boundary_in.Clone());
      bcond_m.back()->SetupObject(container_in);
   }

   if (IsSimulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] initial_in   Initial object for type recognition
\param[in] container_in Data container for initializating the initial object
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::AddInitial(const InitialBase& initial_in, const DataContainer& container_in)
{
// Time condition
   if (BITS_RAISED(initial_in.GetStatus(), INITIAL_TIME)) {
      icond_t = initial_in.Clone();
      icond_t->ConnectRNG(rng);
      icond_t->SetupObject(container_in);
   }

// Spatial condition
   else if (BITS_RAISED(initial_in.GetStatus(), INITIAL_SPACE)) {
      icond_s = initial_in.Clone();
      icond_s->ConnectRNG(rng);
      icond_s->SetupObject(container_in);
   }

// Momentum condition
   else if (BITS_RAISED(initial_in.GetStatus(), INITIAL_MOMENTUM)) {
      icond_m = initial_in.Clone();
      icond_m->ConnectRNG(rng);
      icond_m->SetupObject(container_in);
   };

   if (IsSimulationReady()) RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/08/2024

To start a new trajectory its objects must be set to their initial state. This function determines the initial position and momentum from the respective distributions, calculates the fields, initializes the boundaries at the initial poasition, and resets the counters. A time step evaluation is not performed because it is done in"Advance()" at the beginning of each step.
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::SetStart(void)
try {

// Get the starting time from the initial time distribution
   _coords.Time() = icond_t->GetTimeSample();
// Get the starting position from the initial space distribution.
   _coords.Pos() = icond_s->GetPosSample();
// Get a momentum sample along an arbitrary axis (bhat is unknown at this step). Only the momentum magnitude is needed for the first call to CommonFields().
   _coords.Mom() = icond_m->GetMomSample(gv_ones);

// Obtain the fields for the coordinates (this initializes _dmax).
   CommonFields(_coords, _fields);

// Get the starting momentum from the distribution along the correct axis (bhat is now determined).
   _coords.Mom() = icond_m->GetMomSample(_fields.HatMag());
// Get the starting momentum from the distribution along the correct axis (bhat is now determined).
   _coords.Vel() = Vel<specie>(_coords.Mom());

// Record the initial spatial data for distribution purposes.
   fields0 = _fields;
// Initialize the trajectory records class.
   records.SetStart(_fields);

// Adaptive step must be large at first so that "dt" starts with a physical step.
   dt_adaptive = sp_large * _dmax / c_code;

// Store the initial coordinates (needed for distributions)
   coords0 = _coords;
   
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
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::Integrate(void)
{
   bool was_advanced;
   int time_step_adaptations = 0;

// Time loop is very simple - a single call to "Advance()" followed by a global boundary update. It is the responsibility of Advance() to record the distribution on boundary crossing events.
   while (BITS_LOWERED(_status, TRAJ_FINISH) && BITS_LOWERED(_status, TRAJ_DISCARD)) {

// Attempt to advance trajectory by one segment
      was_advanced = Advance();

      if constexpr (HConfig::trajectory_adv_safety_level > 1) {
// Too many steps were taken - terminate
         if (records.Segments() > HConfig::max_trajectory_steps) {
            RAISE_BITS(_status, TRAJ_DISCARD);
            throw ExMaxStepsReached();
         };

// Too many time adaptations were performed - terminate
         if (was_advanced) time_step_adaptations = 0;
         else {
            time_step_adaptations++;
            if (time_step_adaptations > HConfig::max_time_adaptations) {
               RAISE_BITS(_status, TRAJ_DISCARD);
               throw ExMaxTimeAdaptsReached();
            };
         };
      }

      if constexpr (HConfig::trajectory_adv_safety_level > 0) {
// Time step is too small - terminate
         if (dt < sp_tiny * _dmax / c_code) {
            RAISE_BITS(_status, TRAJ_DISCARD);
            throw ExTimeStepTooSmall();
         }
         else if (!std::isnormal(dt)) {
            RAISE_BITS(_status, TRAJ_DISCARD);
            throw ExTimeStepNan();
         };
      }

   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 02/17/2023
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::StopBackground(void)
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
template <typename Background, typename Diffusion>
int TrajectoryBase<Background, Diffusion>::Crossings(unsigned int output, unsigned int bnd) const
{
   if (bnd < 0) return 0;

   if ((output == 0) && (bnd < bcond_t.size())) return bcond_t[bnd]->CrossingsMade();
   else if ((output == 1) && (bnd < bcond_s.size())) return bcond_s[bnd]->CrossingsMade();
   else if ((output == 2) && (bnd < bcond_m.size())) return bcond_m[bnd]->CrossingsMade();
   return 0;
};


/*!
\author Vladimir Florinski
\date 12/27/2021
*/
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::InterpretStatus(void) const
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
template <typename Background, typename Diffusion>
void TrajectoryBase<Background, Diffusion>::PrintInfo(void) const
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
