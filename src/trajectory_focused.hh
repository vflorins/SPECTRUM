/*!
\file trajectory_focused.hh
\brief Declares a class for trajectory based on the focused transport equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_FOCUSED_HH
#define SPECTRUM_TRAJECTORY_FOCUSED_HH

#include "trajectory_base.hh"

namespace Spectrum {

// Switch controlling how to calculate mu. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
#define PPERP_METHOD 1

//! Flag to use gradient and curvature drifts in drift velocity calculation
#define TRAJ_FOCUSED_USE_B_DRIFTS

//! Readable name of the TrajectoryFocused class
const std::string traj_name_focused = "TrajectoryFocused";

//! Default initial size
const unsigned int defsize_focused = 10000;

//! CFL condition for advection
const double cfl_adv_tf = 0.5;

//! Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
const double drift_safety_tf = 0.5;

//! How many time steps to allow before recording a mirror event
const int mirror_thresh_focused = 10;

/*!
\brief Trajectory tracer for the focused transport equation
\author Juan G Alonso Guzman

Components of "traj_mom" are: p_mag (x), mu (y), unused (z)
*/
template <typename Fields_>
class TrajectoryFocused : public TrajectoryBase<Fields_> {

public:

   using Fields = Fields_;
   using TrajectoryBase = TrajectoryBase<Fields>;
   using DistributionBase = DistributionBase<TrajectoryBase>;
   using BackgroundBase = BackgroundBase<TrajectoryBase>;
   using DiffusionBase = DiffusionBase<TrajectoryBase>;

   using TrajectoryBase::_t;
   using TrajectoryBase::_pos;
   using TrajectoryBase::_mom;
   using TrajectoryBase::traj_t;
   using TrajectoryBase::traj_pos;
   using TrajectoryBase::traj_mom;
   using TrajectoryBase::_vel;
   using TrajectoryBase::specie;
   using TrajectoryBase::local_t;
   using TrajectoryBase::local_pos;
   using TrajectoryBase::local_mom;

protected:

//! Magnetic moment (transient)
   double mag_mom;

//! Drift elocity (transient)
   GeoVector drift_vel;

//! mu^2 (transient)
   double ct2;

//! 1.0 - mu^2 (transient)
   double st2;

//! Load the last trajectory point
   void Load(void) override;

//! Load the local trajectory point
   void LoadLocal(void) override;

//! Momentum transformation on reflection at a boundary
   void ReverseMomentum(void) override;

//! Compute the drift coefficient
   virtual void DriftCoeff(void);

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Adjust perp momentum to conserve magnetic moment
   void MomentumCorrection(void) override;

public:

//! Default constructor
   TrajectoryFocused(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryFocused(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Copy constructor (class not copyable)
   TrajectoryFocused(const TrajectoryFocused& other) = delete;

//! Destructor
   ~TrajectoryFocused() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryFocused);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFocused inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/06/2023
*/
template <typename Fields>
inline void TrajectoryFocused<Fields>::Load(void)
{
#ifdef RECORD_TRAJECTORY
   _t = traj_t.back();
   _pos = traj_pos.back();
   _mom = traj_mom.back();
#endif
   _vel[0] = Vel(_mom[0], specie);
   _vel[1] = _mom[1];
   _vel[2] = 0.0;
};
/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/06/2023
*/
template <typename Fields>
inline void TrajectoryFocused<Fields>::LoadLocal(void)
{
   _t = local_t;
   _pos = local_pos;
   _mom = local_mom;
   _vel[0] = Vel(_mom[0], specie);
   _vel[1] = _mom[1];
   _vel[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
*/
template <typename Fields>
inline void TrajectoryFocused<Fields>::ReverseMomentum(void)
{
   _mom[1] = -_mom[1];
};

};

// Something like this is needed for templated classes
#include "trajectory_focused.cc"

#endif
