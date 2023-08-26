/*!
\file trajectory_parker.hh
\brief Declares a class for trajectory based on the Parker Transport Equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef TRAJECTORY_PARKER_HH
#define TRAJECTORY_PARKER_HH

#include "trajectory_base.hh"

namespace Spectrum {

//! Which stochastic method to use for diffusion, 0 = Euler, 1 = Milstein, 2 = RK2
#define STOCHASTIC_METHOD_DIFF 0

//! Readable name of the TrajectoryParker class
const std::string traj_name_parker = "TrajectoryParker";

//! Default initial size
const unsigned int defsize_parker = 100000;

//! CFL condition for advection
const double cfl_adv_tp = 0.5;

//! Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
const double drift_safety_tp = 0.1;

//! CFL condition for diffusion
const double cfl_dif_tp = 0.5;

//! How many time steps to allow before recording a mirror event
const int mirror_thresh_parker = 10;

/*!
\brief A derived class for Parker equation (diffusive simulation)
\author Juan G Alonso Guzman

Components of "traj_mom" are: p_mag (x), unused (y), unused (z)
*/
class TrajectoryParker : public TrajectoryBase {

//! Gradient of |B| (transient)
   GeoVector gradBmag;

//! Drift elocity (transient)
   GeoVector drift_vel;

//! The field aligned coordinate system unit vectors (transient)
   GeoVector fa_basis[3];

//! Perpendicular diffusion coefficient (transient)
   double Dperp;

//! Parallel diffusion coefficient (transient)
   double Dpara;

//! Gradient of diffusion tensor (transient)
   GeoVector divK;

//! Stochastic increment in guiding center
   GeoVector dr_perp;

//! Compute the basis vectors of the field aligned frame
   void FieldAlignedFrame(void);

//! Computes the diffusion coefficients
   void DiffusionCoeff(void);

//! Computes diffusion slopes using Euler method
   void EulerDiffSlopes(void);

//! Compute the drift coefficient
   void DriftCoeff(void);

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   bool IsSimmulationReady(void) const override;

//! Conversion from (p_x,p_y,p_z) to (p,mu,phi)
   GeoVector ConvertMomentum(void) const override;

//! Momentum transformation on reflection at a boundary
   void ReverseMomentum(void) override;

public:

//! Default constructor
   TrajectoryParker(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryParker(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Copy constructor (class not copyable)
   TrajectoryParker(const TrajectoryParker& other) = delete;

//! Destructor
   ~TrajectoryParker() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryParker);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryParker inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/30/2022
\return A vector in the (p,mu,phi) format
\note Not used, but needs to be "overriden" from virtual definition in TrajectoryBase
*/
inline GeoVector TrajectoryParker::ConvertMomentum(void) const
{
   return _mom;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/07/2023
*/
inline void TrajectoryParker::ReverseMomentum(void)
{
};

//! Trajectory type
#if TRAJ_TYPE == TRAJ_PARKER
typedef TrajectoryParker TrajectoryType;
#endif

};

#endif