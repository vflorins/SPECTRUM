/*!
\file trajectory_parker.hh
\brief Declares a class for trajectory based on the Parker Transport Equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef TRAJECTORY_PARKER_HH
#define TRAJECTORY_PARKER_HH

#include "trajectory_base.hh"
#include "common/fields.hh"

namespace Spectrum {

//! Which stochastic method to use for diffusion: 0 = Euler, 1 = Milstein, 2 = RK2
#define TRAJ_PARKER_STOCHASTIC_METHOD_DIFF 0

//! Flag to use gradient and curvature drifts in drift velocity calculation
//#define TRAJ_PARKER_USE_B_DRIFTS

//! Which method of computation to use for divK: 0 = using direct central FD, 1 = using _spdata.grad quantities
#define TRAJ_PARKER_DIVK_METHOD 0

//! Readable name of the TrajectoryParker class
const std::string traj_name_parker = "TrajectoryParker";

//! Default initial size
const unsigned int defsize_parker = 100000;

//! CFL condition for advection
const double cfl_adv_tp = 0.5;

//! CFL condition for diffusion
const double cfl_dif_tp = 0.5;

//! CFL condition for acceleration
const double cfl_acc_tp = 0.5;

//! Maximum allowed fraction of momentum change per step
const double dlnpmax = 0.01;

/*!
\brief A derived class for Parker equation (diffusive simulation)
\author Juan G Alonso Guzman

Components of "traj_mom" are: p_mag (x), unused (y), unused (z)
*/
template <typename Fields_>
class TrajectoryParker : public TrajectoryBase<TrajectoryParker<Fields_>, Fields_> {
public:

   using Fields = Fields_;
   using BackgroundBase = BackgroundBase<Fields>;
   using TrajectoryBase = TrajectoryBase<TrajectoryParker<Fields_>, Fields_>;

   using TrajectoryBase::_t;
   using TrajectoryBase::_pos;
   using TrajectoryBase::_vel;
   using TrajectoryBase::_mom;
   using TrajectoryBase::_status;
   using TrajectoryBase::dt;
   using TrajectoryBase::rng;
//   using TrajectoryBase::traj_t;
//   using TrajectoryBase::traj_pos;
//   using TrajectoryBase::traj_mom;
   using TrajectoryBase::specie;

//   using TrajectoryBase::local_t;
//   using TrajectoryBase::local_pos;
//   using TrajectoryBase::local_mom;
   using TrajectoryBase::diffusion;
   using TrajectoryBase::background;

   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt_physical;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::slope_pos;
   using TrajectoryBase::slope_mom;
   // methods:
   using TrajectoryBase::Load;
   using TrajectoryBase::Store;
   using TrajectoryBase::StoreLocal;
   using TrajectoryBase::TimeBoundaryProximityCheck;
   using TrajectoryBase::RKSlopes;
   using TrajectoryBase::RKStep;
   using TrajectoryBase::HandleBoundaries;
   using TrajectoryBase::CommonFields;
   using TrajectoryBase::ConnectRNG;

// todo: the list of checks is not exhausive
   static_assert(!Fields::template found<HatMag_t>(), "HatMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(!Fields::template found<AbsMag_t>(), "AbsMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(!Fields::template found<DelAbsMag_t>(), "DelAbsMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");

protected:

//! Drift velocity (transient)
   GeoVector drift_vel;

//! The field aligned coordinate system unit vectors (transient)
   GeoVector fa_basis[3];

//! Perpendicular diffusion coefficient (transient)
   double Kperp;

//! Parallel diffusion coefficient (transient)
   double Kpara;

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
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/07/2023
*/
template <typename Fields>
inline void TrajectoryParker<Fields>::ReverseMomentum(void)
{
};

};

// Something like this is needed for templated classes
#include "trajectory_parker.cc"

#endif
