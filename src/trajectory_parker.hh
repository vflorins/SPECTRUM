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

/*!
\brief A derived class for Parker equation (diffusive simulation)
\author Juan G Alonso Guzman

Components of "traj_mom" are: p_mag (x), unused (y), unused (z)
*/
template <typename Background_, typename Diffusion_>
class TrajectoryParker : public TrajectoryBase<Background_, Diffusion_> {

   //! Readable name
   static constexpr std::string_view traj_name = "TrajectoryParker";

public:

   using Background = Background_;
   using Diffusion = Diffusion_;
   using HConfig = Background::HConfig;
   // todo Traj'yConfig
   using TrajectoryCoordinates = HConfig::TrajectoryCoordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<Background, Diffusion>;
   using HConfig::specie;

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;

   using TrajectoryBase::diffusion;
   using typename TrajectoryBase::DiffusionCoordinates;
   using typename TrajectoryBase::DiffusionFields;
   using typename TrajectoryBase::DiffusionFieldsRemainder;

   using TrajectoryBase::background;
   using TrajectoryBase::rng;
   using TrajectoryBase::slope_pos;
   using TrajectoryBase::slope_mom;

   using TrajectoryBase::Load;
//   using TrajectoryBase::Store;
   using TrajectoryBase::StoreLocal;
   using TrajectoryBase::RKSlopes;
   using TrajectoryBase::RKStep;
   using TrajectoryBase::HandleBoundaries;
   using TrajectoryBase::CommonFields;
   using TrajectoryBase::ConnectRNG;
   using TrajectoryBase::TimeBoundaryProximityCheck;
   using TrajectoryBase::records;

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
   bool IsSimulationReady(void) const override;

//! Momentum transformation on reflection at a boundary
   void ReverseMomentum(void) override;

public:

//! Default constructor
   TrajectoryParker(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryParker(const std::string& name_in, uint16_t status_in);

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
template <typename Background, typename Diffusion>
inline void TrajectoryParker<Background, Diffusion>::ReverseMomentum(void)
{
};

};

// Something like this is needed for templated classes
#include "trajectory_parker.cc"

#endif
