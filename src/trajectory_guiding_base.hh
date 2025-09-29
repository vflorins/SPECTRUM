/*!
\file trajectory_guiding.hh
\brief Declares a class for trajectory based on relativistic guiding center equations
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_BASE_HH
#define SPECTRUM_TRAJECTORY_GUIDING_BASE_HH

#include "trajectory_base.hh"

namespace Spectrum {

/*!
\brief Trajectory tracer for the relativistic guiding center equations
\author Vladimir Florinski

Components of "traj_mom" are: p_perp (x), unused (y), p_para (z)
*/
template <typename Trajectory_, typename HConfig_>
class TrajectoryGuidingBase : public TrajectoryBase<Trajectory_, HConfig_> {

//! Readable name
   static constexpr std::string_view traj_name = "TrajectoryGuidingBase";

public:

   using HConfig = HConfig_;
   using TrajectoryCoordinates = HConfig::TrajectoryCoordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<TrajectoryFocused<HConfig>, HConfig>;
   using HConfig::specie;

protected:

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;
   using TrajectoryBase::RKAdvance;

   using HConfig::q;

//   using TrajectoryBase::_status;
////   using TrajectoryBase::_t;
//   using TrajectoryBase::_vel;
////   using TrajectoryBase::_pos;
//   using TrajectoryBase::_mom;
//   using TrajectoryBase::q;
//   using TrajectoryBase::_fields;
//   using TrajectoryBase::_dmax;
////   using TrajectoryBase::traj_t;
////   using TrajectoryBase::traj_pos;
//   using TrajectoryBase::traj_mom;
//   using TrajectoryBase::specie;
////   using TrajectoryBase::local_t;
////   using TrajectoryBase::local_pos;
////   using TrajectoryBase::local_mom;
//   using TrajectoryBase::dt_physical;

protected:

//! Magnetic moment (transient)
   double mag_mom;

//! Modified electric field (transient)
   GeoVector Evec_star;

//! Modified magnetic field (transient)
   GeoVector Bvec_star;

//! Drift elocity (transient)
   GeoVector drift_vel;

//! Conversion from (p_perp,0,p_para) to (p,mu,0)
   GeoVector ConvertMomentum(void) const override;

//! Momentum transformation on reflection at a boundary
   void ReverseMomentum(void) override;

//! Compute the modified electromagnetic fields
   void ModifiedFields(void);

//! Compute the drift coefficient
   void DriftCoeff(void);

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
   TrajectoryGuidingBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuidingBase(const std::string& name_in, uint16_t status_in);

//! Copy constructor (class not copyable)
   TrajectoryGuidingBase(const TrajectoryGuidingBase& other) = delete;

//! Destructor
   ~TrajectoryGuidingBase() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuidingBase);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuiding inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/11/2022
\return A vector in the (p,mu,0) format
*/
template <typename Trajectory, typename HConfig>
inline GeoVector TrajectoryGuidingBase<Trajectory, HConfig>::ConvertMomentum(void) const
{
   return GeoVector(sqrt(Sqr(_coords.Mom()[0])+Sqr(_coords.Mom()[2])), _coords.Mom()[2] / _coords.Mom().Norm(), 0.0);
};


/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/07/2023
*/
template <typename Trajectory, typename HConfig>
inline void TrajectoryGuidingBase<Trajectory, HConfig>::ReverseMomentum(void)
{
   _coords.Mom()[2] = -_coords.Mom()[2];
};


};

// Something like this is needed for templated classes
#include "trajectory_guiding_base.cc"

#endif
