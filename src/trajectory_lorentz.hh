/*!
\file trajectory_lorentz.hh
\brief Declares a class for trajectory based on Newton-Lorentz equation
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_LORENTZ_HH
#define SPECTRUM_TRAJECTORY_LORENTZ_HH

#include "trajectory_base.hh"
#include "common/fields.hh"

namespace Spectrum {

/*!
\brief Trajectory tracer for the Newton-Lorentz equation (full orbit)
\author Vladimir Florinski
*/
template <typename HConfig_>
class TrajectoryLorentz : public TrajectoryBase<TrajectoryLorentz<HConfig_>, HConfig_> {

//! Readable name
   static constexpr std::string_view traj_name = "TrajectoryLorentz";

public:

   using HConfig = HConfig_;
   using Coordinates = HConfig::Coordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<TrajectoryFocused<HConfig>, HConfig>;
   using HConfig::specie;

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;
//   // methods
   using TrajectoryBase::RKAdvance;

   static_assert(!TrajectoryFields::template found<AbsMag_t>(), "AbsMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(!TrajectoryFields::template found<Elc_t>(), "Elc must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");

protected:

//! Previous value of parallel momentum (transient)
   double p_para;

//! Time of the last mirror point (transient)
   double t_mirror;

//! Conversion from (p_x,p_y,p_z) to (p,mu,phi)
   GeoVector ConvertMomentum(void) const override;

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

public:

//! Default constructor
   TrajectoryLorentz(void);

//! Copy constructor (class not copyable)
   TrajectoryLorentz(const TrajectoryLorentz& other) = delete;

//! Destructor
   ~TrajectoryLorentz() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryLorentz);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryLorentz inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/02/2022
\return A vector in the (p,mu,phi) format
*/
template <typename Fields>
inline GeoVector TrajectoryLorentz<Fields>::ConvertMomentum(void) const
{
   return GeoVector(_coords.Mom().Norm(), (_coords.Mom() * _fields.HatMag()) / _coords.Mom().Norm(), 0.0);
};


};

// Something like this is needed for templated classes
#include "trajectory_lorentz.cc"

#endif
