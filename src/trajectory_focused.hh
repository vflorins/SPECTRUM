/*!
\file trajectory_focused.hh
\brief Declares a class for trajectory based on the focused transport equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_FOCUSED_HH
#define SPECTRUM_TRAJECTORY_FOCUSED_HH

#include "trajectory_base.hh"
#include "common/fields.hh"

namespace Spectrum {

/*!
\brief Trajectory tracer for the focused transport equation
\author Juan G Alonso Guzman

Components of "traj_mom" are: p_mag (x), mu (y), unused (z)
*/
template <typename HConfig_>
class TrajectoryFocused : public TrajectoryBase<TrajectoryFocused<HConfig_>, HConfig_> {

//! Readable name
   static constexpr std::string_view traj_name = "TrajectoryFocused";

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

   // methods:
   using TrajectoryBase::RKAdvance;

   using TrajectoryBase::local_coords;

   static_assert(TrajectoryFields::template found<Vel_t>(), "Vel must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(TrajectoryFields::template found<AbsMag_t>(), "AbsMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(TrajectoryFields::template found<HatMag_t>(), "HatMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(TrajectoryFields::template found<DelAbsMag_t>(), "DelAbsMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(TrajectoryFields::template found<DelMag_t>(), "DelMag must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(TrajectoryFields::template found<DelVel_t>(), "DelVel must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");
   static_assert(TrajectoryFields::template found<DotVel_t>(), "DotVel must be tracked by the Trajectory. Add it to the Fields type defined during configuration.");

protected:

//! Magnetic moment (transient)
// todo Mom() and Mom_t is momentum...
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
   TrajectoryFocused(const std::string& name_in, uint16_t status_in);

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
   _coords.Vel()[0] = Vel(_coords.Mom()[0], specie);
   _coords.Vel()[1] = _coords.Mom()[1];
   _coords.Vel()[2] = 0.0;
};
/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/06/2023
*/
template <typename Fields>
inline void TrajectoryFocused<Fields>::LoadLocal(void)
{
   _coords = local_coords;
   _coords.Vel()[0] = Vel(_coords.Mom()[0], specie);
   _coords.Vel()[1] = _coords.Mom()[1];
   _coords.Vel()[2] = 0.0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/07/2023
*/
template <typename Fields>
inline void TrajectoryFocused<Fields>::ReverseMomentum(void)
{
   _coords.Mom()[1] = -_coords.Mom()[1];
};

};

// Something like this is needed for templated classes
#include "trajectory_focused.cc"

#endif
