/*!
\file trajectory_fieldline.hh
\brief Declares a class for trajectory following a field line
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_FIELDLINE_BASE_HH
#define SPECTRUM_TRAJECTORY_FIELDLINE_BASE_HH

#include "trajectory_base.hh"

namespace Spectrum {

/*!
\brief Field line tracer base class
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

Base class for any Field line tracer class. The existence of this class
allows compile-time checks almost everywhere, without requiring
that any class except the Field line tracer class itself have any knowledge
of which is the traced field.
Components of "traj_mom" are: unused (x), unused (y), p_para (z)
*/
template <typename Trajectory_, typename HConfig_>
class TrajectoryFieldlineBase : public TrajectoryBase<Trajectory_, HConfig_> {

   static constexpr std::string_view  traj_name = "TrajectoryFieldlineBase";

public:

   using HConfig = HConfig_;
   using Coordinates = HConfig::Coordinates;
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

protected:

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
   TrajectoryFieldlineBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryFieldlineBase(const std::string& name_in, uint16_t status_in);

//! Copy constructor (class not copyable)
   TrajectoryFieldlineBase(const TrajectoryFieldlineBase& other) = delete;

//! Destructor
   ~TrajectoryFieldlineBase() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryFieldlineBase);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFieldline inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/30/2022
\return A vector in the (p,mu,phi) format
\note Not used, but needs to be "overriden" from virtual definition in TrajectoryBase
*/
template <typename Trajectory, typename HConfig>
inline GeoVector TrajectoryFieldlineBase<Trajectory, HConfig>::ConvertMomentum(void) const
{
   return GeoVector(fabs(_coords.Mom()[2]), 0.0, 0.0);
};


};

// Something like this is needed for templated classes
#include "trajectory_fieldline_base.cc"

#endif
