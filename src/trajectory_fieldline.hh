/*!
\file trajectory_fieldline.hh
\brief Declares a class for trajectory following a field line
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_FIELDLINE_HH
#define SPECTRUM_TRAJECTORY_FIELDLINE_HH

#include "trajectory_fieldline_base.hh"
#include "common/fields.hh"

namespace Spectrum {

/*!
\brief Field line tracer
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

The type Field_t is a field type that will be traced during the simulation,
i.e., the "field" in Fieldline.
Components of "traj_mom" are: unused (x), unused (y), p_para (z)
*/
template <typename HConfig_, typename Field_t_>
class TrajectoryFieldline : public TrajectoryFieldlineBase<TrajectoryFieldline<HConfig_, Field_t_>, HConfig_> {

//! Readable name
   static constexpr std::string_view traj_name = std::string_view("TrajectoryFieldline" + std::string(Field_t_::name));

public:

   using HConfig = HConfig_;
   using Coordinates = HConfig::Coordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<TrajectoryFocused<HConfig>, HConfig>;
   using HConfig::specie;

   using Field_t = Field_t_;
   using TrajectoryFieldlineBase = TrajectoryFieldlineBase<TrajectoryFieldline<HConfig, Field_t>, HConfig>;

protected:

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;

   using TrajectoryFieldlineBase::SetStart;
   using TrajectoryFieldlineBase::Slopes;
   using TrajectoryFieldlineBase::ConvertMomentum;
   using TrajectoryFieldlineBase::PhysicalStep;
   using TrajectoryFieldlineBase::Advance;

// static assert(s) for Field_t
   static_assert(TrajectoryFields::template found<Field_t>(), "The trace field for TrajectoryFieldline is not a tracked field. Add it to the Fields type defined during configuration.");
   static_assert((std::same_as<Field_t, Vel_t> || std::same_as<Field_t, Mag_t> || std::same_as<Field_t, Elc_t>), "The trace field for TrajectoryFieldline is not supported by the implementation. Choose another field, or else modify the implementation.");

protected:

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

public:

//! Default constructor
   TrajectoryFieldline(void);

//! Copy constructor (class not copyable)
   TrajectoryFieldline(const TrajectoryFieldline& other) = delete;

//! Destructor
   ~TrajectoryFieldline() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryFieldline);

};

};

// Something like this is needed for templated classes
#include "trajectory_fieldline.cc"

#endif
