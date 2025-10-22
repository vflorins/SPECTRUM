/*!
\file trajectory_fieldline.hh
\brief Declares a class for trajectory following a field line
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_FIELDLINE_HH
#define SPECTRUM_TRAJECTORY_FIELDLINE_HH

#include "trajectory_base.hh"
#include "common/fields.hh"


namespace Spectrum {

/*!
\brief Field line tracer base class
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum

The type Field_t is a field type that will be traced during the simulation,
i.e., the "field" in Fieldline.
Components of "traj_mom" are: unused (x), unused (y), p_para (z)
*/
template <typename HConfig_, typename Background_, typename Field_t_>
class TrajectoryFieldline : public TrajectoryBase<HConfig_, Background_> {

   static constexpr std::string_view traj_name = "TrajectoryFieldline";

public:

   using HConfig = HConfig_;
   using Background = Background_;
   using Field_t = Field_t_;
   using TrajectoryFields = Fields<FConfig<>, Field_t>;
   using TrajectoryCoordinates = HConfig::TrajectoryCoordinates;
   using TrajectoryBase = TrajectoryBase<HConfig, Background_>;
   using HConfig::specie;

   static_assert((std::same_as<Field_t, Vel_t> || std::same_as<Field_t, Mag_t> || std::same_as<Field_t, Elc_t>), "The trace field for TrajectoryFieldline is not supported by the implementation. Choose another field, or else modify the implementation.");

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
//   GeoVector ConvertMomentum(void) const override;

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

public:

//! Default constructor
   TrajectoryFieldline(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryFieldline(const std::string& name_in, uint16_t status_in);

//! Copy constructor (class not copyable)
   TrajectoryFieldline(const TrajectoryFieldline& other) = delete;

//! Destructor
   ~TrajectoryFieldline() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryFieldline);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;


};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryFieldline inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

///*!
//\author Vladimir Florinski
//\date 09/30/2022
//\return A vector in the (p,mu,phi) format
//\note Not used, but needs to be "overriden" from virtual definition in TrajectoryBase
//*/
//template <typename HConfig>
//inline GeoVector TrajectoryFieldlineBase<HConfig>::ConvertMomentum(void) const
//{
//   return GeoVector(fabs(_coords.Mom()[2]), 0.0, 0.0);
//};


};

// Something like this is needed for templated classes
#include "trajectory_fieldline.cc"

#endif
