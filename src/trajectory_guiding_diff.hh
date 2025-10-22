/*!
\file trajectory_guiding_diff.hh
\brief Declares a class for trajectory based on guiding center equations with perpendicular diffusion
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_DIFF_HH
#define SPECTRUM_TRAJECTORY_GUIDING_DIFF_HH

#include "trajectory_guiding.hh"
#include "diffusion_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingDiff class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Trajectory tracer for the relativistic guiding center equations plus perpendicular diffusion
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
template <typename HConfig_, typename Background_>
class TrajectoryGuidingDiff : virtual public TrajectoryGuiding<HConfig_, Background_> {

//! Readable name
   static constexpr std::string_view traj_name = "TrajectoryGuidingDiff";

public:

   using HConfig = HConfig_;
   using Background = Background_;
   using TrajectoryCoordinates = HConfig::TrajectoryCoordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<Background, Diffusion>;
   using HConfig::specie;

   using TrajectoryGuiding = TrajectoryGuiding<Background, Diffusion>;

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;

   using typename TrajectoryBase::DiffusionCoordinates;
   using typename TrajectoryBase::DiffusionFields;
   using typename TrajectoryBase::DiffusionFieldsRemainder;
   using CommonFields_Diffusion = TrajectoryBase::template CommonFields<DiffusionCoordinates, DiffusionFields, DiffusionFieldsRemainder>;

   using TrajectoryBase::diffusion;
   using TrajectoryBase::background;
   using TrajectoryBase::rng;
   using TrajectoryBase::records;
   using TrajectoryBase::SpaceTerminateCheck;
   using TrajectoryBase::slope_pos;
   using TrajectoryBase::slope_mom;
   using TrajectoryBase::Load;
//   using TrajectoryBase::Store;
   using TrajectoryBase::StoreLocal;
   using TrajectoryBase::CommonFields;
   using TrajectoryBase::TimeBoundaryProximityCheck;
   using TrajectoryBase::RKSlopes;
   using TrajectoryBase::RKStep;
   using TrajectoryBase::HandleBoundaries;

   using TrajectoryBase::local_coords;
   using TrajectoryGuiding::drift_vel;

//   using TrajectoryGuiding::ConvertMomentum;
   using TrajectoryGuiding::MomentumCorrection;


protected:

//! The field aligned coordinate system unit vectors (transient)
   GeoVector fa_basis[3];

//! Perpendicular diffusion coefficient (transient)
   double Dperp;

//! Gradient of Dperp (transient)
   GeoVector Vperp;

//! Stochastic increment in guiding center
   GeoVector dr_perp;

//! Compute the basis vectors of the field aligned frame
   void FieldAlignedFrame(void);

//! Computes the perpendicular diffusion coefficient
   virtual void DiffusionCoeff(void);

//! Computes perpendicular diffusion slopes using Euler method
   void EulerPerpDiffSlopes(void);

//! Computes perpendicular diffusion slopes using Milstein method
   void MilsteinPerpDiffSlopes(void);

//! Computes perpendicular diffusion slopes using RK2 method
   bool RK2PerpDiffSlopes(void);

//! Compute the RK(?) slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   bool IsSimulationReady(void) const override;

public:

//! Default constructor
   TrajectoryGuidingDiff(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuidingDiff(const std::string& name_in, uint16_t status_in);

//! Destructor
   ~TrajectoryGuidingDiff() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuidingDiff);

};


};

// Something like this is needed for templated classes
#include "trajectory_guiding_diff.cc"

#endif
