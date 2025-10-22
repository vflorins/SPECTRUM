/*!
\file trajectory_guiding_diff_scatt.hh
\brief Declares a class for trajectory based on guiding center equations with perpendicular diffusion and pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_DIFF_SCATT_HH
#define SPECTRUM_TRAJECTORY_GUIDING_DIFF_SCATT_HH

#include "trajectory_guiding_diff.hh"
#include "trajectory_guiding_scatt.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Trajectory tracer for the relativistic guiding center equations plus perpendicular diffusion and pitch angle scattering
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
template <typename HConfig_, typename Background_>
class TrajectoryGuidingDiffScatt : public TrajectoryGuiding<HConfig_, Background_>, TrajectoryGuidingDiff<HConfig_, Background_>, TrajectoryGuidingScatt<HConfig_, Background_> {

//! Readable name
   static constexpr std::string_view traj_name = "TrajectoryGuidingDiffScatt";

public:

   using HConfig = HConfig_;
   using Background = Background_;
   using TrajectoryCoordinates = HConfig::TrajectoryCoordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<Background, Diffusion>;
   using HConfig::specie;

   using TrajectoryGuiding = TrajectoryGuiding<Background, Diffusion>;
   using TrajectoryGuidingDiff = TrajectoryGuidingDiff<Background, Diffusion>;
   using TrajectoryGuidingScatt = TrajectoryGuidingScatt<Background, Diffusion>;

   using typename TrajectoryBase::DiffusionCoordinates;
   using typename TrajectoryBase::DiffusionFields;
   using typename TrajectoryBase::DiffusionFieldsRemainder;
   using CommonFields_Diffusion = TrajectoryBase::template CommonFields<DiffusionCoordinates, DiffusionFields, DiffusionFieldsRemainder>;

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;

//   using TrajectoryBase::_status;
////   using TrajectoryBase::_t;
//   using TrajectoryBase::_pos;
////   using TrajectoryBase::_mom;
////   using TrajectoryBase::_vel;
//   using TrajectoryBase::dt;
//   using TrajectoryBase::dt_physical;
//   using TrajectoryBase::dt_adaptive;
////   using TrajectoryBase::rng;
////   using TrajectoryBase::_fields;
////   using TrajectoryBase::_dmax;
//////   using TrajectoryBase::traj_t;
//////   using TrajectoryBase::traj_pos;
//////   using TrajectoryBase::traj_mom;
////   using TrajectoryBase::specie;
//////   using TrajectoryBase::local_t;
//////   using TrajectoryBase::local_pos;
//////   using TrajectoryBase::local_mom;
////   using TrajectoryBase::diffusion;
   using TrajectoryBase::slope_pos;
   using TrajectoryBase::slope_mom;
////   // methods:
////   using TrajectoryBase::ConvertMomentum;
   using TrajectoryBase::Load;

////   using TrajectoryBase::TimeBoundaryProximityCheck;
   using TrajectoryBase::StoreLocal;
   using TrajectoryBase::RKSlopes;
   using TrajectoryBase::RKStep;
   using TrajectoryBase::HandleBoundaries;
//   using TrajectoryBase::MomentumCorrection;
   using TrajectoryBase::SpaceTerminateCheck;
   using TrajectoryGuiding::DriftCoeff;
   using TrajectoryGuidingDiff::dr_perp;
   using TrajectoryBase::records;

protected:

//! Additionally computes the pitch angle diffusion coefficient
   void DiffusionCoeff(void) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   bool IsSimulationReady(void) const override;

public:

//! Default constructor
   TrajectoryGuidingDiffScatt(void);

//! Destructor
   ~TrajectoryGuidingDiffScatt() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuidingDiffScatt);
};


};

// Something like this is needed for templated classes
#include "trajectory_guiding_diff_scatt.cc"

#endif
