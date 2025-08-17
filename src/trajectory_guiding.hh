/*!
\file trajectory_guiding.hh
\brief Declares a class for trajectory based on relativistic guiding center equations
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_HH
#define SPECTRUM_TRAJECTORY_GUIDING_HH

#include "trajectory_guiding_base.hh"

namespace Spectrum {

//! Readable name of the TrajectoryGuiding class
const std::string traj_name_guiding = "TrajectoryGuiding";

//! Default initial size
const unsigned int defsize_guiding = 10000;

/*!
\brief Trajectory tracer for the relativistic guiding center equations
\author Vladimir Florinski
\author Lucius Schoenbaum

Components of "traj_mom" are: p_perp (x), unused (y), p_para (z)
*/
template <typename Fields_>
class TrajectoryGuiding : public TrajectoryGuidingBase<TrajectoryGuiding<Fields_>, Fields_> {

public:

   using Fields = Fields_;
   using TrajectoryGuidingBase = TrajectoryGuidingBase<TrajectoryGuiding<Fields>, Fields_>;
   using TrajectoryBase = TrajectoryBase<TrajectoryGuiding<Fields>, Fields>;

   using TrajectoryGuidingBase::_status;
//   using TrajectoryGuidingBase::_t;
   using TrajectoryGuidingBase::_vel;
//   using TrajectoryGuidingBase::_pos;
   using TrajectoryGuidingBase::_mom;
   using TrajectoryGuidingBase::q;
   using TrajectoryGuidingBase::_fields;
   using TrajectoryGuidingBase::_ddata;
//   using TrajectoryGuidingBase::traj_t;
//   using TrajectoryGuidingBase::traj_pos;
   using TrajectoryGuidingBase::traj_mom;
   using TrajectoryGuidingBase::specie;
//   using TrajectoryGuidingBase::local_t;
//   using TrajectoryGuidingBase::local_pos;
//   using TrajectoryGuidingBase::local_mom;
   using TrajectoryGuidingBase::dt_physical;
   using TrajectoryGuidingBase::RKAdvance;

   using TrajectoryGuidingBase::mag_mom;
   using TrajectoryGuidingBase::Evec_star;
   using TrajectoryGuidingBase::Bvec_star;
   using TrajectoryGuidingBase::drift_vel;
   // methods
   using TrajectoryGuidingBase::ConvertMomentum;
   using TrajectoryGuidingBase::ReverseMomentum;
   using TrajectoryGuidingBase::ModifiedFields;
   using TrajectoryGuidingBase::DriftCoeff;
   using TrajectoryGuidingBase::Slopes;
   using TrajectoryGuidingBase::PhysicalStep;
   using TrajectoryGuidingBase::Advance;
   using TrajectoryGuidingBase::MomentumCorrection;
   using TrajectoryGuidingBase::SetStart;

public:

//! Default constructor
   TrajectoryGuiding(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuiding(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Copy constructor (class not copyable)
   TrajectoryGuiding(const TrajectoryGuiding& other) = delete;

//! Destructor
   ~TrajectoryGuiding() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuiding);

};

};

// Something like this is needed for templated classes
#include "trajectory_guiding.cc"

#endif
