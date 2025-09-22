/*!
\file trajectory_guiding_scatt.hh
\brief Declares a class for trajectory based on guiding center equations with pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_SCATT_HH
#define SPECTRUM_TRAJECTORY_GUIDING_SCATT_HH

#include "trajectory_guiding_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Trajectory tracer for the relativistic guiding center equations plus pitch angle scattering
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
template <typename HConfig_>
class TrajectoryGuidingScatt : public TrajectoryGuidingBase<TrajectoryGuidingScatt<HConfig_>, HConfig_> {

//! Readable name of the TrajectoryGuidingScatt class
   static constexpr std::string_view traj_name = "TrajectoryGuidingScatt";

public:

   using HConfig = HConfig_;
   using Coordinates = HConfig::Coordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<TrajectoryFocused<HConfig>, HConfig>;
   using HConfig::specie;

   using TrajectoryGuidingBase = TrajectoryGuidingBase<TrajectoryGuidingScatt<HConfig>, HConfig>;

protected:

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;

//   using TrajectoryBase::_status;
//   using TrajectoryBase::_t;
//   using TrajectoryBase::_pos;
//   using TrajectoryBase::_mom;
//   using TrajectoryBase::_vel;
//   using TrajectoryBase::dt;
//   using TrajectoryBase::dt_physical;
//   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::rng;
//   using TrajectoryBase::_fields;
//   using TrajectoryBase::_dmax;
////   using TrajectoryBase::traj_t;
////   using TrajectoryBase::traj_pos;
////   using TrajectoryBase::traj_mom;
//   using TrajectoryBase::specie;
////   using TrajectoryBase::local_t;
////   using TrajectoryBase::local_pos;
////   using TrajectoryBase::local_mom;
   using TrajectoryBase::diffusion;
   using TrajectoryBase::slope_pos;
   using TrajectoryBase::slope_mom;
//   // methods:
   using TrajectoryBase::Load;
   using TrajectoryBase::Store;
   using TrajectoryBase::TimeBoundaryProximityCheck;
   using TrajectoryBase::StoreLocal;
   using TrajectoryBase::RKSlopes;
   using TrajectoryBase::RKStep;
   using TrajectoryBase::HandleBoundaries;
   using TrajectoryBase::CommonFields;
   using TrajectoryBase::MomentumCorrection;
   using TrajectoryBase::SpaceTerminateCheck;

   using TrajectoryGuidingBase::DriftCoeff;
   using TrajectoryGuidingBase::Slopes;
   using TrajectoryGuidingBase::ConvertMomentum;


protected:

//! Pitch angle scattering coefficient (transient)
   double Dmumu;

//! Rate of change of Dmumu with mu (transient)
   double Vmu;

#ifdef GEO_DEBUG
//! Number of times |mu| > 1
   int Nabsmugt1 = 0;
#endif

//! Computes the pitch angle diffusion coefficient
   virtual void DiffusionCoeff(void);

//! Performs Euler pitch angle scattering (weak order 1, strong order 1/2)
   void EulerPitchAngleScatt(bool second);

//! Performs Milstein pitch angle scattering (weak order 1, strong order 1)
   void MilsteinPitchAngleScatt(bool second);

//! Performs RK2 pitch angle scattering (weak order 2, strong order 1?)
   void RK2PitchAngleScatt(bool second);

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   bool IsSimulationReady(void) const override;

public:

//! Default constructor
   TrajectoryGuidingScatt(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuidingScatt(const std::string& name_in, uint16_t status_in);

//! Destructor
   ~TrajectoryGuidingScatt() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuidingScatt);

#ifdef GEO_DEBUG

//! Return number of |mu| = 1 crossings
   int fabsmugt1(void) const;

//! Reset number of |mu| = 1 crossings
   void ResetNabsmugt1(void);

#endif
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef GEO_DEBUG

/*!
\author Juan G Alonso Guzman
\date 05/03/2022
\return Number of |mu| = 1 crossings
*/
template <typename Fields>
inline int TrajectoryGuidingScatt<Fields>::fabsmugt1(void) const
{
   return Nabsmugt1;
};

/*!
\author Juan G Alonso Guzman
\date 05/03/2022
*/
template <typename Fields>
inline void TrajectoryGuidingScatt<Fields>::ResetNabsmugt1(void)
{
   Nabsmugt1 = 0;
};

#endif

};

// Something like this is needed for templated classes
#include "trajectory_guiding_scatt.cc"

#endif
