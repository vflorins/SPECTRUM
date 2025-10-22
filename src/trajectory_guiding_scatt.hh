/*!
\file trajectory_guiding_scatt.hh
\brief Declares a class for trajectory based on guiding center equations with pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_SCATT_HH
#define SPECTRUM_TRAJECTORY_GUIDING_SCATT_HH

#include "trajectory_guiding.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Trajectory tracer for the relativistic guiding center equations plus pitch angle scattering
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
template <typename Background_, typename Diffusion_>
class TrajectoryGuidingScatt : public TrajectoryGuiding<HConfig_, Background_> {

//! Readable name of the class
   static constexpr std::string_view traj_name = "TrajectoryGuidingScatt";

public:

   using HConfig = HConfig_;
   using Background = Background_;
   using TrajectoryCoordinates = HConfig::TrajectoryCoordinates;
   using TrajectoryFields = HConfig::TrajectoryFields;
   using TrajectoryBase = TrajectoryBase<Background, Diffusion>;
   using HConfig::specie;
   using HConfig::TrajectoryConfig::split_scatt;
   using HConfig::TrajectoryConfig::split_scatt_fraction;
   using HConfig::TrajectoryConfig::stochastic_method_mu;
   using HConfig::TrajectoryConfig::const_dmumax;
   using HConfig::TrajectoryConfig::cfl_pitchangle;

   using TrajectoryGuiding = TrajectoryGuiding<Background, Diffusion>;

protected:

   using TrajectoryBase::_status;
   using TrajectoryBase::_coords;
   using TrajectoryBase::_fields;
   using TrajectoryBase::_dmax;
   using TrajectoryBase::dt;
   using TrajectoryBase::dt_adaptive;
   using TrajectoryBase::dt_physical;
   using TrajectoryBase::records;

   using typename TrajectoryBase::DiffusionCoordinates;
   using typename TrajectoryBase::DiffusionFields;
   using typename TrajectoryBase::DiffusionFieldsRemainder;
   using CommonFields_Diffusion = TrajectoryBase::template CommonFields<DiffusionCoordinates, DiffusionFields, DiffusionFieldsRemainder>;

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
//   using TrajectoryBase::Store;
   using TrajectoryBase::TimeBoundaryProximityCheck;
   using TrajectoryBase::StoreLocal;
   using TrajectoryBase::RKSlopes;
   using TrajectoryBase::RKStep;
   using TrajectoryBase::HandleBoundaries;
   using TrajectoryBase::MomentumCorrection;
   using TrajectoryBase::SpaceTerminateCheck;

   using TrajectoryGuiding::DriftCoeff;
   using TrajectoryGuiding::Slopes;
//   using TrajectoryGuiding::ConvertMomentum;


protected:

//! Pitch angle scattering coefficient (transient)
   double Dmumu;

//! Rate of change of Dmumu with mu (transient)
   double Vmu;

//! Number of times |mu| > 1
   int Nabsmugt1 = 0;

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

//! Return number of |mu| = 1 crossings
   int fabsmugt1(void) const;

//! Reset number of |mu| = 1 crossings
   void ResetNabsmugt1(void);

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
template <typename HConfig>
inline int TrajectoryGuidingScatt<HConfig>::fabsmugt1(void) const
{
   return Nabsmugt1;
};

/*!
\author Juan G Alonso Guzman
\date 05/03/2022
*/
template <typename HConfig>
inline void TrajectoryGuidingScatt<HConfig>::ResetNabsmugt1(void)
{
   Nabsmugt1 = 0;
};

#endif

};

// Something like this is needed for templated classes
#include "trajectory_guiding_scatt.cc"

#endif
