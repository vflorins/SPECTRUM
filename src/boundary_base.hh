/*!
\file boundary_base.hh
\brief Declares a base class to capture boundary crossing events
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BOUNDARY_BASE_HH
#define SPECTRUM_BOUNDARY_BASE_HH

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "config.h"
#include "common/params.hh"
#include "common/physics.hh"
#include <memory>

namespace Spectrum {

//! Boundary is in time
const uint16_t BOUNDARY_TIME = 0x0010;

//! Boundary is spatial
const uint16_t BOUNDARY_SPACE = 0x0020;

//! Boundary is momentum boundary
const uint16_t BOUNDARY_MOMENTUM = 0x0040;

//! Reflecting boundary, trajectory will continue
const uint16_t BOUNDARY_REFLECT = 0x0080;

//! Trajectory will stop after next crossing
const uint16_t BOUNDARY_TERMINAL = 0x0100;

//! Boundary was crossed since last saved position (delta's have opposite signs)
const uint16_t BOUNDARY_CROSSED = 0x0200;

//! Boundary is recurrent (moves by a certain amount after each crossing)
const uint16_t BOUNDARY_RECURRENT = 0x0400;

//! Types of boundary conditions
const std::string boundary_type_names[3] = {"time", "space", "momentum"};

//! Clone function pattern
#define CloneFunctionBoundary(T) std::unique_ptr<BoundaryBase> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if boundary evaluation failed
\author Vladimir Florinski
*/
class ExBoundaryError : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\return Text describing the error
*/
inline const char* ExBoundaryError::what(void) const noexcept
{
   return "Boundary evaluation error";
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Base class to perform a test for a boundary crossing event
\author Vladimir Florinski

The "BoundaryXXXXX" classes describe generic boundary condisions in time, space, or momentum. An absorbing boundary should cause the trajectory to stop. A reflecting or event boundary records the times of each reflection or crossing for further processing. One specific case is the mirror event counter that records each reversal of the parallel momentum component.

Parameters: int max_crossings, VEC_INT actions
*/
class BoundaryBase : public Params {

protected:

//! Total number of crossings allowed, -1 means unlimited (persistent)
   int max_crossings;

//! Actions to perform when recording events, -1 for no action (persistent)
   std::vector<int> actions;

//! A typical value of "_delta", typically the expected dt, dx, or dp (transient)
   double delta_scale;

//! Normal vector, used to determine velocity after reflection (transient)
   GeoVector _normal;

//! Number of crossings left (transient)
   int _crossings_left;

//! Some measure, specific to a derived class, of the distance to the boundary; value can be positive or negative (transient)
   double _delta;

//! Previous value of delta to monitor for change in sign (transient)
   double _delta_old;

//! Crossing times stored here (transient)
   std::vector<double> _cross_t;

//! Magnetic field direction (transient)
   GeoVector bhat;

//! Region vector (transient)
   GeoVector region;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryBase(const BoundaryBase& other);

//! Set up the boundary crossing tester based on "params"
   virtual void SetupBoundary(bool construct);

//! Compute the distance to the boundary
   virtual void EvaluateBoundary(void);

public:

//! Destructor
   virtual ~BoundaryBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<BoundaryBase> Clone(void) const = 0;

//! Set up the class parameters
   void SetupObject(const DataContainer& cont_in);

//! Assign action
   void SetAction(int action_in);

//! Retrieve action
   int GetAction(int distro) const;

//! Set the time, space, or momentum scale of the problem
   void SetScale(double scale_in);

//! Return number of crossings that occurred
   int CrossingsMade(void) const;

//! Return number of crossings remaining
   int CrossingsLeft(void) const;

//! Decrement number of crossings left by 1
   void DecrCrossingsLeft(void);

//! Reset the crossings count and set the initial delta. The function should be called once, before starting a new trajectory.
   void ResetBoundary(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const GeoVector& bhat_in, const GeoVector& region_in);

//! Set the state and evaluate the boundary. Call this near the end of each time step, before the BC handler.
   void ComputeBoundary(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const GeoVector& bhat_in, const GeoVector& region_in);

// Update the state of the boundary with the current PSCs. Should be used at the very end of each time step.
   void RecordBoundary(void);

//! Return current value of "_delta"
   double GetDelta(void) const;

//! Return current value of "_normal"
   GeoVector GetNormal(void) const;

//! Return the timestamp of a selected past crossing
   double GetRecord(int crs) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBase inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/25/2022
\param[in] action_in Action to assign
*/
inline void BoundaryBase::SetAction(int action_in)
{
   actions.push_back(action_in);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/25/2022
\return Current action
*/
inline int BoundaryBase::GetAction(int distro) const
{
   return actions[distro];
};

/*!
\author Vladimir Florinski
\date 04/13/2022
\params[in] scale_in Typical time, space or momentum scale of the problem
*/
inline void BoundaryBase::SetScale(double scale_in)
{
   delta_scale = fabs(scale_in);
};

/*!
\author Vladimir Florinski
\date 01/26/2021
\return Number of crossings
*/
inline int BoundaryBase::CrossingsMade(void) const
{
   return max_crossings - _crossings_left;
};

/*!
\author Vladimir Florinski
\date 01/26/2021
\return Number of crossings remaining

This function could be used to figure out if the traectory should terminate ("_crossings_left" is zero). For negative "_max_crossings" the value of "_crossings_left" will be negative until it wraps around after some 2 billion crossings (so a test for zero might succeed, however unlikely).
*/
inline int BoundaryBase::CrossingsLeft(void) const
{
   return _crossings_left;
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
*/
inline void BoundaryBase::DecrCrossingsLeft(void)
{
   _crossings_left--;
};

/*!
\author Vladimir Florinski
\date 02/10/2022
\return Distance to the boundary
*/
inline double BoundaryBase::GetDelta(void) const
{
   return _delta;
};

/*!
\author Vladimir Florinski
\date 02/10/2022
\return Normal to the boundary
*/
inline GeoVector BoundaryBase::GetNormal(void) const
{
   return _normal;
};

/*!
\author Vladimir Florinski
\date 02/08/2021
\params[in] crs Index of the crossing
\return Time of crossing
*/
inline double BoundaryBase::GetRecord(int crs) const
{
   if ((crs < 0) || (std::abs(crs) >= _cross_t.size())) return 0.0;
   else return _cross_t[crs];
};

};

#endif
