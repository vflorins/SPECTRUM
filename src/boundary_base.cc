/*!
\file boundary_base.cc
\brief Implements a base class to capture boundary crossing events
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "boundary_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/01/2020
*/
template <typename Trajectory>
BoundaryBase<Trajectory>::BoundaryBase(void)
            : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 12/14/2020
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
template <typename Trajectory>
BoundaryBase<Trajectory>::BoundaryBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
            : Params(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 12/27/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBoundary()" with the argument of "true".
*/
template <typename Trajectory>
BoundaryBase<Trajectory>::BoundaryBase(const BoundaryBase& other)
            : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 12/27/2021
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupBoundary()" method.
*/
template <typename Trajectory>
void BoundaryBase<Trajectory>::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupBoundary(false);
};

/*!
\author Vladimir Florinski
\date 04/13/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Trajectory>
void BoundaryBase<Trajectory>::SetupBoundary(bool construct)
{
// Only needed in the parent version
   container.Reset();
   container.Read(max_crossings);
   container.Read(actions);
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);
   LOWER_BITS(_status, BOUNDARY_TERMINAL);
   LOWER_BITS(_status, BOUNDARY_CROSSED);
};

/*!
\author Vladimir Florinski
\date 12/27/2021
\note This is only a stub.
*/
template <typename Trajectory>
void BoundaryBase<Trajectory>::EvaluateBoundary(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/10/2025
\param[in] t_in    Time
\param[in] pos_in  Position
\param[in] mom_in  Momentum
\param[in] fields_in Fields
\note This is a common routine that the derived classes should not change.
*/
template <typename Trajectory>
void BoundaryBase<Trajectory>::ComputeBoundary(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const Fields& fields_in)
{
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

// Set internal coordinates and evaluate boundary
   SetState(t_in, pos_in, mom_in);
   if constexpr (Fields::AbsMag_found())
      _fields.AbsMag() = fields_in.AbsMag();
   if constexpr (Fields::Iv0_found())
      _fields.Iv0() = fields_in.Iv0();
   if constexpr (Fields::Iv1_found())
      _fields.Iv1() = fields_in.Iv1();
   if constexpr (Fields::Iv2_found())
      _fields.Iv2() = fields_in.Iv2();
   EvaluateBoundary();
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExBoundaryError();

// Change the "BOUNDARY_CROSSED" flag if needed
   if (_delta * _delta_old < 0.0) RAISE_BITS(_status, BOUNDARY_CROSSED);
   else LOWER_BITS(_status, BOUNDARY_CROSSED);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/10/2025
\param[in] t_in   Time
\param[in] pos_in Position
\param[in] mom_in Momentum
\param[in] fields_in Fields
\note This is a common routine that the derived classes should not change.
*/
template <typename Trajectory>
void BoundaryBase<Trajectory>::ResetBoundary(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const Fields& fields_in)
{
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

// Erase all crossing records and the crossing count
   _crossings_left = max_crossings;
   if (_crossings_left == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
   else LOWER_BITS(_status, BOUNDARY_TERMINAL);
   _cross_t.clear();

// Determine the state of the boundary for the given position and set the initial "_delta_old"
   SetState(t_in, pos_in, mom_in);
   if constexpr (Fields::AbsMag_found())
      _fields.AbsMag() = fields_in.AbsMag();
   if constexpr (Fields::Iv0_found())
      _fields.Iv0() = fields_in.Iv0();
   if constexpr (Fields::Iv1_found())
      _fields.Iv1() = fields_in.Iv1();
   if constexpr (Fields::Iv2_found())
      _fields.Iv2() = fields_in.Iv2();
   EvaluateBoundary();
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExBoundaryError();
   _delta_old = _delta;

// The initial point may be right on the boundary, so the code will not be able to determine whether a crossing occurred. For this resaon "_delta_old" is set to a small negative value. This will work for an _external_ boundary, which is the most common case. TODO 
   if (_delta_old == 0.0) _delta_old = -sp_tiny * delta_scale;
   LOWER_BITS(_status, BOUNDARY_CROSSED);
};
   
/*!
\author Vladimir Florinski
\date 04/15/2022
*/
template <typename Trajectory>
void BoundaryBase<Trajectory>::RecordBoundary(void)
{
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

// This is the only place where crossings are recorded. The flag "BOUNDARY_CROSSED" is temporary, and must be now lowered (the event, if it occurred, is now considered to have passed).
   LOWER_BITS(_status, BOUNDARY_CROSSED);
   if (_delta * _delta_old < 0.0) {
      _crossings_left--;
      _cross_t.push_back(_t);

// For a recurrent boundary, we want "_delta_old" to always be negative, so an update doesn't happen on a change of sign of "_delta"
      if (BITS_RAISED(_status, BOUNDARY_RECURRENT)) EvaluateBoundary();
      else _delta_old = _delta;
   }
   else _delta_old = _delta;

   if (_crossings_left == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
   else LOWER_BITS(_status, BOUNDARY_TERMINAL);
};

};
