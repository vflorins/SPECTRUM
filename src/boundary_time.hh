/*!
\file boundary_time.hh
\brief Declares several classes representing temporal boundaries
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BOUNDARY_TIME_HH
#define SPECTRUM_BOUNDARY_TIME_HH

#include "boundary_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTime class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Time base class
\author Vladimir Florinski

Parameters: (BoundaryBase), double timemark
*/
template <typename HConfig_>
class BoundaryTime : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryTime";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryBase::_delta;

protected:

//! Time from the start (persistent)
   double timemark;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryTime(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryTime(const std::string& name_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryTime(const BoundaryTime& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundaryTime() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTimeExpire class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Time expiration boundary
\author Vladimir Florinski

Parameters: (BoundaryTime)
*/
template <typename HConfig_>
class BoundaryTimeExpire : public BoundaryTime<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryTimeExpire";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryTime = BoundaryTime<HConfig_>;

   using BoundaryBase::max_crossings;

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryTimeExpire(void);

//! Copy constructor
   BoundaryTimeExpire(const BoundaryTimeExpire& other);

//! Destructor
   ~BoundaryTimeExpire() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryTimeExpire);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTimePass class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Describes a fixed time marker
\author Vladimir Florinski

Parameters: (BoundaryTime)
*/
template <typename HConfig_>
class BoundaryTimePass : public BoundaryTime<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryTimePass";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryTime = BoundaryTime<HConfig_>;

   using BoundaryBase::max_crossings;

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryTimePass(void);

//! Copy constructor
   BoundaryTimePass(const BoundaryTimePass& other);

//! Destructor
   ~BoundaryTimePass() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryTimePass);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTimeRecurrent class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Describes a boundary that recurs, i.e., moves forward by a fixed time after each crossing
\author Vladimir Florinski

Parameters: (BoundaryTime)
*/
template <typename HConfig_>
class BoundaryTimeRecurrent : public BoundaryTime<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryTimeRecurrent";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryTime = BoundaryTime<HConfig_>;

   using BoundaryBase::_delta;
   using BoundaryBase::CrossingsMade;

   using BoundaryTime::timemark;

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Default constructor
   BoundaryTimeRecurrent(void);

//! Copy constructor
   BoundaryTimeRecurrent(const BoundaryTimeRecurrent& other);

//! Destructor
   ~BoundaryTimeRecurrent() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryTimeRecurrent);
};

};

// Something like this is needed for templated classes
#include "boundary_time.cc"

#endif
