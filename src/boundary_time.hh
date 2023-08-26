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
class BoundaryTime : public BoundaryBase {

protected:

//! Time from the start (persistent)
   double timemark;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryTime(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryTime(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

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

//! Readable name of the BoundaryTimeExpire class
const std::string bnd_name_time_expire = "BoundaryTimeExpire";

/*!
\brief Time expiration boundary
\author Vladimir Florinski

Parameters: (BoundaryTime)
*/
class BoundaryTimeExpire : public BoundaryTime {

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

//! Readable name of the BoundaryTimePass class
const std::string bnd_name_time_pass = "BoundaryTimePass";

/*!
\brief Describes a fixed time marker
\author Vladimir Florinski

Parameters: (BoundaryTime)
*/
class BoundaryTimePass : public BoundaryTime {

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

//! Readable name of the BoundaryTimeRecurrent class
const std::string bnd_name_time_recur = "BoundaryTimeRecurrent";

/*!
\brief Describes a boundary that recurs, i.e., moves forward by a fixed time after each crossing
\author Vladimir Florinski

Parameters: (BoundaryTime)
*/
class BoundaryTimeRecurrent : public BoundaryTime {

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

#endif
