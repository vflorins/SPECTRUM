/*!
\file initial_base.hh
\brief Declares a base class to specify initial (starting) conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INITIAL_BASE_HH
#define SPECTRUM_INITIAL_BASE_HH

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "config.h"
#include "common/params.hh"
#include "common/physics.hh"
#include "common/random.hh"
#include <memory>

namespace Spectrum {

//! Condition is in time
const uint16_t INITIAL_TIME = 0x0010;

//! Condition is in space
const uint16_t INITIAL_SPACE = 0x0020;

//! Condition is in momentum
const uint16_t INITIAL_MOMENTUM = 0x0040;

//! Condition is a single point
const uint16_t INITIAL_POINT = 0x0080;

//! Condition is a curve
const uint16_t INITIAL_CURVE = 0x0100;

//! Condition is a surface
const uint16_t INITIAL_SURFACE = 0x0200;

//! Condition is a volume
const uint16_t INITIAL_VOLUME = 0x0400;

//! Clone function pattern
#define CloneFunctionInitial(T) std::unique_ptr<InitialBase> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Base class to produce a distribution of starting points
\author Vladimir Florinski

The InitialXXXXX" classes describe initial phase space coordinates for trajectories. Mixed position-momentum conditions are not supported at this time. The controlling program should have two objects, one for position, and the other for momentum.

Parameters:
*/
class InitialBase : public Params {

protected:

//! Preferred direction (transient). We don't want a persistent argument because the direction could change if we have variable initial position.
   GeoVector axis;

//! Default constructor (protected, class not designed to be instantiated)
   InitialBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   InitialBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   InitialBase(const InitialBase& other);

//! Set up the initial condition generator based on "params"
   virtual void SetupInitial(bool construct);

//! Compute the internal position or momentum
   virtual void EvaluateInitial(void);

public:

//! Destructor
   virtual ~InitialBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<InitialBase> Clone(void) const = 0;

//! Set up the object's persistent class data members (generic)
   void SetupObject(const DataContainer& cont_in);

//! Produce a realization of initial conditions based on an internal distribution function
   GeoVector GetSample(const GeoVector& axis_in);

//! Tell if the class is for time
   bool IsInitialTime(void) const;

//! Tell if the class is for space coordinate
   bool IsInitialSpace(void) const;

//! Tell if the class is for momentum coordinate
   bool IsInitialMomentum(void) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTime class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialTime class
const std::string init_name_time = "InitialTime";

/*!
\brief Initial time condition. This class doubles as the "final" time for backward simulations. This class will not be used in the same way as position or momentum initial conditions, but is included here for completion.
\author Juan G Alonso Guzman

Parameters: (InitialBase)
*/
class InitialTime : public InitialBase {

public:

//! Default constructor
   InitialTime(void);

//! Copy constructor
   InitialTime(const InitialTime& other);

//! Clone function
   CloneFunctionInitial(InitialTime);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTable class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Starting points from a table
\author Juan G Alonso Guzman

Parameters: (InitialBase), std::string init_file_name, double scale, bool random
*/
class InitialTable : public InitialBase {

protected:

//! Flag to iterate through initial positions randomly (true) or in sequence (false)
   bool random;

//! Array with initial positions
   std::vector <GeoVector> initvec;

//! Table entry counter
   int table_counter;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

public:

//! Default constructor
   InitialTable(void);

//! Constructor with arguments (to speed up construction of derived classes)
   InitialTable(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   InitialTable(const InitialTable& other);

//! Clone function
   CloneFunctionInitial(InitialTable);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialBase inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/01/2024
\return True if the class describes initial condition in time
*/
inline bool InitialBase::IsInitialTime(void) const
{
   return BITS_RAISED(_status, INITIAL_TIME);
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\return True if the class describes initial condition in space
*/
inline bool InitialBase::IsInitialSpace(void) const
{
   return BITS_RAISED(_status, INITIAL_SPACE);
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\return True if the class describes initial condition in momentum
*/
inline bool InitialBase::IsInitialMomentum(void) const
{
   return BITS_RAISED(_status, INITIAL_MOMENTUM);
};

};

#endif
