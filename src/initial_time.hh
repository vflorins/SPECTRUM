/*!
\file initial_space.hh
\brief Declares several classes to specify time initial conditions
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INITIAL_TIME_HH
#define SPECTRUM_INITIAL_TIME_HH

#include "initial_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeFixed class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialTimeFixed class
const std::string init_name_time_fixed = "InitialTimeFixed";

/*!
\brief Starting points at a fixed time
\author Juan G Alonso Guzman

Parameters: (InitialBase), double inittime
*/
template <typename Trajectory_>
class InitialTimeFixed : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_t;

protected:

//! Time (persistent)
   double inittime;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialTimeFixed(void);

//! Copy constructor
   InitialTimeFixed(const InitialTimeFixed& other);

//! Clone function
   CloneFunctionInitial(InitialTimeFixed);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeInterval class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialTimeInterval class
const std::string init_name_time_interval = "InitialTimeInterval";

/*!
\brief Uniformly distributed starting times on an interval
\author Juan G Alonso Guzman

Parameters: (InitialBase), double starttime, double endtime, int n_intervals
*/
template <typename Trajectory_>
class InitialTimeInterval : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_t;
   using InitialBase::rng;

protected:

//! Random or evenly distributed
   bool randomtime;

//! First point (persistent)
   double starttime;

//! Second point (persistent)
   double endtime;

//! Increment (persistent)
   double increment;

//! Constructor with arguments (to speed up construction of derived classes)
   InitialTimeInterval(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialTimeInterval(void);

//! Copy constructor
   InitialTimeInterval(const InitialTimeInterval& other);

//! Clone function
   CloneFunctionInitial(InitialTimeInterval);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeTable class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialTimeTable class
const std::string init_name_time_table = "InitialTimeTable";

/*!
\brief Starting times from a table
\author Juan G Alonso Guzman

Parameters: (InitialTable)
*/
template <typename Trajectory_>
class InitialTimeTable : public InitialTable<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;
   using InitialTable = InitialTable<Trajectory, double>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_t;
   using InitialBase::rng;

   using InitialTable::random;
   using InitialTable::table_counter;
   using InitialTable::initquant;
   using InitialTable::SetupInitial;

protected:

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialTimeTable(void);

//! Copy constructor
   InitialTimeTable(const InitialTimeTable& other);

//! Clone function
   CloneFunctionInitial(InitialTimeTable);
};

};

// Something like this is needed for templated classes
#include "initial_time.cc"

#endif
