/*!
\file initial_space.cc
\brief Implements several classes to specify time initial conditions
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_time.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeFixed methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename HConfig>
InitialTimeFixed<HConfig>::InitialTimeFixed(void)
                : InitialBase(init_name, INITIAL_TIME | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialTimeFixed<HConfig>::InitialTimeFixed(const InitialTimeFixed& other)
                : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialTimeFixed<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(inittime);

// Pre-assign "_coords.Time()" so that it never needs to change
   _coords.Time() = inittime;
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename HConfig>
void InitialTimeFixed<HConfig>::EvaluateInitial(void)
{
// Nothing to do - the value of "_coords.Time()" was assigned in "Setupinitial()"
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeInterval methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename HConfig>
InitialTimeInterval<HConfig>::InitialTimeInterval(void)
                   : InitialBase(init_name, INITIAL_TIME | INITIAL_CURVE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
InitialTimeInterval<HConfig>::InitialTimeInterval(const std::string& name_in, uint16_t status_in)
                   : InitialBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialTimeInterval<HConfig>::InitialTimeInterval(const InitialTimeInterval& other)
                   : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   RAISE_BITS(_status, INITIAL_CURVE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialTimeInterval<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);

   int n_intervals;
   container.Read(starttime);
   container.Read(endtime);
   container.Read(n_intervals);

   if (n_intervals <= 0) randomtime = true;
   else {
      randomtime = false;
      increment = (endtime - starttime) / (double)n_intervals;
      _coords.Time() = starttime - 0.5 * increment;
   };
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename HConfig>
void InitialTimeInterval<HConfig>::EvaluateInitial(void)
{
   if (randomtime) _coords.Time() = starttime + (endtime - starttime) * rng->GetUniform();
   else _coords.Time() += increment;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename HConfig>
InitialTimeTable<HConfig>::InitialTimeTable(void)
                : InitialTable(init_name, INITIAL_TIME | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialTimeTable<HConfig>::InitialTimeTable(const InitialTimeTable& other)
                : InitialTable(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
template <typename HConfig>
void InitialTimeTable<HConfig>::EvaluateInitial(void)
{
   if (random) {
// Generate random integer between 0 and initquant.size() - 1
      table_counter = rng->GetUniform() * initquant.size();
// Pull time in randomly selected place on the table
      _coords.Time() = initquant[table_counter];
   }
   else {
// Pull next time on the table
      _coords.Time() = initquant[table_counter++];
// If all positions have been sampled, reset the counter
      if (table_counter == initquant.size()) table_counter = 0;
   };
};

};
