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
InitialTimeFixed::InitialTimeFixed(void)
                : InitialBase(init_name_time_fixed, 0, INITIAL_TIME | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialTimeFixed::InitialTimeFixed(const InitialTimeFixed& other)
                : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   RAISE_BITS(_status, INITIAL_POINT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialTimeFixed::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);
   container.Read(inittime);

// Pre-assign "_t" so that it never needs to change
   _t = inittime;
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
void InitialTimeFixed::EvaluateInitial(void)
{
// Nothing to do - the value of "_t" was assigned in "Setupinitial()"
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeInterval methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
InitialTimeInterval::InitialTimeInterval(void)
                   : InitialBase(init_name_time_interval, 0, INITIAL_TIME | INITIAL_CURVE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
InitialTimeInterval::InitialTimeInterval(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                   : InitialBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialTimeInterval::InitialTimeInterval(const InitialTimeInterval& other)
                   : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   RAISE_BITS(_status, INITIAL_CURVE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialTimeInterval::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);

   int n_intervals;
   container.Read(starttime);
   container.Read(endtime);
   container.Read(n_intervals);

   if(n_intervals <= 0) randomtime = true;
   else {
      randomtime = false;
      increment = (endtime - starttime) / (double)n_intervals;
      _t = starttime - 0.5 * increment;
   };
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
void InitialTimeInterval::EvaluateInitial(void)
{
   if(randomtime) _t = starttime + (endtime - starttime) * rng->GetUniform();
   else _t += increment;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTimeTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
InitialTimeTable::InitialTimeTable(void)
                : InitialTable(init_name_time_table, 0, INITIAL_TIME | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialTimeTable::InitialTimeTable(const InitialTimeTable& other)
                : InitialTable(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   RAISE_BITS(_status, INITIAL_POINT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/08/2024
*/
void InitialTimeTable::EvaluateInitial(void)
{
   if(random) {
// Generate random integer between 0 and initquant.size() - 1
      table_counter = rng->GetUniform() * initquant.size();
// Pull time in randomly selected place on the table
      _t = initquant[table_counter];
   }
   else {
// Pull next time on the table
      _t = initquant[table_counter++];
// If all positions have been sampled, reset the counter
      if(table_counter == initquant.size()) table_counter = 0;
   };
};

};
