/*!
\file initial_base.cc
\brief Implements a base class to specify initial (starting) conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_base.hh"
#include <iostream>
#include <fstream>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/27/2021
*/
InitialBase::InitialBase(void)
           : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
InitialBase::InitialBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
           : Params(name_in, specie_in, status_in)
{
};

/*!
\brief Copy constructor (protected, class not designed to be instantiated)
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialBase::InitialBase(const InitialBase& other)
           : Params(other)
{
// Params' constructor sets the state to "STATE_NONE"
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupInitial()" method.
*/
void InitialBase::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupInitial(false);
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialBase::SetupInitial(bool construct)
{
// Only needed in the parent version
   container.Reset();
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 12/27/2021
\note This is only a stub
*/
void InitialBase::EvaluateInitial(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 03/25/2021
\param[in] axis_in Preferred direction
\note This is a common routine that the derived classes should not change.
*/
GeoVector InitialBase::GetSample(const GeoVector& axis_in)
{
// The only possible error is if "axis" is zero, so we don't throw an exception here.
   axis = UnitVec(axis_in);
   EvaluateInitial();

// Return the internal position or momentum. Initial time will never be assigned this way.
   if(BITS_RAISED(_status, INITIAL_SPACE)) return _pos;
   else return _mom;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTime methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/01/2024
*/
InitialTime::InitialTime(void)
           : InitialBase(init_name_time, 0, INITIAL_TIME)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialTime::InitialTime(const InitialTime& other)
           : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_TIME);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
*/
InitialTable::InitialTable(void)
            : InitialBase("", 0, INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
InitialTable::InitialTable(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
            : InitialBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialTable::InitialTable(const InitialTable& other)
            : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_POINT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialTable::SetupInitial(bool construct)
{
   std::string init_file_name, coord_type;
   std::ifstream init_file;
   int i, size;
   double scale;
   GeoVector entry;

// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);
   container.Read(&init_file_name);
   container.Read(&scale);
   container.Read(&random);

// Input initial positions from file
   init_file.open(init_file_name.c_str());

   init_file >> coord_type;
   init_file >> size;
#ifdef GEO_DEBUG
   if(coord_type == "RTP") {
      std::cerr << "Reading initial positions file in spherical coordinates." << std::endl;
   }
   else if(coord_type == "XYZ") {
      std::cerr << "Reading initial positions file in cartesian coordinates." << std::endl;
   }
   else {
      std::cerr << "Reading initial positions file in unrecognized coordinate type. "
                << "Defaulting to cartesian coordinates." << std::endl;
   };
#endif

   initvec.resize(size);
   for(i = 0; i < size; i++) {
      init_file >> entry[0];
      init_file >> entry[1];
      init_file >> entry[2];
      if(coord_type == "RTP") {
         entry[0] *= scale;
         entry.RTP_XYZ();
      }
      else entry = entry * scale;
      initvec[i] = entry;
   };
   init_file.close();

// Initialize table counter
   table_counter = 0;
};

};
