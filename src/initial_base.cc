/*!
\file initial_base.cc
\brief Implements a base class to specify initial (starting) conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_base.hh"
#include <iostream> // std::cout
#include <fstream> // std::ifsteam

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
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
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
\author Juan G Alonso Guzman
\date 08/10/2024
\return Initial time from internal distribution
\note This is a common routine that the derived classes should not change.
*/
double InitialBase::GetTimeSample(void)
{
// Evaluate internal distribution.
   EvaluateInitial();

// Return the internal time.
   return _t;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/10/2024
\return Initial position from internal distribution
\note This is a common routine that the derived classes should not change.
*/
GeoVector InitialBase::GetPosSample(void)
{
// Evaluate internal distribution.
   EvaluateInitial();

// Return the internal position.
   return _pos;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/10/2024
\param[in] axis_in Preferred direction
\return Initial momentum from internal distribution
\note This is a common routine that the derived classes should not change.
*/
GeoVector InitialBase::GetMomSample(const GeoVector& axis_in)
{
// Evaluate internal distribution. The only possible error is if "axis" is zero, so we don't throw an exception here.
   axis = UnitVec(axis_in);
   EvaluateInitial();

// Return the internal momentum.
   return _mom;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
*/
template <class tableClass>
InitialTable<tableClass>::InitialTable(void)
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
template <class tableClass>
InitialTable<tableClass>::InitialTable(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                        : InitialBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <class tableClass>
InitialTable<tableClass>::InitialTable(const InitialTable& other)
                        : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <class tableClass>
void InitialTable<tableClass>::SetupInitial(bool construct)
{
   std::string init_file_name, coord_type;
   std::ifstream init_file;
   int i, size;
   double scale;
   tableClass entry;

// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(init_file_name);
   container.Read(scale);
   container.Read(random);

// Input initial quantities from file
   init_file.open(init_file_name.c_str());

   init_file >> coord_type;
   init_file >> size;

#ifdef GEO_DEBUG

   if (coord_type == "S") {
      std::cerr << "Reading initial times from file." << std::endl;
   }
   else if (coord_type == "RTP") {
      std::cerr << "Reading initial vectors file in spherical coordinates." << std::endl;
   }
   else if (coord_type == "XYZ") {
      std::cerr << "Reading initial vectors file in cartesian coordinates." << std::endl;
   }
   else {
      std::cerr << "Reading file of unrecognized initial quantities. "
                << "Defaulting to scalars." << std::endl;
   };

#endif

   initquant.resize(size);
   for (i = 0; i < size; i++) {
      init_file >> entry;
// The redundant conversion to GeoVector type is necessary due to the templated nature of this class.
      if (coord_type == "RTP") GeoVector(entry).RTP_XYZ();
      initquant[i] = entry * scale;
   };
   init_file.close();

// Initialize table counter
   table_counter = 0;
};

};
