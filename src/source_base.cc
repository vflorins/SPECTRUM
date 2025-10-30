/*!
\file diffusion_base.cc
\brief Implements a base class to compute source terms
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "source_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
*/
SourceBase::SourceBase(void)
          : Params("", 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
SourceBase::SourceBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
          : Params(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupSource()" with the argument of "true".
*/
SourceBase::SourceBase(const SourceBase& other)
          : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupSource(true);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupSource()" method.
*/
void SourceBase::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupSource(false);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void SourceBase::SetupSource(bool construct)
{
// Only needed in the parent version
   container.Reset();
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\note This is only a stub.
*/
void SourceBase::EvaluateSource(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] t_in      Time
\param[in] pos_in    Position
\param[in] mom_in    Momentum (p,mu,phi) coordinates
\param[in] spdata_in Spatial data at the required location
\return Value of source term
\note This is a common routine that the derived classes should not change.
*/
double SourceBase::GetSource(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const SpatialData& spdata_in, double dt_in)
{
   SetState(t_in, pos_in, mom_in);
   _spdata._mask = spdata_in._mask;
   _spdata = spdata_in;
   _dt = dt_in;

   EvaluateSource();
   return SourceTerm;
};

};
