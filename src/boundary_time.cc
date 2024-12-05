/*!
\file boundary_time.cc
\brief Implements several classes representing temporal boundaries
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "boundary_time.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTime methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/26/2021
*/
BoundaryTime::BoundaryTime(void)
            : BoundaryBase("", 0, BOUNDARY_TIME)
{
};

/*!
\author Vladimir Florinski
\date 01/25/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryTime::BoundaryTime(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
            : BoundaryBase(name_in, 0, status_in)
{
};

/*!
\author Vladimir Florinski
\date 01/20/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBoundary()" with the argument of "true".
*/
BoundaryTime::BoundaryTime(const BoundaryTime& other)
            : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/14/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryTime::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(timemark);
};

/*!
\author Vladimir Florinski
\date 01/26/2021
*/
void BoundaryTime::EvaluateBoundary(void)
{
   _delta = _t - timemark;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTimeExpire methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/26/2021
*/
BoundaryTimeExpire::BoundaryTimeExpire(void)
                  : BoundaryTime(bnd_name_time_expire, 0, BOUNDARY_TIME | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/20/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBoundary()" with the argument of "true".
*/
BoundaryTimeExpire::BoundaryTimeExpire(const BoundaryTimeExpire& other)
                  : BoundaryTime(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/14/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryTimeExpire::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryTime::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTimePass methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/27/2021
*/
BoundaryTimePass::BoundaryTimePass(void)
                : BoundaryTime(bnd_name_time_pass, 0, BOUNDARY_TIME)
{
   max_crossings = -1;
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBoundary()" with the argument of "true".
*/
BoundaryTimePass::BoundaryTimePass(const BoundaryTimePass& other)
                : BoundaryTime(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = -1;
};

/*!
\author Vladimir Florinski
\date 01/14/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryTimePass::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryTime::SetupBoundary(false);
   max_crossings = -1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryTimeRecurrent methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/12/2022
*/
BoundaryTimeRecurrent::BoundaryTimeRecurrent(void)
                     : BoundaryTime(bnd_name_time_pass, 0, BOUNDARY_TIME | BOUNDARY_RECURRENT)
{
};

/*!
\author Vladimir Florinski
\date 05/12/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBoundary()" with the argument of "true".
*/
BoundaryTimeRecurrent::BoundaryTimeRecurrent(const BoundaryTimeRecurrent& other)
                     : BoundaryTime(other)
{
   RAISE_BITS(_status, BOUNDARY_RECURRENT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 05/12/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryTimeRecurrent::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryTime::SetupBoundary(false);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void BoundaryTimeRecurrent::EvaluateBoundary(void)
{
   _delta = _t - (CrossingsMade() + 1) * timemark;
};

};
