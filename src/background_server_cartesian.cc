/*!
\file background_server_cartesian.cc
\brief Implements a background class using data from a uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_server_cartesian.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServerCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename Fields>
BackgroundServerCartesian<Fields>::BackgroundServerCartesian(void)
                         : BackgroundServer(bg_name_server_cartesian, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
template <typename Fields>
BackgroundServerCartesian<Fields>::BackgroundServerCartesian(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                         : BackgroundServer(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundServerCartesian<Fields>::BackgroundServerCartesian(const BackgroundServerCartesian& other)
                         : BackgroundServer(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

};
