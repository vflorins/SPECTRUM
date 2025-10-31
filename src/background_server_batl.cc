/*!
\file background_server_batl.cc
\brief Implements a background class using data from BATL adaptive mesh on distributed memory
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_server_batl.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServerBATL methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/05/2022
*/
BackgroundServerBATL::BackgroundServerBATL(void)
                    : BackgroundServerCartesian(bg_name_server_batl, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 10/05/2022
\param[in] other Object to initialize from

A copy constructor should first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundServerBATL::BackgroundServerBATL(const BackgroundServerBATL& other)
                    : BackgroundServerCartesian(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

};
