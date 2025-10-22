/*!
\file class.cc
\brief Implements a simple class for entering parameters
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "params.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Class methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/08/2025
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
Params::Params(const std::string_view& name_in, uint16_t status_in)
      : class_name(name_in),
        _status(status_in)
{
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/08/2025
\param[in] other Object to initialize from
\note Performs a deep copy of persistent variables only
*/
Params::Params(const Params& other)
      : class_name(other.class_name),
        rng(other.rng),
        container(other.container),
        _status(STATE_NONE)
{
};



};
