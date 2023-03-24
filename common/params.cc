/*!
\file params.cc
\brief Implements a simple class for entering parameters
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "params.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Params methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/02/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
Params::Params(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
      : class_name(name_in),
        specie(specie_in),
        _status(status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/01/2021
\param[in] other Object to initialize from
\note Performs a deep copy of persistent variables only
*/
Params::Params(const Params& other)
      : class_name(other.class_name),
        specie(other.specie),
        rng(other.rng),
        container(other.container),
        _status(STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 11/12/2020
\param[in] t_in   Time
\param[in] pos_in Position
\param[in] mom_in Momentum (optional)
*/
void Params::SetState(double t_in, const GeoVector& pos_in, const GeoVector& mom_in)
{
   _t = t_in;
   _pos = pos_in;
   _mom = mom_in;
};

/*!
\author Vladimir Florinski
\date 11/29/2020
\param[out] t_out   Time
\param[out] pos_out Position
\param[out] mom_out Momentum
*/
void Params::GetState(double& t_out, GeoVector& pos_out, GeoVector& mom_out) const
{
   t_out = _t;
   pos_out = _pos;
   mom_out = _mom;
};

};
