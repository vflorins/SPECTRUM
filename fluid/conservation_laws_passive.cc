/*!
\file conservation_laws_passive.cc
\brief Implements rules for a passively advected scalar (PAS)
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "fluid/conservation_laws_passive.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Passive::PrimitiveState methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] pas_in PAS
*/
SPECTRUM_DEVICE_FUNC Passive::PrimitiveState::PrimitiveState(double pas_in)
{
   pas() = pas_in;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC Passive::PrimitiveState::PrimitiveState(const SimpleArray<double, CL_PASSIVE_NVARS>& other)
                                            : SimpleArray<double, CL_PASSIVE_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Vector of conserved variables
*/
template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC Passive::ConservedState Passive::PrimitiveState::ToConserved(const ActiveConsLaw::PrimitiveState& active_prim) const
{
   return Passive::ConservedState(active_prim.den() * pas());
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Flux vector
*/
template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC Passive::FluxFunction Passive::PrimitiveState::ToFlux(const ActiveConsLaw::PrimitiveState& active_prim) const
{
   return Passive::FluxFunction(active_prim.den() * pas() * active_prim.vel()[0]);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Passive::ConservedState methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] dpas_in primary density*PAS
*/
SPECTRUM_DEVICE_FUNC Passive::ConservedState::ConservedState(double dpas_in)
{
   dpas() = dpas_in;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC Passive::ConservedState::ConservedState(const SimpleArray<double, CL_PASSIVE_NVARS>& other)
                                            : SimpleArray<double, CL_PASSIVE_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Vector of primitive variables
*/
template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC Passive::PrimitiveState Passive::ConservedState::ToPrimitive(const ActiveConsLaw::PrimitiveState& active_prim) const
{
   return Passive::PrimitiveState(dpas() / active_prim.den());
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Flux vector
*/
template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC Passive::FluxFunction Passive::ConservedState::ToFlux(const ActiveConsLaw::PrimitiveState& active_prim) const
{
   return Passive::FluxFunction(dpas() * active_prim.vel()[0]);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Passive::FluxFunction methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] dvarf_in Flux of primary density*PAS
*/
SPECTRUM_DEVICE_FUNC Passive::FluxFunction::FluxFunction(double dpasf_in)
{
   dpasf() = dpasf_in;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC Passive::FluxFunction::FluxFunction(const SimpleArray<double, CL_PASSIVE_NVARS>& other)
                                          : SimpleArray<double, CL_PASSIVE_NVARS>(other)
{
};

};
