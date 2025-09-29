/*!
\file fluid_fields.hh
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPECIES_HH
#define SPECTRUM_SPECIES_HH

#include "fluid_specie.hh"

namespace Spectrum {

/*!
\brief Provides a container combining all conserved fluid variables into single contiguous storage.
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/26/2025
\note A Fields type can hold a FluidSpecie type, but this is perhaps unintuitive, so we rename Fields for this use case.
Note that a general Fields type can hold a mixture of FluidSpecie and Field.
For example, an indicator variable could be defined that applies to the set of all target species as a group.
For this reason, we do not name this type 'FluidSpecies' but instead choose the name 'FluidFields'.
*/
template<typename ... Ts>
using FluidFields = Fields<Ts...>;

}

#endif
