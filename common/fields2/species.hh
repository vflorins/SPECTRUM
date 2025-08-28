/*!
\file namedmhdtuple.hh
\author Vladimir Florinski
\author Lucius Schoenbaum
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPECIES_HH
#define SPECTRUM_SPECIES_HH

#include <iostream>
#include <any>
#include <stdexcept>
#include "generated/field_types.hh"

namespace Spectrum {

template <typename ... Ts>
class Fields;

/*!
\brief Multi-purpose class storing physical data defined at a spatial location.
Used in the data interface with the (pseudo-)particle tracer, and to house
data distributed on the grid.
\author Lucius Schoenbaum
\author Vladimir Florinski
\date 08/26/2025
\note Recursion in protected methods takes place only once, at compile time.
*/
template <Field::Id nameid, typename ... Ts>
class Species: public Fields<Ts...> {

public:

   static constexpr const std::string_view name = Field::Names[nameid];

};

/*!
\author Lucius Schoenbaum
\date 08/26/2025
A Fields type can hold a Species type, but this is perhaps unintuitive, so we rename Fields for this use case.
A Fields type can also hold a mixture of Species and Field, and this can be useful in some rarer cases.
For example, if a SpeciesTuple needs to hold an indicator variable, then that is possible.
*/
template<typename ... Ts>
using SpeciesTuple = Fields<Ts...>;

};

#endif
