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


// todo: deprecated
///*!
//\brief Multi-purpose class storing physical data defined at a spatial location.
//Used in the data interface with the (pseudo-)particle tracer, and to house
//data distributed on the grid.
//\author Lucius Schoenbaum
//\author Vladimir Florinski
//\date 08/26/2025
//*/
//template <Field::Id nameid, typename ... Ts>
//class Species: public Fields<Ts...> {
//
//public:
//
//   static constexpr const std::string_view name = Field::Names[nameid];
//
//};



};

#endif
