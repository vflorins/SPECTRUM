/*!
\file rsconfig.hh
\brief Configuration input for Riemann Solver classes
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_RSCONFIG_HH
#define SPECTRUM_RSCONFIG_HH

#include "../common/compiletime_lists.hh"

namespace Spectrum {

template <
   Model model_,
   SpecieId specieid,
   Passivity passivity_,
//! Enables iteration to calculate maximum extent of the Riemann fan
      bool iterate_HLL_wave_
>
class RiemannSolverConfig {

   static constexpr auto model = model_;
   static constexpr auto specie = Specie<specieid>();
   static constexpr auto passivity = passivity_;
   static constexpr auto iterate_HLL_wave = iterate_HLL_wave_;

};


}

#endif
