/*!
\file fconfig.hh
\brief Configuration for Fields and Coordinates classes
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_FCONFIG_HH
#define SPECTRUM_FCONFIG_HH

#include "../compiletime_lists.hh"

namespace Spectrum {

template <
      SpecieId specieid = SpecieId::proton_core,
      CoordinateSystem Pos_sys_ = CoordinateSystem::cartesian,
      CoordinateSystem Mom_sys_ = CoordinateSystem::cartesian
>
class FConfig {
   static constexpr auto specie = Specie<specieid>();
   static constexpr auto Pos_sys = Pos_sys_;
   static constexpr auto Mom_sys = Mom_sys_;
   static constexpr auto Vel_sys = Mom_sys_;
   // A CoordinateSystem is 'radial' if the 0th component is the radial component.
   static constexpr bool Pos_radial = (Pos_sys == CoordinateSystem::polar || Pos_sys == CoordinateSystem::spherical);
   static constexpr bool Mom_radial = (Mom_sys == CoordinateSystem::polar || Mom_sys == CoordinateSystem::spherical || Mom_sys == CoordinateSystem::pitchangle);
   static  constexpr bool Vel_radial = Mom_radial;
   static constexpr bool Mag_radial = false;
   static constexpr bool Flum_radial = false;
   static constexpr bool Fluv_radial = false;
};

}

#endif
