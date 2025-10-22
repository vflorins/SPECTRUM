/*!
\file distribution.config.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM distribution class
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/


#ifndef SPECTRUM_DISTRIBUTION_CONFIG_HH
#define SPECTRUM_DISTRIBUTION_CONFIG_HH

#include "common/compiletime_lists.hh"

namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM distribution class
\author Lucius Schoenbaum
\date 09/05/2025
*/
template<

   // todo
   int a

>
struct DistributionConfig {

};


using DistributionDefault = DistributionConfig<1, 2, 3>;

}


#endif