/*!
\file server_config.hh
\brief Data server configurator
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SERVER_CONFIG_HH
#define SPECTRUM_SERVER_CONFIG_HH

#include "config.h"

#if SERVER_TYPE == SERVER_SELF
// This needs to be here for proper compilation (e.g. for ExServerError to be defined)
#include "server_base.hh"
#elif SERVER_TYPE == SERVER_CARTESIAN
#include "server_cartesian.hh"
#elif SERVER_TYPE == SERVER_BATL
#include "server_batl.hh"
#else
#error Unsupported Server type
#endif

#endif
