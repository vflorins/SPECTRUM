/*!
\file background_data_base.hh
\brief Declares a background class using data from a grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_BASE_HH
#define SPECTRUM_BACKGROUND_DATA_BASE_HH

#include "server_interface_base__old.hh"
#include "cache_lru.hh"
#include "common/mpi_config.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDataBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

///*!
//\brief Common methods for a plasma background using data from a grid on distributed memory
//\author Vladimir Florinski
//\author Juan G Alonso Guzman
//\author Lucius Schoenbaum
//*/
//template <typename HConfig_>
//class BackgroundDataBase : public InterfaceCartesian<HConfig_> {
//public:
//
//   //! Readable name of the class
//   static constexpr std::string_view name = "BackgroundDataBase";
//
//public:
//
//   using HConfig = HConfig_;
//   using Config = HConfig::BackgroundConfig;
//   using BlockCache = BlockCache<HConfig>;
//   using InterfaceBase = InterfaceBase<HConfig>;
//   using MPI = MPI<HConfig>;
//
//protected:
//
////! Cache line
//   static BlockCache cache_line = BlockCache();
//
//public:
//
//////! Front end set up prior to main loop
////   static void Start(void);
////
//////! Front end clean up tasks after the main loop
////   static void Finish(void);
////
//////! Frontend tasks during the main loop
////   static int ServerFunctions(void);
////
//////! Return the number of cached blocks
////   static int GetNCachedBlocks(void);
////
//////! Empty the cache
////   static void InvalidateCache(void);
////
////   static void Start(void);
////
//////! Signal the backend this client no longer needs its service
////   static void Stop(void);
//
//};


};

// Something like this is needed for templated classes
#include "background_data_base__old.cc"

#endif
