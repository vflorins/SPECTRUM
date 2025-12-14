/*!
\file background_data_base.cc
\brief Declares a base background class using grid data on distributed memory
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_BASE_HH
#define SPECTRUM_BACKGROUND_DATA_BASE_HH

#include "server_interface.hh"
#include "cache_lru.hh"
#include "common/vectors.hh"
#include "common/status.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDataBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background from stored data as a uniform Cartesian grid
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
*/
template <typename HConfig_>
class BackgroundDataBase : public ServerInterface<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDataBase";

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using BlockCache = BlockCache<HConfig>;
   using MPI = MPI<HConfig, HConfig::MPI_enabled()>;
   using Block = Block<HConfig>;
   using BlockPtr = std::shared_ptr<Block>;

   using DataFields = Config::DataFields;

   using ServerInterface = ServerInterface<HConfig>;
   using ServerInterface::_inquiry;
   using ServerInterface::stencil;
   using ServerInterface::MPIInquiryType;
   using ServerInterface::MPIStencilType;
   using ServerInterface::MPIBlockType;

   static constexpr bool allow_server_worker = Config::allow_server_worker;

   static constexpr int num_ghost_cells = Config::num_ghost_cells;

   static constexpr int server_interpolation_order = Config::server_interpolation_order;

   static constexpr double dmax0 = Config::dmax0;

   static_assert(server_interpolation_order <= 1, "Interpolation orders > 1 are not supported.");
   static_assert(!(!HConfig::MPI_enabled() && !allow_server_worker), "Servers and worker duties cannot be divided unless the simulation is parallel.");

protected:

//! locally cached blocks, with LRU policy
   BlockCache cache_line;

//! Smallest position in the domain
   GeoVector domain_min;

//! Largest position in the domain
   GeoVector domain_max;

//! Status of the most recently computed stencil
   int stencil_status = 0;

//! Counts of different stencil outcomes
   int stencil_outcomes[3];

//! Count of total blocks requested
   int num_blocks_requested;

//! Primary block pointer
   BlockPtr block_pri;

//! Secondary block pointer
   BlockPtr block_sec;

   // todo review, I think this should be a temporary in methods and not here
//! Stencil block pointer
   BlockPtr block_stn;

public:

//! Default constructor
   BackgroundDataBase(void) = default;

//! Destructor
   ~BackgroundDataBase() override = default;

public:

//! Load interpolation stencil using interior zones
   void InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi, const GeoVector offset_lo, const GeoVector offset_hi, const GeoVector delta);

//! Find block order in the cache or get a block from the server if not in the cache
   int RequestBlock(void);

//! Generate an interpolation stencil in 3D
   virtual int BuildInterpolationStencil(const GeoVector& pos) = 0;

public: // data background API:

//! Front end set up prior to main loop
   void Start(void);

//! Front end clean up tasks after the main loop
   void Finish(void);

public: // debug:

   //! Print how many times internal/external interpolators were used
   void PrintStencilOutcomes(void) const requires (HConfig::buildmode == BuildMode::debug);

//! Print how many blocks were requested
   void PrintNumBlocksRequested(void) const requires (HConfig::buildmode == BuildMode::debug);

   int GetNCachedBlocks(void) const requires (HConfig::buildmode == BuildMode::debug);

   void InvalidateCache(void) requires (HConfig::buildmode == BuildMode::debug);

};


};

// Something like this is needed for templated classes
#include "background_data_base.cc"

#endif