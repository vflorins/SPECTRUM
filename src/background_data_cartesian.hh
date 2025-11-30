/*!
\file background_data_cartesian.hh
\brief Declares a background class using data from uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_CARTESIAN_HH
#define SPECTRUM_BACKGROUND_DATA_CARTESIAN_HH

#include "server_interface.hh"
#include "common/status.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDataCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background from stored data as a uniform Cartesian grid
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

*/
template <typename HConfig_>
class BackgroundDataCartesian : public ServerInterface<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDataCartesian";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using BlockCache = BlockCache<HConfig>;
   using MPI = MPI<HConfig>;
   using Block = Block<HConfig>;
   using BlockPtr = std::shared_ptr<Block>;

   using DataFields = Config::DataFields;

   using ServerInterface = ServerInterface<HConfig>;
   using ServerInterface::_inquiry;
   using ServerInterface::stencil;
   using ServerInterface::MPIInquiryType;
   using ServerInterface::MPIStencilType;
   using ServerInterface::MPIBlockType;

   static_assert(!(!HConfig::MPI_enabled() && !Config::servers_are_workers), "Servers and worker duties cannot be divided unless the simulation is parallel.");

   static constexpr int num_ghost_cells = Config::num_ghost_cells;

   static constexpr int server_interp_order = Config::server_interp_order;

   static constexpr double dmax0 = Config::dmax0;

   static_assert(server_interp_order <= 1, "Interpolation orders > 1 are not supported.");

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

//! Stencil block pointer
   BlockPtr block_stn;

//! Load interpolation stencil using interior zones
   void InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi, const GeoVector offset_lo, const GeoVector offset_hi, const GeoVector delta);

//! Generate an interpolation stencil in 3D
   int BuildInterpolationStencil(const GeoVector& pos);

//! Find block order in the cache or get a block from the server if not in the cache
   int RequestBlock(void);

//! Get variables directly from data reader
   template <typename Fields, typename RequestedFields>
   status_t Evaluate_FromReader(Fields& fields);

//! Get variables using 0th order interpolation
   template <typename Fields, typename RequestedFields>
   status_t Evaluate_Interp0(const GeoVector& pos, Fields& fields);

//! Get variables using 1st order interpolation
   template <typename Fields, typename RequestedFields>
   status_t Evaluate_Interp1(const GeoVector& pos, Fields& fields);

//! Get gradients using 1st order interpolation
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t Gradients_Interp1(Coordinates& coords, Fields& fields);

public:

//! Default constructor
   BackgroundDataCartesian(void) = default;

//! Destructor
   ~BackgroundDataCartesian() override = default;

   /*
    * TODO - the change still needed is
    *  to provide all backgrounds with Start() and Finish() for init and deinit.
    *  These are only non-trivial in the BackgroundData cases,
    *  where they (of course) start up and close up MPI comms.
    *
    * TODO
    *  Review this class to see what can be removed from public
    *  to establish a completely consistent background class API:
    *  Start
    *  EvaluateBackground
    *  EvaluateBackroundDerivatives
    *  Finish
    *
    */

protected:

//! Obtain the variables
   template <typename Fields>
   void GetVariables(double t, const GeoVector& pos, Fields& fields, double& dmax);

////! Obtain the gradients
//   template <typename Fields>
//   void GetGradients(Fields& fields, DerivativeData& ddata);


public: // data background API:

//! Front end set up prior to main loop
   void Start(void);

//! Front end clean up tasks after the main loop
   void Finish(void);

public: // general background API:

//! Compute the maximum distance per time step
   template <typename Coordinates>
   status_t EvaluateDmax(Coordinates&, double&);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

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
#include "background_data_cartesian.cc"

#endif
