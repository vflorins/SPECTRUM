/*!
\file background_data_cartesian.hh
\brief Declares a background class using data from uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_CARTESIAN_HH
#define SPECTRUM_BACKGROUND_DATA_CARTESIAN_HH

#include "background_data_base__old.hh"
#include "server_interface.hh"

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

   using ServerInterface = ServerInterface<HConfig>;
   using ServerInterface::_inquiry;
   using ServerInterface::stencil;

   static_assert(!(!HConfig::MPI_enabled() && !Config::servers_are_workers), "Servers and worker duties cannot be divided unless the simulation is parallel.");

   static constexpr int num_ghost_cells = Config::num_ghost_cells;

protected:

   void LoadFromReader(BlockPtr&) const;

   void LoadNeighborsFromReader(BlockPtr&) const;

   void LoadFieldsFromReader(BlockPtr&) const;

public:

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

//! Make shared block
   virtual void MakeSharedBlock(BlockPtr &block_new);

//! Load interpolation stencil using interior zones
   void InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi, const GeoVector offset_lo, const GeoVector offset_hi, const GeoVector delta);

//! Generate an interpolation stencil in 3D
   virtual int BuildInterpolationStencil(const GeoVector& pos);

//! Find block order in the cache or get a block from the server if not in the cache
   int RequestBlock(void);

//! Get variables directly from data reader
   template <typename Fields>
   void GetVariablesFromReader(Fields& fields);

//! Get variables using 0th order interpolation
   template <typename Fields>
   void GetVariablesInterp0(const GeoVector& pos, Fields& fields);

//! Get variables using 1st order interpolation
   template <typename Fields>
   void GetVariablesInterp1(const GeoVector& pos, Fields& fields);

//! Get gradients using 1st order interpolation
   template <typename Fields>
   void GetGradientsInterp1(Fields& fields, DerivativeData& ddata);

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

//! Front end set up prior to main loop
   void Start(void) override;

//! Front end clean up tasks after the main loop
   void Finish(void) override;

//! Obtain the variables
   template <typename Fields>
   void GetVariables(double t, const GeoVector& pos, Fields& fields, double& dmax);

//! Obtain the gradients
   template <typename Fields>
   void GetGradients(Fields& fields, DerivativeData& ddata);

//! Print how many times internal/external interpolators were used
   void PrintStencilOutcomes(void) override;

//! Print how many blocks were requested
   void PrintNumBlocksRequested(void);

};


};

// Something like this is needed for templated classes
#include "background_data_cartesian.cc"

#endif
