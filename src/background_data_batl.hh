/*!
\file background_data_batl.hh
\brief Declares a background class using data from BATL adaptive mesh on distributed memory
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_BATL_HH
#define SPECTRUM_BACKGROUND_DATA_BATL_HH

#include "background_data_base.hh"
#include "server_interface.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDataBATL class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background delivered by a BATL server
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

*/
template <typename HConfig_>
class BackgroundDataBATL : public BackgroundDataBase<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDataBATL";

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

   using BackgroundDataBase = BackgroundDataBase<HConfig>;
   using Block = BackgroundDataBase::Block;
   using BlockPtr = BackgroundDataBase::BlockPtr;
   using MPI = BackgroundDataBase::MPI;
   using DataFields = BackgroundDataBase::DataFields;

   using ServerInterface = BackgroundDataBase::ServerInterface;
   using ServerInterface::_inquiry;
   using ServerInterface::stencil;
   using ServerInterface::MPIInquiryType;
   using ServerInterface::MPIStencilType;
   using ServerInterface::MPIBlockType;

   using BackgroundDataBase::cache_line;
   using BackgroundDataBase::block_pri;
   using BackgroundDataBase::block_sec;
   using BackgroundDataBase::block_stn;
   using BackgroundDataBase::stencil_status;
   using BackgroundDataBase::stencil_outcomes;

// base methods:
   using BackgroundDataBase::RequestBlock;
   using BackgroundDataBase::InteriorInterpolationStencil;

// constexpr values:
   using BackgroundDataBase::allow_server_worker;
   using BackgroundDataBase::num_ghost_cells;
   using BackgroundDataBase::server_interpolation_order;
   using BackgroundDataBase::dmax0;

// debug methods:
   using BackgroundDataBase::PrintStencilOutcomes;
   using BackgroundDataBase::PrintNumBlocksRequested;
   using BackgroundDataBase::GetNCachedBlocks;
   using BackgroundDataBase::InvalidateCache;

   static constexpr bool request_stencil_from_batl = Config::request_stencil_from_batl;

/*!
\brief Convert a node multi-index to a level multi-index
\param[in] node_idx Node index (0-3)
\return Level index: 0->0, 1->1, 2->1, 3->2
*/
   static constexpr MultiIndex NodeToLevel(const MultiIndex& node_idx)
   {
      return (node_idx + 1) / 2;
   };

/*!
\brief Convert a level multi-index to a node multi-index
\param[in] level_idx Level index (0-2)
\return Node index: 0->0, 1->2, 2->3
*/
   static constexpr MultiIndex LevelToNode(const MultiIndex& level_idx)
   {
      return 2 * level_idx - level_idx / 2;
   };


protected:

   //! Obtain an interpolation stencil from the server
   int RequestStencil(const GeoVector& pos);

//! Generate an interpolation stencil for one plane
   int BuildInterpolationPlane(const GeoVector& pos, int plane, int half);

//! Generate an interpolation stencil in 3D
   int BuildInterpolationStencil(const GeoVector& pos) final;

public:

//! Default constructor
   BackgroundDataBATL(void) = default;

//! Destructor
   ~BackgroundDataBATL() override = default;

public: // data background API:

//! Front end set up prior to main loop
   void Start(void);

//! Front end clean up tasks after the main loop
   void Finish(void);

public: // general background API:

//! Compute the maximum distance per time step
   template <typename Coordinates>
   status_t EvaluateDmax(Coordinates&, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};


};


// Something like this is needed for templated classes
#include "background_data_batl.cc"

#endif
