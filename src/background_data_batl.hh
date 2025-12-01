/*!
\file background_data_batl.hh
\brief Declares a background class using data from BATL adaptive mesh on distributed memory
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_BATL_HH
#define SPECTRUM_BACKGROUND_DATA_BATL_HH

#include "background_data_cartesian.hh"
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

Parameters: (BackgroundServerCartesian)
*/
template <typename HConfig_>
class BackgroundDataBATL : public BackgroundDataCartesian<HConfig_>, public ServerInterface<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDataBATL";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using Block = Block<HConfig>;
   using MPI = MPI<HConfig, HConfig::MPI_enabled()>;

   // todo as the following notes indicate, both classes can inherit independently from a BackgroundDataBase class.

   // todo base class
   using BackgroundDataCartesian = BackgroundDataCartesian<HConfig>;
   // todo base class
   using BackgroundDataCartesian::RequestBlock;
   // todo base class
   using BackgroundDataCartesian::cache_line;
   // todo base class
   using BackgroundDataCartesian::block_pri;
   // todo base class
   using BackgroundDataCartesian::block_sec;
   // todo base class
   using BackgroundDataCartesian::server_interp_order;
   // todo base class
   using BackgroundDataCartesian::num_ghost_cells;
   // todo base class
   using BackgroundDataCartesian::InteriorInterpolationStencil;

   using ServerInterface = ServerInterface<HConfig>;
   using ServerInterface::_inquiry;
   using ServerInterface::stencil;
   using ServerInterface::MPIInquiryType;
   using ServerInterface::MPIStencilType;
   using ServerInterface::MPIBlockType;
   using BlockPtr = ServerInterface::BlockPtr;

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

   // todo this is BATL-specific

   //! Obtain an interpolation stencil from the server
   int RequestStencil(const GeoVector& pos);

//! Generate an interpolation stencil for one plane
   int BuildInterpolationPlane(const GeoVector& pos, int plane, int half);

   // todo this is virtual in base class

//! Generate an interpolation stencil in 3D
   int BuildInterpolationStencil(const GeoVector& pos);

public:

//! Default constructor
   BackgroundDataBATL(void) = default;

//! Destructor
   ~BackgroundDataBATL() override = default;

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

};


};


// Something like this is needed for templated classes
#include "background_data_batl.cc"

#endif
