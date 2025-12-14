/*!
\file background_data_cartesian.hh
\brief Declares a background class using data from uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DATA_CARTESIAN_HH
#define SPECTRUM_BACKGROUND_DATA_CARTESIAN_HH

#include "background_data_base.hh"
#include "cache_lru.hh"
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
class BackgroundDataCartesian : public BackgroundDataBase<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDataCartesian";

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

// constexpr values:
   using BackgroundDataBase::allow_server_worker;
   using BackgroundDataBase::num_ghost_cells;
   using BackgroundDataBase::server_interpolation_order;
   using BackgroundDataBase::dmax0;

   using BackgroundDataBase::cache_line;
   using BackgroundDataBase::domain_min;
   using BackgroundDataBase::domain_max;
   using BackgroundDataBase::stencil_status;
   using BackgroundDataBase::stencil_outcomes;
   using BackgroundDataBase::block_pri;
   using BackgroundDataBase::block_sec;
   using BackgroundDataBase::block_stn;

// debug methods:
   using BackgroundDataBase::PrintStencilOutcomes;
   using BackgroundDataBase::PrintNumBlocksRequested;
   using BackgroundDataBase::GetNCachedBlocks;
   using BackgroundDataBase::InvalidateCache;

   // secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;


public:

//! Default constructor
   BackgroundDataCartesian(void) = default;

//! Destructor
   ~BackgroundDataCartesian() override = default;

//! Generate an interpolation stencil in 3D
   int BuildInterpolationStencil(const GeoVector& pos) final;

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
#include "background_data_cartesian.cc"

#endif
