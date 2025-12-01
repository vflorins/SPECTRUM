/*!
\file server_cartesian.hh
\brief Defines a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SERVER_CARTESIAN_HH
#define SPECTRUM_SERVER_CARTESIAN_HH

#include "server_base.hh"
#include "server_interface.hh"
#include "reader_cartesian.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <typename HConfig_>
class ServerCartesian : public ServerBase<HConfig_>, public ServerInterface<HConfig_> {
public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using ServerInterface = ServerInterface<HConfig>;
   using Block = Block<HConfig>;
   using Stencil = Stencil<Config::stencil_n_elements>;
   using MPI = MPI<HConfig, HConfig::MPI_enabled>;

   using DataFields = Config::DataFields;
   using BlockPtr = ServerInterface::BlockPtr;

   using ServerBase = ServerBase<HConfig>;
   using ServerBase::block_served;
   using ServerBase::buf_needblock;
   using ServerBase::buf_needstencil;
   using ServerBase::buf_needvars;
   using ServerBase::req_needblock;
   using ServerBase::req_needstencil;
   using ServerBase::req_needvars;
   using ServerBase::req_stopserve;
   using ServerBase::index_needblock;
   using ServerBase::index_needstencil;
   using ServerBase::index_needvars;
   using ServerBase::index_stopserve;
   using ServerBase::domain_max;
   using ServerBase::domain_min;
   using ServerBase::file_name_pattern;

   using ServerInterface::_inquiry;
   using ServerInterface::stencil;
   using ServerInterface::MPIInquiryType;
   using ServerInterface::MPIStencilType;
   using ServerInterface::MPIBlockType;

   static constexpr int server_interp_order = Config::server_interp_order;
   static constexpr int num_ghost_cells = Config::num_ghost_cells;

protected:

//! Cartesian data structure
   ReaderCartesian<HConfig> reader;

public: // API:

//! Default constructor
   ServerCartesian(void) = default;

//! Destructor
   ~ServerCartesian() = default;

//! Back end set up prior to main loop
   void ServerStart(void);

//! Back end clean up tasks after the main loop
   void ServerFinish(void);

//! Backend tasks during the main loop
   int ServerFunctions(void);

protected: // internal virtual dispatch API

   void LoadFromReader(BlockPtr&) final;

   void LoadNeighborsFromReader(BlockPtr&) final;

   void LoadFieldsFromReader(BlockPtr&) final;

//! Get block data from reader
   virtual void GetBlockData(const double* pos, double* vars, int* found) final;

//! Get block from reader
   virtual void GetBlock(const double* pos, int* node) final;

//! Serve stencil
   void GetStencil(const double* pos, Stencil* stencil) final;

private: // helpers

   //! Read data file
   virtual void ReadData(const std::string data_file);

//! Clean reader
   virtual void CleanReader(void);

};


};

#include "server_cartesian.cc"

#endif
