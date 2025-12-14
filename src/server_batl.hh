/*!
\file server_batl.hh
\brief Defines a class of a data server from BATL
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SERVER_BATL_HH
#define SPECTRUM_SERVER_BATL_HH

#include "server_base.hh"
#include "server_interface.hh"

namespace Spectrum {


#ifdef __cplusplus
extern "C" {
#endif

void wrapamr_clean();
void wrapamr_read_header(const char *, int, int);
void wrapamr_read_file(const char *, int, int, int);
void wrapamr_read_file_partial(const char *, int, int, int, int, int *);
void wrapamr_get_ndim(int *);
void wrapamr_get_nvar(int *);
void wrapamr_get_domain(double *, double *);
void wrapamr_get_block_size(int *, int *, int *, int *);
void wrapamr_get_data_serial(const double *pos, double *variables, int *);

void spectrum_init_mpi(int comm);
void spectrum_get_node(const double *pos, int *node);
void spectrum_get_neighbor_node(int node, int i, int j, int k, int *neighbor_node, int *neighbor_level);
void spectrum_get_interpolation_stencil(const double *pos, int *n_nodes, int *stencil_nodes, int *stencil_zones,
                                        double *stencil_weights);

// FIXME debug only
void spectrum_read_header(const char *, int, int);

#ifdef __cplusplus
}
#endif


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Interface for the Fortran routines in "server_batl.f90"
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

/*!
\brief Obtain the pointers to the neighbor arrays for a node
\param[in]  node            Node
\param[out] neighbor_nodes  Entry into the BATL array of neighbor nodes
\param[out] neighbor_levels Entry into the BATL array of neighbor levels
*/
//void spectrum_get_all_neighbor_nodes(int node, int** neighbor_nodes, int** neighbor_levels);

/*!
\brief Obtain copies of the neighbor arrays for a node
\param[in]  node            Node
\param[out] neighbor_nodes  Array of neighbor nodes
\param[out] neighbor_levels Array of neighbor levels
*/
void spectrum_get_all_neighbor_copies(int node, int* neighbor_nodes, int* neighbor_levels);

/*!
\brief Return the coordinates of two opposite corners of the block
\param[in]  node     Node
\param[out] face_min Coordinates of the lower corner
\param[out] face_max Coordinates of the upper corner
*/
void spectrum_get_block_corners(int node, double* face_min, double* face_max);

/*!
\brief Obtain copies of the variables associated with a node
\param[in]  node      Node
\param[out] variables Array of variables
*/
void spectrum_get_block_data(int node, double* variables);

#ifdef __cplusplus
}
#endif


//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATL class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <typename HConfig_>
class ServerBATL : public ServerBase<HConfig_>, public ServerInterface<HConfig_> {
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

   static constexpr int server_interpolation_order = Config::server_interpolation_order;
   static constexpr int num_ghost_cells = Config::num_ghost_cells;
   static constexpr bool request_stencil_from_batl = Config::request_stencil_from_batl;

public: // API:

//! Default constructor
   ServerBATL(void) = default;

//! Destructor
   ~ServerBATL() override = default;

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
   void GetBlockData(const double* pos, double* vars, int* found) final;

//! Serve block function
   void GetBlock(const double* pos, int* node) final;

//! Serve stencil
   void GetStencil(const double* pos, Stencil* stencil) final;

private: // helpers

//! Read data file
   void ReadData(const std::string data_file);

//! Clean reader
   void CleanReader(void);

};

};

#include "server_batl.cc"

#endif
