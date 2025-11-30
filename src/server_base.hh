/*!
\file server_base.hh
\brief Defines a base class of a data server from an external source
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SERVER_BASE_HH
#define SPECTRUM_SERVER_BASE_HH

#include "common/mpi_config.hh"
#include "common/vectors.hh"
#include "common/physics.hh"
#include "server_interface.hh"
#include <memory>

#include "common/status.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Server backend (remote) specific functionality
\author Vladimir Florinski
\author Lucius Schoenbaum
*/
template <typename HConfig_>
class ServerBase : virtual public ServerInterface<HConfig_> {
public:

   using HConfig = HConfig_;
   using ServerInterface = ServerInterface<HConfig>;
   using Block  = Block<HConfig>;
   using MPI = MPI<HConfig>;

   using ServerInterface::_inquiry;
   using ServerInterface::stencil;
   using ServerInterface::MPIInquiryType;
   using ServerInterface::MPIStencilType;
   using ServerInterface::MPIBlockType;


protected:

//! File name pattern
   std::string file_name_pattern;

//! Indices of processes returned by "Testsome"
   int* index_needblock = nullptr;
   int* index_needstencil = nullptr;
   int* index_needvars = nullptr;
   int* index_stopserve = nullptr;

//! Request arrays
   MPI_Request* req_needblock = nullptr;
   MPI_Request* req_needstencil = nullptr;
   MPI_Request* req_needvars = nullptr;
   MPI_Request* req_stopserve = nullptr;

//! Buffers (not required for stopserve because the message is of zero length)
   Inquiry* buf_needblock = nullptr;
   Inquiry* buf_needstencil = nullptr;
   Inquiry* buf_needvars = nullptr;

//! Buffer for the block to be served
   Block* block_served = nullptr;

//! Smallest position in the domain
   GeoVector domain_min;

//! Largest position in the domain
   GeoVector domain_max;

//! Default constructor
   ServerBase(void) = default;

//! Constructor with arguments
   ServerBase(const std::string& file_name_pattern_in);

public:

//! Destructor
   ~ServerBase(void) = default;

//! Backend set up prior to main loop
   void ServerStart(void) override;

//! Backend clean up tasks after the main loop
   void ServerFinish(void) override;

//! Backend tasks during the main loop
   int ServerFunctions(void) override;

};

};

#include "server_base.cc"

#endif
