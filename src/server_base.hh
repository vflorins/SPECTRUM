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
#include "common/spatial_data.hh"
#include "cache_lru.hh"
#include <memory>

namespace Spectrum {

//! Index of the mass density variable
// #define SERVER_VAR_INDEX_RHO 0

//! Index of the number density variable
// #define SERVER_VAR_INDEX_DEN 0

//! Index of the momentum variable
//#define SERVER_VAR_INDEX_MOM 1

//! Index of the bulk flow variable
#define SERVER_VAR_INDEX_FLO 0

//! Index of the magnetic field variable
#define SERVER_VAR_INDEX_MAG 3

//! Index of the electric field variable
#define SERVER_VAR_INDEX_ELE 6

//! Index and number of the regions variable
// #define SERVER_VAR_INDEX_REG 7
// #define SERVER_NUM_INDEX_REG 2

//! Index of thermal pressure
// #define SERVER_VAR_INDEX_PTH 9

//! Unit of length
const double unit_length_server = unit_length_fluid;
// const double unit_length_server = 1.4959787e+13;

//! Unit of number density
const double unit_number_density_server = 1.0;

//! Unit of velocity
const double unit_velocity_server = unit_velocity_fluid;
// const double unit_velocity_server = 1.0E5;

//! Unit of magnetic field
const double unit_magnetic_server = unit_magnetic_fluid;
// const double unit_magnetic_server = 1.0E-5;

//! Unit of electric field
const double unit_electric_server = unit_electric_fluid;
// const double unit_electric_server = unit_velocity_server * unit_magnetic_server / unit_velocity_fluid;

//! Unit of pressure
const double unit_pressure_server = 1.0;

//! MPI tag for "need block" message (W->B)
const int tag_needblock = 1001;

//! MPI tag for "send block" message (B->W)
const int tag_sendblock = 1002;

//! MPI tag for "need stencil" message (W->B)
const int tag_needstencil = 1003;

//! MPI tag for "send stencil" message (B->W)
const int tag_sendstencil = 1004;

//! MPI tag for "need vars" message (W->B)
const int tag_needvars = 1005;

//! MPI tag for "send vars" message (B->W)
const int tag_sendvars = 1006;

//! MPI tag for "stop serve" message (W->B)
const int tag_stopserve = 1007;

/*!
\brief Data inquiry type
\author Vladimir Florinski
*/
struct Inquiry {

//! Type of inquiry: "0" is by ID, "1" is by position
   int type;

//! Node index if requesting by ID
   int node;

//! Position if requesting by position
   GeoVector pos;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if server functions failed
\author Juan G Alonso Guzman
*/
class ExServerError : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 23/08/2023
\return Text describing the error
*/
inline const char* ExServerError::what(void) const noexcept
{
   return "Server error";
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Common functions of the server frontend and backend
\author Vladimir Florinski

The Trajectory objects on the worker processes request state variables at arbitrary locations. The Server object on the worker processes fulfills these requests from its cache. If the blocks needed are not cached, they are requested from the Server object on the Boss process, which obtains them from its external interface and sends them to requesting processes via MPI. Each worker proceess has its own cache line, except for the server process, which has no cache.
*/
class ServerBase {

protected:

//! MPI configuration object
   std::shared_ptr<MPI_Config> mpi_config;

//! Number of variables
   int n_variables = 0;

//! Smallest position in the domain
   GeoVector domain_min;

//! Largest position in the domain
   GeoVector domain_max;

//! Current inquiry
   Inquiry _inquiry;

//! MPI data type for the "Inquiry" class
   MPI_Datatype MPIInquiryType;

//! MPI data type for the "Block" class
   MPI_Datatype MPIBlockType;

//! MPI data type for the "Stencil" class
   MPI_Datatype MPIStencilType;

//! Default constructor
   ServerBase(void) = default;

public:

//! Destructor
   virtual ~ServerBase(void) = default;

//! Connect to an existing MPI_Config object
   void ConnectMPIConfig(const std::shared_ptr<MPI_Config> mpi_config_in);

//! Common set up prior to main loop
   virtual void ServerStart(void);

//! Common clean up tasks after the main loop
   virtual void ServerFinish(void);

//! Common tasks during the main loop
   virtual int ServerFunctions(void);

//! Return the vector to one of the corners of the block
   GeoVector GetDomainMin(void) const;

//! Return the vector to the corner opposite to that returned in "GetDomainMin()"
   GeoVector GetDomainMax(void) const;
};

/*!
\author Vladimir Florinski
\date 06/19/2020
\return Coordinates closest to the origin
*/
inline GeoVector ServerBase::GetDomainMin(void) const
{
   return domain_min;
};

/*!
\author Vladimir Florinski
\date 06/19/2020
\return Coordinates farthest from the origin
*/
inline GeoVector ServerBase::GetDomainMax(void) const
{
   return domain_max;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBaseFront class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Server frontend (local) specific functionality
\author Vladimir Florinski
*/
class ServerBaseFront : virtual public ServerBase {

protected:

//! Cache line
   BlockCache cache_line;

//! Default constructor
   ServerBaseFront(void) = default;

#ifdef NEED_SERVER
//! Generate an interpolation stencil in 3D
   virtual int BuildInterpolationStencil(const GeoVector& pos) = 0;

//! Find block order in the cache or get a block from the server if not in the cache
   virtual int RequestBlock(void) = 0;
#endif

public:

//! Destructor
   ~ServerBaseFront(void) = default;

//! Front end set up prior to main loop
   void ServerStart(void) override;

//! Front end clean up tasks after the main loop
   void ServerFinish(void) override;

//! Frontend tasks during the main loop
   int ServerFunctions(void) override;

//! Return the number of cached blocks
   int GetNCachedBlocks(void) const;

//! Empty the cache
   void InvalidateCache(void);

#ifdef NEED_SERVER
//! Obtain the variables
   virtual void GetVariables(double t, const GeoVector& pos, SpatialData& spdata) = 0;

//! Obtain the gradients
   virtual void GetGradients(SpatialData& spdata) = 0;
#endif

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBaseBack class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Server backend (remote) specific functionality
\author Vladimir Florinski
*/
class ServerBaseBack : virtual public ServerBase {

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
   BlockBase* block_served = nullptr;

//! Default constructor
   ServerBaseBack(void) = default;

//! Constructor with arguments
   ServerBaseBack(const std::string& file_name_pattern_in);

public:

//! Destructor
   ~ServerBaseBack(void) = default;

//! Backend set up prior to main loop
   void ServerStart(void) override;

//! Backend clean up tasks after the main loop
   void ServerFinish(void) override;

//! Backend tasks during the main loop
   int ServerFunctions(void) override;
};

// FIXME: perhaps there is a better way to treat this case
//! Server types
#if SERVER_TYPE == SERVER_SELF
typedef ServerBase ServerType;
typedef ServerBaseFront ServerFrontType;
typedef ServerBaseBack ServerBackType;
#endif

};

#endif
