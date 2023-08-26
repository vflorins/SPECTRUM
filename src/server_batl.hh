/*!
\file server_batl.hh
\brief Defines a class of a data server from BATL
\author Vladimir Florinski
*/

#ifndef SPECTRUM_SERVER_BATL_HH
#define SPECTRUM_SERVER_BATL_HH

#include "server_cartesian.hh"

namespace Spectrum {

#ifdef __cplusplus
extern "C" {
#endif

void wrapamr_clean();
void wrapamr_read_header(const char*, int, int);
void wrapamr_read_file(const char*, int, int, int);
void wrapamr_get_ndim(int*);
void wrapamr_get_nvar(int*);
void wrapamr_get_domain(double*, double*);
void wrapamr_get_block_size(int*, int*, int*, int*);
void wrapamr_get_data_serial(const double* pos, double* variables, int*);

void spectrum_init_mpi(int comm);
void spectrum_get_node(const double* pos, int* node);
void spectrum_get_neighbor_node(int node, int i, int j, int k, int* neighbor_node, int* neighbor_level);
void spectrum_get_interpolation_stencil(const double* pos, int* n_nodes, int* stencil_nodes, int* stencil_zones, double* stencil_weights);

// FIXME debug only
void spectrum_read_header(const char*, int, int);

#ifdef __cplusplus
}
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATL class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

class ServerBATL : virtual public ServerCartesian {

protected:

//! Initialize block pointer
   void InitializeBlockPtr(BlockBase* &block_ptr) override;

//! Default constructor - disabled
   ServerBATL(void) = default;

public:

//! Destructor
   virtual ~ServerBATL() = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATLFront class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

class ServerBATLFront : virtual public ServerBATL, virtual public ServerCartesianFront {

protected:

//! Generate an interpolation stencil for one plane
   int BuildInterpolationPlane(const GeoVector& pos, int plane, int half);

//! Generate an interpolation stencil in 3D
   int BuildInterpolationStencil(const GeoVector& pos) override;

//! Make shared block
   void MakeSharedBlock(BlockPtrType &block_new) override;

//! Obtain an interpolation stencil from the server
   void RequestStencil(void) override;

//! Obtain all variables at the specified location from the server
   void RequestVariables(double* vars) override;

public:

//! Default constructor
   ServerBATLFront(void) = default;

//! Destructor
   ~ServerBATLFront() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATLBack class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

class ServerBATLBack : virtual public ServerBATL, virtual public ServerCartesianBack {

public:

//! Default constructor
   ServerBATLBack(void) = default;

//! Constructor with arguments
   ServerBATLBack(const std::string& file_name_pattern_in);

//! Destructor
   ~ServerBATLBack() override = default;

//! Read data file
   void ReadData(const std::string data_file) override;

//! Clean reader
   void CleanReader(void) override;

//! Backend tasks during the main loop
   int ServerFunctions(void) override;

//! Serve block function
   void GetBlock(const double* pos, int* node) override;

//! Handle "needstencil" requests
   void HandleNeedStencilRequests(void);

//! Handle "needvars" requests
   void HandleNeedVarsRequests(void);
};

//! Server types
#if SERVER_TYPE == SERVER_BATL
typedef ServerBATL ServerType;
typedef ServerBATLFront ServerFrontType;
typedef ServerBATLBack ServerBackType;
#endif

};

#endif
