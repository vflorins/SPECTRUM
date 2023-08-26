/*!
\file server_cartesian.hh
\brief Defines a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman
*/

#ifndef SPECTRUM_SERVER_CARTESIAN_HH
#define SPECTRUM_SERVER_CARTESIAN_HH

#include "server_base.hh"

namespace Spectrum {

/*!
\brief Interpolation stencil for Cartesian
\author Vladimir Florinski
\author Juan G Alonso Guzman
*/
struct StencilCartesian {

//! Number of elements in the stencil
   int n_elements = 8;

//! List of blocks (can also be used to store global nodes)
   int blocks[8];

//! List of zones
   MultiIndex zones[8];

//! Weights
   double weights[8];

//! Partial derivatives of weights
   double derivatives[24];

//! Print the list of blocks, zones, and weights
   void Print(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

class ServerCartesian : virtual public ServerBase {

protected:

//! A working stencil variable
   StencilCartesian stencil;

//! Default constructor - disabled
   ServerCartesian(void) = default;

//! Initialize block pointer
   virtual void InitializeBlockPtr(BlockBase* &block_ptr);

//! Create block datatype
   void CreateBlockDatatype(void);

//! Create stencil datatype
   void CreateStencilDatatype(void);

//! Create Cartesian specific datatypes
   void CreateMPIDatatypes(void);

public:

//! Destructor
   virtual ~ServerCartesian() = default;

//! Common set up prior to main loop
   void ServerStart(void) override;

//! Common clean up tasks after the main loop
   void ServerFinish(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianFront class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

class ServerCartesianFront : virtual public ServerCartesian, virtual public ServerBaseFront {

protected:

//! Status of the most recently computed stencil
   int stencil_status = 0;

//! Counts of different stencil outcomes
   int stencil_outcomes[2];

//! Load interpolation stencil using interior zones
   void InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi, const GeoVector offset_lo, const GeoVector offset_hi, const GeoVector delta);

//! Generate an interpolation stencil in 3D
   int BuildInterpolationStencil(const GeoVector& pos) override;

//! Make shared block
   virtual void MakeSharedBlock(BlockPtrType &block_new);

//! Find block order in the cache or get a block from the server if not in the cache
   int RequestBlock(void) override;

//! Obtain an interpolation stencil from the server
   void RequestStencil(void) override;

//! Obtain all variables at the specified location from the server
   void RequestVariables(double* vars) override;

public:

//! Default constructor
   ServerCartesianFront(void) = default;

//! Destructor
   ~ServerCartesianFront() override = default;

//! Front end set up prior to main loop
   void ServerStart(void) override;

//! Front end clean up tasks after the main loop
   void ServerFinish(void) override;

//! Obtain the variables
   void GetVariables(double t, const GeoVector& pos, SpatialData& spdata) override;

//! Obtain the gradients
   void GetGradients(SpatialData& spdata) override;

//! Print how many times internal/external interpolators were used
   void PrintStencilOutcomes(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianBack class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

class ServerCartesianBack : virtual public ServerCartesian, virtual public ServerBaseBack {

public:

//! Default constructor
   ServerCartesianBack(void) = default;

//! Constructor with arguments
   ServerCartesianBack(const std::string& file_name_pattern_in);

//! Destructor
   ~ServerCartesianBack() override = default;

//! Read data file
   virtual void ReadData(const std::string data_file);

//! Back end set up prior to main loop
   void ServerStart(void) override;

//! Clean reader
   virtual void CleanReader(void);

//! Back end clean up tasks after the main loop
   void ServerFinish(void) override;

//! Backend tasks during the main loop
   int ServerFunctions(void) override;

//! Serve block function
   virtual void GetBlock(const double* pos, int* node);

//! Handle "needblock" requests
   void HandleNeedBlockRequests(void);

//! Handle "stopserve" requests
   int HandleStopServeRequests(void);
};

//! Server types
#if SERVER_TYPE == SERVER_CARTESIAN
typedef ServerCartesian ServerType;
typedef ServerCartesianFront ServerFrontType;
typedef ServerCartesianBack ServerBackType;
#endif

};

#endif