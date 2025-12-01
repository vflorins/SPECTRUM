/*!
\file server_interface.hh
\brief Defines a server-worker interface for a data-defined grid
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INTERFACE_HH
#define SPECTRUM_INTERFACE_HH

#include "server_types.hh"

namespace Spectrum {



//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerInterface class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief common methods for a server-worker interface
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

The Trajectory objects on the worker processes request state variables at arbitrary locations.
 The Server object on the worker processes fulfills these requests from its cache.
 If the blocks needed are not cached, they are requested from the Server object on the Server process,
 which obtains them from its external interface and sends them to requesting processes via MPI.
 Each worker proceess has its own cache line, except for the server process, which has no cache.
*/
template <typename HConfig_>
class ServerInterface {
public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using Block = Block<HConfig>;
   using Stencil = Stencil<Config::stencil_n_elements>;
   using BlockPtr = std::shared_ptr<Block>;

protected:

//! Current inquiry
   Inquiry _inquiry;

//! A working stencil variable
   Stencil stencil;

//! MPI data type for the "Inquiry" class
   MPI_Datatype MPIInquiryType;

//! MPI data type for the "Block" class
   MPI_Datatype MPIBlockType;

//! MPI data type for the "Stencil" class
   MPI_Datatype MPIStencilType;

   //! Default constructor - disabled
   ServerInterface(void) = default;

//! Create block datatype
   void CreateBlockDatatype(void);

//! Create inquiry datatype
   void CreateInquiryDatatype(void);

//! Create stencil datatype
   void CreateStencilDatatype(void);

public:

//! Destructor
   virtual ~ServerInterface() = default;

//! Common set up prior to main loop
   void ServerInterfaceStart(void);

//! Common clean up tasks after the main loop
   void ServerInterfaceFinish(void);

};

};

#include "server_interface.cc"

#endif
