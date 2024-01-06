/*!
\file server_base.cc
\brief Implements a base class of a data server from an external source
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_base.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/25/2022
\param[in] mpi_config_in A pointer to an MPI_Config object
*/
void ServerBase::ConnectMPIConfig(const std::shared_ptr<MPI_Config> mpi_config_in)
{
   mpi_config = mpi_config_in;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
*/
void ServerBase::ServerStart(void)
{
// Set up MPI data type for "Inquiry"
   MPI_Datatype inquiry_types[] = {MPI_INT, MPI_INT, MPI_DOUBLE};
   int inquiry_lengths[] = {1, 1, 3};
   MPI_Aint inquiry_displ[3];

// Figure out field displacements using "_inquiry" as template
   MPI_Get_address(&_inquiry.type, &inquiry_displ[0]);
   MPI_Get_address(&_inquiry.node, &inquiry_displ[1]);
   MPI_Get_address(&_inquiry.pos , &inquiry_displ[2]);
   for(auto i = 2; i >= 0; i--) inquiry_displ[i] -= inquiry_displ[0];

// Commit the type
   MPI_Type_create_struct(3, inquiry_lengths, inquiry_displ, inquiry_types, &MPIInquiryType);
   MPI_Type_commit(&MPIInquiryType);
};

/*!
\author Vladimir Florinski
\date 10/27/2022
*/
void ServerBase::ServerFinish(void)
{
   MPI_Type_free(&MPIInquiryType);
};

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
int ServerBase::ServerFunctions(void)
{
   return 0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBaseFront methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
void ServerBaseFront::ServerStart(void)
{
   ServerBase::ServerStart();
   cache_line.Empty();
};

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
void ServerBaseFront::ServerFinish(void)
{
   MPI_Send(nullptr, 0, MPI_BYTE, 0, tag_stopserve, mpi_config->node_comm);
   ServerBase::ServerFinish();
};

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
int ServerBaseFront::ServerFunctions(void)
{
   return 0;
};

/*!
\author Vladimir Florinski
\date 06/19/2020
\return Number of blocks in the cache
*/
int ServerBaseFront::GetNCachedBlocks(void) const
{
   return cache_line.size();
};

/*!
\author Vladimir Florinski
\date 01/27/2023
*/
void ServerBaseFront::InvalidateCache(void)
{
   cache_line.Empty();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBaseBack methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/26/2022
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
ServerBaseBack::ServerBaseBack(const std::string& file_name_pattern_in) 
              : ServerBase(),
                file_name_pattern(file_name_pattern_in)
{
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
void ServerBaseBack::ServerStart(void)
{
   ServerBase::ServerStart();

// Indices of processes returned by "Testsome"
   index_needblock   = new int[mpi_config->node_comm_size];
   index_needstencil = new int[mpi_config->node_comm_size];
   index_needvars    = new int[mpi_config->node_comm_size];
   index_stopserve   = new int[mpi_config->node_comm_size];

// Request arrays
   req_needblock   = new MPI_Request[mpi_config->node_comm_size];
   req_needstencil = new MPI_Request[mpi_config->node_comm_size];
   req_needvars    = new MPI_Request[mpi_config->node_comm_size];
   req_stopserve   = new MPI_Request[mpi_config->node_comm_size];

// Message buffers
   buf_needblock   = new Inquiry[mpi_config->node_comm_size];
   buf_needstencil = new Inquiry[mpi_config->node_comm_size];
   buf_needvars    = new Inquiry[mpi_config->node_comm_size];

// Post initial receives for all request types
   req_needblock[0] = MPI_REQUEST_NULL;
   req_needstencil[0] = MPI_REQUEST_NULL;
   req_needvars[0] = MPI_REQUEST_NULL;
   req_stopserve[0] = MPI_REQUEST_NULL;
   for(int cpu = 1; cpu < mpi_config->node_comm_size; cpu++) {
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, tag_needblock, mpi_config->node_comm, &req_needblock[cpu]);
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, tag_needstencil, mpi_config->node_comm, &req_needstencil[cpu]);
      MPI_Irecv(&buf_needvars[cpu], 1, MPIInquiryType, cpu, tag_needvars, mpi_config->node_comm, &req_needvars[cpu]);
      MPI_Irecv(nullptr, 0, MPI_BYTE, cpu, tag_stopserve, mpi_config->node_comm, &req_stopserve[cpu]);
   };
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
void ServerBaseBack::ServerFinish(void)
{
// Deallocate index arrays
   delete[] index_needblock;
   delete[] index_needstencil;
   delete[] index_needvars;
   delete[] index_stopserve;

// Deallocate request arrays
   delete[] req_needblock;
   delete[] req_needstencil;
   delete[] req_needvars;
   delete[] req_stopserve;
   
// Deallocate buffers
   delete[] buf_needblock;
   delete[] buf_needstencil;
   delete[] buf_needvars;

   ServerBase::ServerFinish();
};

/*!
\author Vladimir Florinski
\date 10/27/2022
\return Unused
*/
int ServerBaseBack::ServerFunctions()
{
   return 0;
};

};
