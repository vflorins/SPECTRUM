/*!
\file server_base.cc
\brief Implements a base class of a data server from an external source
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

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
\date 10/26/2022
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
template <typename HConfig>
ServerBase<HConfig>::ServerBase(const std::string& file_name_pattern_in)
              : ServerInterface(),
                file_name_pattern(file_name_pattern_in)
{
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
template <typename HConfig>
void ServerBase<HConfig>::ServerStart(void)
{
   ServerInterface::ServerInterfaceStart();

// Indices of processes returned by "Testsome"
   index_needblock   = new int[MPI::node_comm_size];
   index_needstencil = new int[MPI::node_comm_size];
   index_needvars    = new int[MPI::node_comm_size];
   index_stopserve   = new int[MPI::node_comm_size];

// Request arrays
   req_needblock   = new MPI_Request[MPI::node_comm_size];
   req_needstencil = new MPI_Request[MPI::node_comm_size];
   req_needvars    = new MPI_Request[MPI::node_comm_size];
   req_stopserve   = new MPI_Request[MPI::node_comm_size];

// Message buffers
   buf_needblock   = new Inquiry[MPI::node_comm_size];
   buf_needstencil = new Inquiry[MPI::node_comm_size];
   buf_needvars    = new Inquiry[MPI::node_comm_size];

// Post initial receives for all request types
   req_needblock[0] = MPI_REQUEST_NULL;
   req_needstencil[0] = MPI_REQUEST_NULL;
   req_needvars[0] = MPI_REQUEST_NULL;
   req_stopserve[0] = MPI_REQUEST_NULL;
   for (int cpu = 1; cpu < MPI::node_comm_size; cpu++) {
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, MPI::tag::needblock, MPI::node_comm, &req_needblock[cpu]);
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, MPI::tag::needstencil, MPI::node_comm, &req_needstencil[cpu]);
      MPI_Irecv(&buf_needvars[cpu], 1, MPIInquiryType, cpu, MPI::tag::needvars, MPI::node_comm, &req_needvars[cpu]);
      MPI_Irecv(nullptr, 0, MPI_BYTE, cpu, MPI::tag::stopserve, MPI::node_comm, &req_stopserve[cpu]);
   };
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
template <typename HConfig>
void ServerBase<HConfig>::ServerFinish(void)
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

   ServerInterface::ServerInterfaceFinish();
};

/*!
\author Vladimir Florinski
\date 10/27/2022
\return Unused
*/
template <typename HConfig>
int ServerBase<HConfig>::ServerFunctions()
{
   return 0;
};

};
