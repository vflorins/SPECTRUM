/*!
\file server_batl.cc
\brief Implements a class of a data server from BATL
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_batl.hh"
#include "common/print_warn.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBATL public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/25/2022
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
template <typename HConfig>
ServerBATL<HConfig>::ServerBATL(const std::string& file_name_pattern_in)
              : ServerCartesian(file_name_pattern_in),
                ServerInterface(file_name_pattern_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/01/2023
*/
template <typename HConfig>
void ServerBATL<HConfig>::ReadData(const std::string data_file)
{
// Initialize the BATL library and read the data into memory
   MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_SELF);
   spectrum_init_mpi(fcomm);

   wrapamr_read_header(data_file.c_str(), line_width, 1);
   int vars_inds_fortran[DataFields::size()];
   for (int i = 0; i < DataFields::size(); i++) vars_inds_fortran[i] = i+1;
   wrapamr_read_file_partial(data_file.c_str(), line_width, 0, 0, DataFields::size(), vars_inds_fortran);
   wrapamr_get_domain(domain_min.Data(), domain_max.Data());
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
template <typename HConfig>
void ServerBATL<HConfig>::CleanReader(void)
{
   wrapamr_clean();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\return Number of clients that completed their tasks during this cycle

\note "needvars" requests are only handled when SERVER_INTERP_ORDER = -1. "needstencil" requests are only handled when using 1st order interpolation and either there are no ghost cells or stencils are always requested from BATL (for debugging purposes)
*/
template <typename HConfig>
int ServerBATL<HConfig>::ServerFunctions(void)
{
   if constexpr (server_interp_order == -1) {
// Handle "needvars" requests
      HandleNeedVarsRequests();
   }
   else if constexpr (server_interp_order == 1 && (num_ghost_cells == 0 || request_stencil_from_batl)) {
// Handle "needstencil" requests
      HandleNeedStencilRequests();
   }
// Handle "needblock" requests
   HandleNeedBlockRequests();
// Handle "stopserve" requests
   return HandleStopServeRequests();
};

/*!
\author Juan G Alonso Guzman
\date 08/03/2023
\param[in] pos    position array in reader coordinates
\param[out] vars  1D array containing all variables at pos
\param[out] found flag indicating whether variables where found or not for given position
*/
template <typename HConfig>
void ServerBATL<HConfig>::GetBlockData(const double* pos, double* vars, int* found)
{
   wrapamr_get_data_serial(pos, vars, found);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[in] pos    position array in reader coordinates
\param[out] node  ID tag of block containing pos within its physical boundaries
*/
template <typename HConfig>
void ServerBATL<HConfig>::GetBlock(const double* pos, int* node)
{
   spectrum_get_node(pos, node);
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
template <typename HConfig>
void ServerBATL<HConfig>::HandleNeedStencilRequests(void)
{
   GeoVector pos_batl;
   int cpu, cpu_idx, count_needstencil = 0;

// Service the "needstencil" requests
   MPI_Testsome(MPI::node_comm_size, req_needstencil, &count_needstencil, index_needstencil, MPI_STATUSES_IGNORE);

   for (cpu_idx = 0; cpu_idx < count_needstencil; cpu_idx++) {
      cpu = index_needstencil[cpu_idx];

// Obtain the stencil requested.
      pos_batl = buf_needstencil[cpu].pos / unit_length_server * unit_length_fluid;
      spectrum_get_interpolation_stencil(pos_batl.Data(), &stencil.n_elements, stencil.blocks, &stencil.zones[0][0], stencil.weights);

// Send the stencil to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(&stencil, 1, MPIStencilType, cpu, MPI::tag::sendstencil, MPI::node_comm);

// Post the receive for the next stencil request from this worker.
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, MPI::tag::needstencil, MPI::node_comm, &req_needstencil[cpu]);
   };
};

};
