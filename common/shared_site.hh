/*!
\file shared_site.hh
\brief Implements a communication object controlling a shared site
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SHARED_SITE_HH
#define SPECTRUM_SHARED_SITE_HH

#include <common/communication_site.hh>

#ifdef USE_MPI

namespace Spectrum {

/*
\brief Minimum MPI operation template
\author Vladimir Florinski
\date 01/17/2025
\param[in]     first        First array
\param[in,out] second       Second array
\param[in]     length       Length of array
\param[in]     mpi_datatype MPI data type
*/
template <typename datatype>
void MPI_User_Min(datatype* first, datatype* second, int* length, MPI_Datatype* mpi_datatype)
{
   for(auto i = 0; i < *length; i++) {
      *second = std::min(*first, *second);
      first++;
      second++;
   };
};

/*
\brief Maximum MPI operation template
\author Vladimir Florinski
\date 01/17/2025
\param[in]     first        First array
\param[in,out] second       Second array
\param[in]     length       Length of array
\param[in]     mpi_datatype MPI data type
*/
template <typename datatype>
void MPI_User_Max(datatype* in, datatype* inout, int* len, MPI_Datatype* mpi_datatype)
{
   for(auto i = 0; i < *len; i++) {
      *inout = std::max(*in, *inout);
      in++;
      inout++;
   };
};

/*
\brief Maximum MPI operation template
\author Vladimir Florinski
\date 01/17/2025
\param[in]     first        First array
\param[in,out] second       Second array
\param[in]     length       Length of array
\param[in]     mpi_datatype MPI data type
*/
template <typename datatype>
void MPI_User_Sum(datatype* in, datatype* inout, int* len, MPI_Datatype* mpi_datatype)
{
   for(auto i = 0; i < *len; i++) {
      *inout = *in + *inout;
      in++;
      inout++;
   };
};

/*!
\brief A class representing an MPI shared site
\author Vladimir Florinski
*/
template <typename datatype, int op>
struct SharedSite : public CommunicationSite<datatype>
{
protected:

   using CommunicationSite<datatype>::site_comm;
   using CommunicationSite<datatype>::site_comm_size;
   using CommunicationSite<datatype>::buffer;
   using CommunicationSite<datatype>::sendcounts;
   using CommunicationSite<datatype>::sdispls;
   using CommunicationSite<datatype>::mpi_datatype;
   using CommunicationSite<datatype>::n_parts;
   using CommunicationSite<datatype>::buf_size;

//! MPI reduction operation
   MPI_Op mpi_red_func = MPI_OP_NULL;

//! Create an operation for the user-defined datatype
   void CreateOp(void);

public:

//! Default constructor
   SharedSite(void);

//! Copy constructor - deleted because each site must be unique due to MPI restrictions
   SharedSite(const SharedSite& other) = delete;

//! Move constructor
   SharedSite(SharedSite&& other);

//! Assignment operator - deleted because we cannot have multiple copies of the MPI objects
   SharedSite& operator =(const SharedSite& other) = delete;

//! Perform the sharing
   void Share(void);
};

//! Dummy specialization for void datatype
template <> class SharedSite<void> {};

/*!
\author Vladimir Florinski
\date 01/16/2025
*/
template <typename datatype, int op>
inline SharedSite<datatype, op>::SharedSite(void)
                           : CommunicationSite<datatype>()
{
   switch(op) {

   case 0:
      MPI_Op_create((MPI_User_function*)MPI_User_Min<datatype>, true, &mpi_red_func);
      break;

   case 1:
      MPI_Op_create((MPI_User_function*)MPI_User_Max<datatype>, true, &mpi_red_func);
      break;
   
   case 2:
      MPI_Op_create((MPI_User_function*)MPI_User_Sum<datatype>, true, &mpi_red_func);
      break;

   default:
      break;
   };
};

/*!
\author Vladimir Florinski
\date 01/16/2025
\param[in] other Object to move into this
*/
template <typename datatype, int op>
inline SharedSite<datatype, op>::SharedSite(SharedSite<datatype, op>&& other)
                               : CommunicationSite<datatype>(std::move(static_cast<CommunicationSite<datatype>&&>(other)))
{
   mpi_red_func = other.mpi_red_func;
   other.mpi_red_func = MPI_OP_NULL;
};

/*!
\author Vladimir Florinski
\date 01/17/2025
*/
template <typename datatype, int op>
inline void SharedSite<datatype, op>::Share(void)
{
// If this process is the only one participating, the data is accessible directly from the buffer.
   if ((site_comm == MPI_COMM_NULL) || (site_comm_size <= 1)) return;

   int my_rank, n_parts_rank;
   datatype* buf_ptr;
   datatype* buf_reduced = new datatype[buf_size];
   MPI_Comm_rank(site_comm, &my_rank);

   n_parts_rank = sendcounts[my_rank] / buf_size;
   
// Perform local reduction first
   for (auto i = 0; i < buf_size; i++) {

// Very first entry - we cannot use 0 or some other number because the starting value depends on "op".
      buf_reduced[i] = *(buffer + sdispls[my_rank]);
      for (auto part_rank = 1; part_rank < n_parts_rank; part_rank++) {
         buf_ptr = buffer + sdispls[my_rank] + part_rank * buf_size;

         switch(op) {

// minimum
         case 0:
            buf_reduced[i] = std::min(buf_reduced[i], buf_ptr[i]);
            break;

// maximum
         case 1:
            buf_reduced[i] = std::max(buf_reduced[i], buf_ptr[i]);
            break;

// average
         case 2:
            buf_reduced[i] += buf_ptr[i];
            break;
            
         default:
            break;
         };
      };
   };
            
// Perform global reduction
   MPI_Allreduce(MPI_IN_PLACE, buf_reduced, buf_size, mpi_datatype, mpi_red_func, site_comm);

// Put the numbers back in "buffer"
   for (auto part_rank = 0; part_rank < n_parts_rank; part_rank++) {
      buf_ptr = buffer + sdispls[my_rank] + part_rank * buf_size;
      std::memcpy(buf_ptr, buf_reduced, buf_size * sizeof(datatype));

// For averages, an extra step is requires
      if(op == 2) {
         for (auto i = 0; i < buf_size; i++) buf_ptr[i] /= n_parts;
      };
   };

   delete[] buf_reduced;
};

#ifdef GEO_DEBUG

/*!
\brief Test the functionality of this class
\author Vladimir Florinski
\date 01/17/2025
\param[in] test_rank Process rank that will do the testing
*/
inline void TestShared(int test_rank)
{
   const int n_parts = 8;
   const int buf_size = 1;
   int ranks[n_parts];
   int labels[n_parts];

// Make sure the user-supplied "test_rank" is sane
   if ((test_rank < 0) || (test_rank >= MPI_Config::glob_comm_size)) return;

   if (test_rank == MPI_Config::glob_comm_rank) {
      std::cerr << std::endl;
      std::cerr << "Testing exchange for global rank " << test_rank << std::endl;
      std::cerr << "Global communicator size: " << MPI_Config::glob_comm_size << std::endl;
      std::cerr << "--------------------------------------------------------------------------------\n";
      std::cerr << "Printing input values before site initialization\n";
   };

// Generate a test topology (some processes will not be in it)
   srand48(382667);
   for (auto part = 0; part < n_parts; part++) {

// Pick ranks at random - this would typically create a few duplicate entries, which is good for testing.
      ranks[part] = drand48() * MPI_Config::glob_comm_size;

// Make labels complimentary to parts
      labels[part] = n_parts - part - 1;
      if (test_rank == MPI_Config::glob_comm_rank) {
         std::cerr << "Part" << std::setw(3) << part << "   Global rank" << std::setw(3) << ranks[part] << "   Label" << std::setw(3) << labels[part] << std::endl;
      };
   };

// Create and set up an instance of ExchangeSite
   SharedSite<int, 0> exch_site;
   int* buf_ptr;
   exch_site.SetUpProperties(0, n_parts, buf_size, ranks, labels);

   if (test_rank == MPI_Config::glob_comm_rank) {
      std::cerr << "--------------------------------------------------------------------------------\n";
      exch_site.PrintSiteInfo(test_rank);
   };

// Fill the buffers with test data - can only do that on the ranks that own the respective part
   if (exch_site.GetCommSize() > 0) {
      for (auto part = 0; part < n_parts; part++) {
         buf_ptr = exch_site.BufferAddress(part);
         if (MPI_Config::glob_comm_rank == ranks[part]) *buf_ptr = 10 * drand48() * MPI_Config::glob_comm_size;
         else *buf_ptr = -1;
      };
   }
   else if (test_rank == MPI_Config::glob_comm_rank) {
      std::cerr << "Process not in the site communicator\n";
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Print buffers before the exchange
   if ((test_rank == MPI_Config::glob_comm_rank) && (exch_site.GetCommSize() > 0)) {
      std::cerr << "--------------------------------------------------------------------------------\n";
      std::cerr << "Buffer before exchange:\n";
      for (auto part = 0; part < exch_site.GetPartCount(); part++) {
         buf_ptr = exch_site.BufferAddress(part);
         std::cerr << "Part" << std::setw(3) << part << "   Value" << std::setw(3) << *buf_ptr << std::endl;
      };
//      std::cerr << std::endl;
   };      

// To properly test we must call this on all processes, not just those in "site_comm".
   exch_site.Share();

// Print buffers after the exchange
   if ((test_rank == MPI_Config::glob_comm_rank) && (exch_site.GetCommSize() > 0)) {
      std::cerr << "--------------------------------------------------------------------------------\n";
      std::cerr << "Buffer after exchange:\n";
      for (auto part = 0; part < exch_site.GetPartCount(); part++) {
         buf_ptr = exch_site.BufferAddress(part);
         std::cerr << "Part" << std::setw(3) << part << "   Value" << std::setw(3) << *buf_ptr << std::endl;
      };
      std::cerr << std::endl;
   };      
};

#endif

};

#endif

#endif
