/*!
\file exchange_site.hh
\brief Implements a communication object controlling one exchange site
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_EXCHANGE_SITE_HH
#define SPECTRUM_EXCHANGE_SITE_HH

#include <common/communication_site.hh>

#ifdef USE_MPI

#include <utility>

namespace Spectrum {

/*!
\brief A class representing an MPI exchange site
\author Vladimir Florinski
*/
template <typename datatype>
class ExchangeSite : public CommunicationSite<datatype>
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

public:

//! Default constructor
   ExchangeSite(void);

//! Copy constructor - deleted because each site must be unique due to MPI restrictions
   ExchangeSite(const ExchangeSite& other) = delete;

//! Move constructor
   ExchangeSite(ExchangeSite&& other);

//! Assignment operator - deleted because we cannot have multiple copies of the MPI objects
   ExchangeSite& operator =(const ExchangeSite& other) = delete;

//! Perform the exchange
   void Exchange(void);

#ifdef GEO_DEBUG

//! Fill all buffers with the same value
   void FillBuffers(datatype val);

//! Print the participant information
   void PrintParts(void) const;

#endif

};

/*!
\author Vladimir Florinski
\date 01/16/2025
*/
template <typename datatype>
inline ExchangeSite<datatype>::ExchangeSite(void)
                             : CommunicationSite<datatype>()
{
};

/*!
\author Vladimir Florinski
\date 01/16/2025
\param[in] other Object to move into this
*/
template <typename datatype>
inline ExchangeSite<datatype>::ExchangeSite(ExchangeSite<datatype>&& other)
                             : CommunicationSite<datatype>(std::move(static_cast<CommunicationSite<datatype>&&>(other)))
{
};

/*!
\author Vladimir Florinski
\date 06/26/2024
*/
template <typename datatype>
inline void ExchangeSite<datatype>::Exchange(void)
{
// If this process is the only one participating, the data is accessible directly from the buffer.
   if ((site_comm == MPI_COMM_NULL) || (site_comm_size <= 1)) return;

   MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buffer, sendcounts, sdispls, mpi_datatype, site_comm);
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 01/11/2025
\param[in] val Value to fill the buffers with
*/
template <typename datatype>
inline void ExchangeSite<datatype>::FillBuffers(datatype val)
{
   for(auto bufidx = 0; bufidx < buf_size * n_parts; bufidx++) buffer[bufidx] = val;
};

/*!
\brief Test the functionality of this class
\author Vladimir Florinski
\date 01/16/2025
\param[in] test_rank Process rank that will do the testing
*/
inline void TestExchange(int test_rank)
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
   ExchangeSite<int> exch_site;
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
         if (MPI_Config::glob_comm_rank == ranks[part]) *buf_ptr = part;
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
   exch_site.Exchange();

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
