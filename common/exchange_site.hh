/*!
\file exchange_site.hh
\brief Implements a communication object controlling one exchange site
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_EXCHANGE_SITE_HH
#define SPECTRUM_EXCHANGE_SITE_HH

#include <common/mpi_config.hh>

#ifdef USE_MPI

#include <map>
#include <vector>
#include <memory>
#include <algorithm>

#ifdef GEO_DEBUG
#include <common/print_warn.hh>
#endif

namespace Spectrum {

// Example:
//          -------------------------               
//          |           |           |     ranks_in:   11,  10,  11,  12 - in participating order
//          | part: 2   | part: 3   |
//          | rank: 11  | rank: 12  |     site_comm: 0->11, 1->10, 2->12
//          |           |           |     
//          ------------+------------     sendcounts: 2, 1, 1 (x buf_size)
//          |           |           |     sdispls:    0, 2, 3 (x buf_size)
//          | part: 0   | part: 1   |
//          | rank: 11  | rank: 10  |     buffer_entry_map: 0->0, 1->2, 2->1, 3->3
//          |           |           |
//          -------------------------

/*!
\brief A class representing an MPI exchange site
\author Vladimir Florinski

An object of this class contains an MPI communicator, a data buffer, and metadata to enable an exchange with MPI_Allgatherv. The array "buffer" has "n_part" slots corresponding to each participant. Because a process could host multiple participants, "ExchangeSite" will aggregate those into larger slots in the array (in MPI_Allgatherv each process sends and receives _one_ message, but the message size could be different). The intended use is for the calling function to create an array of ExchangeSite objects, one per site. This is redundant because each process may only need access to a subset of all sites. However, MPI communicator creation routines must be called on all processes, and it is desirable to perform all such operations within the class. The unused sites will remain, but since the communicator and buffers are not allocated, the memory lost will be small.
*/
template <typename datatype>
struct ExchangeSite
{
//! Index of this site
   int site_index;

//! Number of participants (actual)
   int n_parts;

//! Site communicator
   MPI_Comm site_comm = MPI_COMM_NULL;

//! Size of the site communicator
   int site_comm_size = 0;

//! MPI data type corresponding to the template type
   MPI_Datatype mpi_datatype = MPI_DATATYPE_NULL;

//! Shared buffer - should not be accessed directly by the caller, but via "buffer_entry"
   datatype* buffer = nullptr;

//! Send counts
   int* sendcounts = nullptr;

//! Displacements
   int* sdispls = nullptr;

//! Part index lookup map, accessed via labels
   std::map<int, int> part_lookup;

//! Entry points in the buffer, accessed via the part index
   std::map<int, datatype*> buffer_entry;

//! Default constructor
   ExchangeSite(void) = default;

//! Destructor
   ~ExchangeSite();

//! Import lists of ranks and part labels
   void SetUpProperties(int index_in, int n_parts_in, int buf_size_in, const int* ranks_in,
                        const int* labels_in, const std::shared_ptr<MPI_Config> mpi_config_in);

//! Perform the exchange
   void Exchange(void);
};

/*!
\author Vladimir Florinski
\date 06/26/2024
\param[in] index_in   Index of this exchange site
\param[in] n_parts_in Number of participants
\param[in] buf_size   Size of a single buffer in units of "datatype"
\param[in] ranks_in   List of ranks (with possible repeats)
\param[in] labels_in  List of labels
\param[in] mpi_config MPI configuration object
*/
template <typename datatype>
inline void ExchangeSite<datatype>::SetUpProperties(int index_in, int n_parts_in, int buf_size, const int* ranks_in,
                                                    const int* labels_in, const std::shared_ptr<MPI_Config> mpi_config)
{
   int part, rank, newrank;
   std::vector<int> unique_ranks;
   std::vector<std::vector<int>> part_lists;
   std::vector<int>::iterator it;

   site_index = index_in;
   n_parts = n_parts_in;

// Generate a list of unique ranks
   for (part = 0; part < n_parts; part++) {
      newrank = ranks_in[part];

// Generate the "parts_per_rank" table
      it = std::find(unique_ranks.begin(), unique_ranks.end(), newrank);
      if (it == unique_ranks.end()) {
         unique_ranks.push_back(newrank);

// This rank was not encountered before, so we create a new element in the part list vector for it
         part_lists.emplace_back();
         part_lists.back().push_back(part);
      }
      else {
         part_lists[it - unique_ranks.begin()].push_back(part);
      };
   };

// Create the site communicator. Every process must call "MPI_Comm_create", even if not in the group, so the code will not deadlock.
   MPI_Group parent_group, site_group;
   MPI_Comm_group(mpi_config->glob_comm, &parent_group);
   MPI_Group_incl(parent_group, unique_ranks.size(), unique_ranks.data(), &site_group);
   MPI_Comm_create(mpi_config->glob_comm, site_group, &site_comm);
   MPI_Group_free(&site_group);
   MPI_Group_free(&parent_group);

// Not in the communicator - do not allocate storage
   if (site_comm == MPI_COMM_NULL) return;

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Create the MPI datatype
   MPI_Type_contiguous(sizeof(datatype), MPI_BYTE, &mpi_datatype);
   MPI_Type_commit(&mpi_datatype);

// Allocate the shared buffer
   buffer = new datatype[buf_size * n_parts];

// Calculate the counts and displacements for MPI
   MPI_Comm_size(site_comm, &site_comm_size);
   sendcounts = new int[site_comm_size];
   sdispls = new int[site_comm_size];
   for (rank = 0; rank < site_comm_size; rank++) {
      sendcounts[rank] = part_lists[rank].size() * buf_size;
      sdispls[rank] = (rank == 0 ? 0 : sdispls[rank - 1]) + sendcounts[rank - 1];
   };

// Create map to pointers in buffer space
   int ppr, idx = 0;
   for (rank = 0; rank < site_comm_size; rank++) {
      for (ppr = 0; ppr < part_lists[rank].size(); ppr++) {
         buffer_entry.insert(std::make_pair(part_lists[rank][ppr], &buffer[buf_size * idx]));
         idx++;
      };
   };

// Create a map to find a part based on the label
   for (part = 0; part < n_parts; part++) {
      part_lookup.insert(std::make_pair(labels_in[part], part));
   };
};

/*!
\author Vladimir Florinski
\date 06/26/2024
*/
template <typename datatype>
inline ExchangeSite<datatype>::~ExchangeSite()
{
// Not safe to free a null MPI object
   if (site_comm != MPI_COMM_NULL) MPI_Comm_free(&site_comm);
   if (mpi_datatype != MPI_DATATYPE_NULL) MPI_Type_free(&mpi_datatype);

// Always safe to delete a nullptr
   delete[] sendcounts;
   delete[] sdispls;
   delete[] buffer;
};

/*!
\author Vladimir Florinski
\date 06/26/2024
*/
template <typename datatype>
inline void ExchangeSite<datatype>::Exchange(void)
{
// If this process is the only one participating, the data is accessible directly from the buffer.
   if ((site_comm != MPI_COMM_NULL) && (site_comm_size > 1)) {
      MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buffer, sendcounts, sdispls, mpi_datatype, site_comm);
   };
};

#ifdef GEO_DEBUG

/*!
\brief Test the functionality of this class
\author Vladimir Florinski
\date 06/26/2024
\param[in] mpi_config MPI configuraton object
\param[in] test_rank  Process rank that will do the testing
*/
inline void TestExchange(const std::shared_ptr<MPI_Config> mpi_config, int test_rank)
{
   const int n_parts = 8;
   const int buf_size = 1;
   int part;
   int ranks[n_parts];
   int labels[n_parts];

   if ((test_rank < 0) || (test_rank >= mpi_config->glob_comm_size)) return;

   if (test_rank == mpi_config->glob_comm_rank) {
      std::cerr << "Testing exchange for rank " << test_rank << " (exchage has " << n_parts << " parts)\n";
   };

// Generate a test topology (some processes will not be in it)
   srand48(382667);
   for (part = 0; part < n_parts; part++) {

// Pick ranks at random - this would typically create a few duplicate entries, which is good for testing.
      ranks[part] = drand48() * mpi_config->glob_comm_size;

// Make labels complimentary to parts
      labels[part] = n_parts - part - 1;
      if (test_rank == mpi_config->glob_comm_rank) {
         std::cerr << "Part" << std::setw(3) << part << "   Rank" << std::setw(3) << ranks[part] << "   Label" << std::setw(3) << labels[part] << std::endl;
      };
   };
   if (test_rank == mpi_config->glob_comm_rank) std::cerr << std::endl;

// Create and set up an instance of ExchangeSite
   ExchangeSite<int> exch_site;
   exch_site.SetUpProperties(0, n_parts, buf_size, ranks, labels, mpi_config);

// Fill the buffers with test data - can only do that on the ranks that own the respective part
   if (exch_site.site_comm != MPI_COMM_NULL) {
      for (part = 0; part < n_parts; part++) {
         if (mpi_config->glob_comm_rank == ranks[part]) *exch_site.buffer_entry[part] = part;
         else *exch_site.buffer_entry[part] = -1;
      };
   }
   else if (test_rank == mpi_config->glob_comm_rank) {
      std::cerr << "Process not in the site communicator\n";
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Test part lookup from label
   if ((test_rank == mpi_config->glob_comm_rank) && (exch_site.site_comm != MPI_COMM_NULL)) {
      std::cerr << "Part lookup test for this rank:\n";
      for (part = 0; part < exch_site.n_parts; part++) {
         if (mpi_config->glob_comm_rank == ranks[part]) {
            std::cerr << "Label" << std::setw(3) << labels[part] << "   Part" << std::setw(3) << exch_site.part_lookup[labels[part]] << std::endl;
         };
      };
      std::cerr << std::endl;
   };

// Print buffers before the exchange
   if ((test_rank == mpi_config->glob_comm_rank) && (exch_site.site_comm != MPI_COMM_NULL)) {
      std::cerr << "Buffer before exchange:\n";
      for (part = 0; part < exch_site.n_parts; part++) {
         std::cerr << "Part" << std::setw(3) << part << "   Value" << std::setw(3) << *exch_site.buffer_entry[part] << std::endl;
      };
      std::cerr << std::endl;
   };      

// To properly test we must call this on all processes, not just those in "site_comm".
   exch_site.Exchange();

// Print buffers after the exchange
   if ((test_rank == mpi_config->glob_comm_rank) && (exch_site.site_comm != MPI_COMM_NULL)) {
      std::cerr << "Buffer after exchange:\n";
      for (part = 0; part < exch_site.n_parts; part++) {
         std::cerr << "Part" << std::setw(3) << part << "   Value" << std::setw(3) << *exch_site.buffer_entry[part] << std::endl;
      };
      std::cerr << std::endl;
   };      
};

#endif

};

#endif

#endif
