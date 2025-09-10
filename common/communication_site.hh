/*!
\file communication_site.hh
\brief Implements a generic MPI communication object
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COMMUNICATION_SITE_HH
#define SPECTRUM_COMMUNICATION_SITE_HH

#include "common/mpi_config.hh"

#ifdef USE_MPI

#include <cstring>
#include <map>
#include <vector>
#include <algorithm>

#ifdef GEO_DEBUG
#include "common/print_warn.hh"
#endif

namespace Spectrum {

// Example (n_parts=4):
//          -------------------------               
//          |  part 2   |  part 3   |     ranks_in:   11,  10,  11,  15
//          |           |           |     labels_in: 381, 167, 254, 118
//          | rank: 11  | rank: 15  |
//          | bidx: 254 | bidx: 118 |     site_comm ranks: 11->0, 10->1, 15->2
//          ------------+------------
//          |  part 0   |  part 1   |     sendcounts: 2, 1, 1 (x buf_size)
//          |           |           |     sdispls:    0, 2, 3 (x buf_size)
//          | rank: 11  | rank: 10  |
//          | bidx: 381 | bidx: 167 |     buffer_entry: 0->0, 1->2, 2->1, 3->3
//          -------------------------                      ^           ^
//                                                         |           |
//                                                         ------------- contiguous

/*!
\brief A class representing a general MPI communication site
\author Vladimir Florinski

An object of this class contains an MPI communicator, a data buffer, and metadata to enable data exchange or reduction. The array "buffer" has "n_part" slots corresponding to each participant. Because a process could host multiple participants, "ExchangeSite" will aggregate those into larger slots in the array. The intended use is for the calling function to create an array of derived class objects, one per site. This is redundant because each process may only need access to a subset of all sites. However, MPI communicator creation routines must be called on all processes, and it is desirable to perform all such operations within the class. The unused sites will remain, but since the communicator and buffers are not allocated, the memory lost will be small.
*/
template <typename datatype>
class CommunicationSite
{
protected:

//! Index of this site
   int site_index = -1;

//! Number of participants (actual)
   int n_parts = 0;

//! Buffer size in units of "datatype". It is mainly stored for verification purposes (the class doesn't need it).
   int buf_size = 0;

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

//! Default constructor - protected (this is a base class)
   CommunicationSite(void) {};

public:

//! Copy constructor - deleted because each site must be unique due to MPI restrictions
   CommunicationSite(const CommunicationSite& other) = delete;

//! Move constructor
   CommunicationSite(CommunicationSite&& other);

//! Destructor
   ~CommunicationSite();

//! Assignment operator - deleted because we cannot have multiple copies of the MPI objects
   CommunicationSite& operator =(const CommunicationSite& other) = delete;

//! Retrieve the site index
   int GetIndex(void) const {return n_parts;};

//! Retrieve the part number
   int GetPartCount(void) const {return n_parts;};

//! Retrieve the size of the communicator
   int GetCommSize(void) const {return site_comm_size;};

//! Return the buffer entry point for this part
   datatype* BufferAddress(int part) {return buffer_entry[part];};

//! Return the part for this label
   int PartOfLabel(int label) {return part_lookup[label];};

//! Import lists of ranks and part labels
   void SetUpProperties(int index_in, int n_parts_in, int buf_size_in, const int* ranks_in, const int* labels_in);

#ifdef GEO_DEBUG
//! Print the information about this site from the point of view of one rank
   void PrintSiteInfo(int test_rank) const;
#endif

};

//! Dummy specialization for void datatype
template <> class CommunicationSite<void> {};

/*!
\author Vladimir Florinski
\date 01/16/2025
\param[in] other Object to move into this
*/
template <typename datatype>
inline CommunicationSite<datatype>::CommunicationSite(CommunicationSite<datatype>&& other)
{
   site_index = other.site_index;
   n_parts = other.n_parts;
   buf_size = other.buf_size;
   other.site_index = -1;
   other.n_parts = 0;
   other.buf_size = 0;

// Copy the communicator pointer and set the pointer in "other" to a null comm
   site_comm = other.site_comm;
   site_comm_size = other.site_comm_size;
   other.site_comm = MPI_COMM_NULL;
   other.site_comm_size = 0;

// Copy the data type pointer and set the pointer in "other" to a null type
   mpi_datatype = other.mpi_datatype;
   other.mpi_datatype = MPI_DATATYPE_NULL;

// Move the data buffer
   buffer = other.buffer;
   other.buffer = nullptr;

// Move the count and displacement arrays
   sendcounts = other.sendcounts;
   sdispls = other.sdispls;
   other.sendcounts = nullptr;
   other.sdispls = nullptr;

// Transfer the maps
   part_lookup = std::move(other.part_lookup);
   buffer_entry = std::move(other.buffer_entry);
};

/*!
\author Vladimir Florinski
\date 01/16/2025
\param[in] index_in    Index of this exchange site
\param[in] n_parts_in  Number of participants
\param[in] buf_size_in Size of a single buffer in units of "datatype"
\param[in] ranks_in    List of ranks (with possible repeats)
\param[in] labels_in   List of labels
*/
template <typename datatype>
inline void CommunicationSite<datatype>::SetUpProperties(int index_in, int n_parts_in, int buf_size_in, const int* ranks_in, const int* labels_in)
{
   int newrank;
   std::vector<int> unique_ranks;
   std::vector<std::vector<int>> part_lists;

   site_index = index_in;
   n_parts = n_parts_in;
   buf_size = buf_size_in;

// Generate a list of unique ranks
   for (auto part = 0; part < n_parts; part++) {
      newrank = ranks_in[part];

// Generate the "parts_per_rank" table
      auto it = std::find(unique_ranks.begin(), unique_ranks.end(), newrank);
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
   MPI_Comm_group(MPI_Config::glob_comm, &parent_group);
   MPI_Group_incl(parent_group, unique_ranks.size(), unique_ranks.data(), &site_group);
   MPI_Comm_create(MPI_Config::glob_comm, site_group, &site_comm);
   MPI_Group_free(&site_group);
   MPI_Group_free(&parent_group);

// Not in the communicator - do not allocate storage
   if (site_comm == MPI_COMM_NULL) return;

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Create the MPI datatype
   MPI_Type_contiguous(sizeof(datatype), MPI_BYTE, &mpi_datatype);
   MPI_Type_commit(&mpi_datatype);

// Allocate the shared buffer
   buffer = new datatype[buf_size_in * n_parts];

// Calculate the counts and displacements for MPI
   MPI_Comm_size(site_comm, &site_comm_size);
   sendcounts = new int[site_comm_size];
   sdispls = new int[site_comm_size];
   for (auto rank = 0; rank < site_comm_size; rank++) {
      sendcounts[rank] = part_lists[rank].size() * buf_size_in;
      sdispls[rank] = (rank == 0 ? 0 : sdispls[rank - 1]) + sendcounts[rank - 1];
   };

// Create map to pointers in buffer space
   int idx = 0;
   for (auto rank = 0; rank < site_comm_size; rank++) {
      for (auto ppr = 0; ppr < part_lists[rank].size(); ppr++) {
         buffer_entry.insert(std::make_pair(part_lists[rank][ppr], &buffer[buf_size_in * idx]));
         idx++;
      };
   };

// Create a map to find a part based on the label
   for (auto part = 0; part < n_parts; part++) {
      part_lookup.insert(std::make_pair(labels_in[part], part));
   };
};

/*!
\author Vladimir Florinski
\date 01/16/2025
*/
template <typename datatype>
inline CommunicationSite<datatype>::~CommunicationSite()
{
// Not safe to free a null MPI object
   if (site_comm != MPI_COMM_NULL) MPI_Comm_free(&site_comm);
   if (mpi_datatype != MPI_DATATYPE_NULL) MPI_Type_free(&mpi_datatype);

// Always safe to delete a nullptr
   delete[] sendcounts;
   delete[] sdispls;
   delete[] buffer;
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 01/16/2025
*/
template <typename datatype>
void CommunicationSite<datatype>::PrintSiteInfo(int test_rank) const
{
// Either the communicator wasn't set up or this process is not included
   if (site_comm == MPI_COMM_NULL) return;
   if (test_rank != MPI_Config::glob_comm_rank) return;

   int datatype_size, n_parts_rank, part, label, my_site_rank;
   datatype* buf_ptr;
   MPI_Type_size(mpi_datatype, &datatype_size);
   MPI_Comm_rank(site_comm, &my_site_rank);

   std::cerr << "Printing information for site number " << site_index << std::endl;
   std::cerr << "Site communicator size: " << site_comm_size << std::endl;
   std::cerr << "Datatype size, in bytes: " << datatype_size << std::endl;
   std::cerr << "Buffer size per part, in datatype units: " << buf_size << std::endl;
   std::cerr << "My rank in the site communicator: " << my_site_rank << std::endl;

   for (auto rank = 0; rank < site_comm_size; rank++) {
      std::cerr << "Site rank" << std::setw(3) << rank << ":\n";

// This rank has this many local parts
      n_parts_rank = sendcounts[rank] / buf_size;
      for (auto part_rank = 0; part_rank < n_parts_rank; part_rank++) {
         std::cerr << "          Part";
         buf_ptr = buffer + sdispls[rank] + part_rank * buf_size;

// Find the key (part) for this address. A key not found indicates an error.
         auto it1 = std::find_if(buffer_entry.begin(), buffer_entry.end(),
                                 [buf_ptr](const std::pair<const int, datatype*>& p)
                                 {
                                    return p.second == buf_ptr;
                                 });
         if(it1 == buffer_entry.end()) std::cerr << "  -   Label  -\n";

// Address was found, proceed to find the label
         else {
            part = it1->first;
            std::cerr << std::setw(3) << part << "   Label";

// Find the key (label) for this part
            auto it2 = std::find_if(part_lookup.begin(), part_lookup.end(),
                                    [part](const std::pair<const int, int>& p)
                                    {
                                       return p.second == part;
                                    });

// Print the label if the part was found
            if(it2 == part_lookup.end()) std::cerr << "  -\n";
            else {
               label = it2->first;
               std::cerr << std::setw(3) << label << std::endl;
            };
         };
      };
   };
};

#endif

};

#endif

#endif
