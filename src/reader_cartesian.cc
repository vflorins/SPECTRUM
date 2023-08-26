/*!
\file reader_cartesian.hh
\brief Implements a class of data reader for a uniform Cartesian grid
\author Juan G Alonso Guzman
*/

#include "reader_cartesian.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

namespace Spectrum {

//! Global Cartesian data structure
ReaderCartesian CartData;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ReaderCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[in] data_filename null terminated character array containing the name of the data file
\param[in] fname_len     maximum length of data_filename array (maybe not necessary)
\param[in] verbose       flag to output status messages (1) or not (0)
*/
void ReadCartesianHeader(const char* data_filename, int fname_len, int verbose)
{
   int var;          // physical variable index
   std::string header_filename(data_filename);

// Read contents of header file
   header_filename = header_filename.substr(0,header_filename.size()-4) + ".info";
   std::ifstream header_file (header_filename.c_str());

   if(verbose) {
      std::setprecision(6);
      std::cerr << "Reading Cartesian header file: " << header_filename << std::endl;
   };

// Read number of blocks per dimension
   header_file >> CartData.Nblocks[0] >> CartData.Nblocks[1] >> CartData.Nblocks[2];
   if(verbose) {
      std::cerr << std::setw(8) << "Nx"
                << std::setw(8) << "Ny"
                << std::setw(8) << "Nz"
                << std::endl;
      std::cerr << std::setw(8) << CartData.Nblocks[0]
                << std::setw(8) << CartData.Nblocks[1]
                << std::setw(8) << CartData.Nblocks[2]
                << std::endl;
   };

// Read minimum domain coordinates
   header_file >> CartData.domain_min[0] >> CartData.domain_min[1] >> CartData.domain_min[2];
   if(verbose) {
      std::cerr << std::setw(16) << "DomMin_x"
                << std::setw(16) << "DomMin_y"
                << std::setw(16) << "DomMin_z"
                << std::endl;
      std::cerr << std::setw(16) << CartData.domain_min[0]
                << std::setw(16) << CartData.domain_min[1]
                << std::setw(16) << CartData.domain_min[2]
                << std::endl;
   };

// Read maximum domain coordinates
   header_file >> CartData.domain_max[0] >> CartData.domain_max[1] >> CartData.domain_max[2];
   if(verbose) {
      std::cerr << std::setw(16) << "DomMax_x"
                << std::setw(16) << "DomMax_y"
                << std::setw(16) << "DomMax_z"
                << std::endl;
      std::cerr << std::setw(16) << CartData.domain_max[0]
                << std::setw(16) << CartData.domain_max[1]
                << std::setw(16) << CartData.domain_max[2]
                << std::endl;
   };

// Read number of variables
   header_file >> CartData.Nvar;
   if(verbose) {
      std::cerr << std::setw(8) << "Nvar"
                << std::endl;
      std::cerr << std::setw(8) << CartData.Nvar
                << std::endl;
   };

// TODO: Adapt to work with ghost cells

// Calculate additional quantities and allocate memory
   CartData.data_block_size = CartData.Nvar * block_size_cartesian.Prod();
   CartData.block_length = (CartData.domain_max - CartData.domain_min) / CartData.Nblocks;
   CartData.zone_length = (CartData.block_length) / block_size_cartesian;
   delete[] CartData.variables;
   CartData.variables = new double[CartData.data_block_size * CartData.Nblocks.Prod()];

   header_file.close();
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[in] data_filename null terminated character array containing the name of the data file
\param[in] fname_len     maximum length of data_filename array (maybe not necessary)
\param[in] header        flag to also read header (1) or not (0)
\param[in] verbose       flag to output status messages (1) or not (0)
*/
void ReadCartesianData(const char* data_filename, int fname_len, int read_header, int verbose)
{
   int var;          // physical variable index
   int ib,jb,kb;     // block indices
   int iz,jz,kz;     // zone indices
   long var_idx = 0; // variables array index
   MultiIndex block_idx, zone_idx; // block and zone indices

   if(read_header) ReadCartesianHeader(data_filename, fname_len, verbose);

   std::ifstream data_file (data_filename);

   if(verbose) {
      std::cerr << "Reading Cartesian data file: " << data_filename << std::endl;
      std::cerr << "Total number of variables to read: " << CartData.data_block_size * CartData.Nblocks.Prod();
      std::cerr << "Progress:     ";
   };

   data_file.read((char*)CartData.variables, CartData.data_block_size * CartData.Nblocks.Prod() * sizeof(double));
   if(verbose) std::cerr << "\e[4D100%\n";
   data_file.close();

// TODO: Adapt to work with ghost cells
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
*/
void ReadCartesianClean(void)
{
   delete[] CartData.variables;
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[out] domain_min Minimum domain coordinate
\param[out] domain_max Maximum domain coordinate
*/
void ReadCartesianGetDomain(double* domain_min_out, double* domain_max_out)
{
   memcpy(domain_min_out, CartData.domain_min.Data(), 3 * sizeof(double));
   memcpy(domain_max_out, CartData.domain_max.Data(), 3 * sizeof(double));
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[in]  pos  current position
\param[out] node node ID containing pos
*/
void ReadCartesianGetNode(const double* pos, int* node_id)
{
   MultiIndex node_idx; // Node indices
   
   if(   (pos[0] < CartData.domain_min[0]) || (pos[0] > CartData.domain_max[0])
      || (pos[1] < CartData.domain_min[1]) || (pos[1] > CartData.domain_max[1])
      || (pos[2] < CartData.domain_min[2]) || (pos[2] > CartData.domain_max[2]) ) {
      std::cerr << "ReadCartesianGetNode Error: Position out of domain bounds.\n";
      *node_id = -1;
      return;
   };

   node_idx = (pos - CartData.domain_min) / CartData.block_length;
   *node_id = node_idx.i + CartData.Nblocks.i * (node_idx.j + CartData.Nblocks.j * node_idx.k);
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[in]  node     node ID
\param[out] face_min coordinates of minimum block corner
\param[out] face_max coordinates of maximum block corner
*/
void ReadCartesianGetBlockCorners(int node, double* face_min, double* face_max)
{
   MultiIndex node_idx;

   node_idx.i = node % CartData.Nblocks.i;
   node_idx.j = ( (node - node_idx.i) / CartData.Nblocks.i ) % CartData.Nblocks.j;
   node_idx.k = ( (node - node_idx.i) / CartData.Nblocks.i - node_idx.j ) / CartData.Nblocks.j;

   for(auto xyz = 0; xyz < 3; xyz++) {
      face_min[xyz] = CartData.domain_min[xyz] + CartData.block_length[xyz] * node_idx[xyz];
      face_max[xyz] = face_min[xyz] + CartData.block_length[xyz];
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[in]  node            node ID
\param[out] neighbor_nodes  array with node ID's for all neighboring nodes
\param[out] neighbor_levels array with relative refinement level for all neighboring nodes
*/

//                         -Nx*Ny <- z -> Nx*Ny
//                ------------------------------------->
//     +==============+      +==============+      +==============+
//  Nx |  6 |  7 |  8 |   Nx | 15 | 16 | 17 |   Nx | 24 | 25 | 26 |
//   ^ |----|----|----|    ^ |----|----|----|    ^ |----|----|----|
//   y |  3 |  4 |  5 |    y | 12 | 13 | 14 |    y | 21 | 22 | 23 |
//   ! |----|----|----|    ! |----|----|----|    ! |----|----|----|
// -Nx |  0 |  1 |  2 |  -Nx |  9 | 10 | 11 |  -Nx | 18 | 19 | 20 |
//     +==============+      +==============+      +==============+
//       -1 <- x -> 1          -1 <- x -> 1          -1 <- x -> 1

void ReadCartesianGetNodeNeighbors(int node, int* neighbor_nodes, int* neighbor_levels)
{
   int i,j,k,n = 0;

   for(k = 0; k < max_neighbors_per_dim_cartesian; k++) { // iterate over z
      for(j = 0; j < max_neighbors_per_dim_cartesian; j++) { // iterate over y
         for(i = 0; i < max_neighbors_per_dim_cartesian; i++) { // iterate over x
            neighbor_nodes[n] = node + (i - 1) + CartData.Nblocks.i * ( (j - 1) + CartData.Nblocks.j * (k - 1) ); // neighbor[13] = node
            neighbor_levels[0] = 0;
            n++; // Running index
         };
      };
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
\param[in]  node       node ID
\param[out] block_vars variables in the block
*/
void ReadCartesianGetBlockData(int node, double* block_vars)
{
   memcpy(block_vars, CartData.variables + CartData.data_block_size * node, CartData.data_block_size * sizeof(double));
};

};
