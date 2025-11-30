/*!
\file server_interface.cc
\brief server-worker interface for a data-defined grid
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/mpi_config.hh"
#include "server_interface.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerInterface methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/27/2025
*/
template <typename HConfig>
void ServerInterface<HConfig>::CreateInquiryDatatype(void)
{
   MPI_Datatype inquiry_types[] = {MPI_INT, MPI_INT, MPI_DOUBLE};
   int inquiry_lengths[] = {1, 1, 3};
   MPI_Aint inquiry_displ[3];

// Figure out field displacements using "_inquiry" as template
   MPI_Get_address(&_inquiry.type, &inquiry_displ[0]);
   MPI_Get_address(&_inquiry.node, &inquiry_displ[1]);
   MPI_Get_address(&_inquiry.pos , &inquiry_displ[2]);
   for (auto i = 2; i >= 0; i--) inquiry_displ[i] -= inquiry_displ[0];

// Commit the type
   MPI_Type_create_struct(3, inquiry_lengths, inquiry_displ, inquiry_types, &MPIInquiryType);
   MPI_Type_commit(&MPIInquiryType);
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void ServerInterface<HConfig>::CreateBlockDatatype(void)
{
   int i;

// Temporary object for displacement calculations
   auto block_tmp = new Block();

// Set up MPI data type for Block (without the dynamic storage)
   MPI_Datatype block_types[] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
   int block_lengths[] = {1, 3, 3, 3, 3};
   MPI_Aint block_start, block_displ[5];

// Figure out field displacements using "block_tmp" as template
   MPI_Get_address(block_tmp, &block_start);
   MPI_Get_address(block_tmp->GetNodeAddress(), &block_displ[0]);
   MPI_Get_address(block_tmp->GetFaceMinAddress(), &block_displ[1]);
   MPI_Get_address(block_tmp->GetFaceMaxAddress(), &block_displ[2]);
   MPI_Get_address(block_tmp->GetFaceMinPhysAddress(), &block_displ[3]);
   MPI_Get_address(block_tmp->GetFaceMaxPhysAddress(), &block_displ[4]);
   for (i = 4; i >= 0; i--) block_displ[i] -= block_start;

// Commit the type
   MPI_Type_create_struct(5, block_lengths, block_displ, block_types, &MPIBlockType);
   MPI_Type_commit(&MPIBlockType);
   delete block_tmp;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
*/
template <typename HConfig>
void ServerInterface<HConfig>::CreateStencilDatatype(void)
{
   int i;

// Set up MPI data type for "InterpolationStencil"
   MPI_Datatype stencil_types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
   int stencil_lengths[] = {1, 8, 24, 8, 24};
   MPI_Aint stencil_start, stencil_displ[5];

// Figure out field displacements using "stencil" for template
   MPI_Get_address(&stencil, &stencil_start);
   MPI_Get_address(&stencil.n_elements , &stencil_displ[0]);
   MPI_Get_address(&stencil.blocks     , &stencil_displ[1]);
   MPI_Get_address(&stencil.zones      , &stencil_displ[2]);
   MPI_Get_address(&stencil.weights    , &stencil_displ[3]);
   MPI_Get_address(&stencil.derivatives, &stencil_displ[4]);
   for (i = 4; i >= 0; i--) stencil_displ[i] -= stencil_start;

// Commit the type
   MPI_Type_create_struct(5, stencil_lengths, stencil_displ, stencil_types, &MPIStencilType);
   MPI_Type_commit(&MPIStencilType);
};


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 07/19/2023
*/
template <typename HConfig>
void ServerInterface<HConfig>::ServerInterfaceStart(void)
{
// Set up MPI data type for "Inquiry"
   CreateInquiryDatatype();
// Create stencil datatype
   CreateStencilDatatype();
// Create block datatype
   CreateBlockDatatype();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig>
void ServerInterface<HConfig>::ServerInterfaceFinish(void)
{
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);
   MPI_Type_free(&MPIInquiryType);
};



/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Request the bounding dimensions
*/
template <typename HConfig>
void ServerInterface<HConfig>::LoadFromReader(BlockPtr &blockptr) const
{
   // todo review for header/include/linking issues
   if constexpr (HConfig::background == Config::Background::DataCartesian)
      ReadCartesianGetBlockCorners(blockptr->node, blockptr->face_min.Data(), blockptr->face_max.Data());
   else if constexpr (HConfig::background == Config::Background::DataBATL)
      spectrum_get_block_corners(blockptr->node, blockptr->face_min.Data(), blockptr->face_max.Data());
}


/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Request all neighbors
*/
template <typename HConfig>
void ServerInterface<HConfig>::LoadNeighborsFromReader(BlockPtr &blockptr) const
{
   if constexpr (HConfig::background == Config::Background::DataCartesian)
      ReadCartesianGetNodeNeighbors(blockptr->node, blockptr->neighbor_nodes, blockptr->neighbor_levels);
   else if constexpr (HConfig::background == Config::Background::DataBATL)
      spectrum_get_all_neighbor_copies(blockptr->node, blockptr->neighbor_nodes, blockptr->neighbor_levels);
}



/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/26/2025
Load all fields into the block
*/
template <typename HConfig>
void ServerInterface<HConfig>::LoadFieldsFromReader(BlockPtr &blockptr) const
{
   if constexpr (HConfig::background == Config::Background::DataCartesian)
      ReadCartesianGetBlockData(blockptr->node, blockptr->fields[0].Array());
   else if constexpr (HConfig::background == Config::Background::DataBATL)
      spectrum_get_block_data(blockptr->node, blockptr->fields[0].Array());
}




};

