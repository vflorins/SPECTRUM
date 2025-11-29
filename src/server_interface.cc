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
\author Juan G Alonso Guzman
\date 07/28/2023
\param[out] block_ptr pointer to block type
*/
// TODO: look up all calls to this function for possible replacement with a simple "new"
template <typename HConfig>
void ServerInterface<HConfig>::InitializeBlockPtr(Block* &block_ptr)
{
   block_ptr = new Block();
};


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
   Block* block_tmp;
   InitializeBlockPtr(block_tmp);
//   n_variables = block_tmp->GetVariableCount();

// Set up MPI data type for "BlockXXXX" (without the dynamic storage)
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
// Create block datatype
   CreateBlockDatatype();
// Create stencil datatype
   CreateStencilDatatype();
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

};

