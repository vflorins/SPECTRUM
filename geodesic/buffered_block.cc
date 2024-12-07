/*!
\file buffered_block.cc
\brief Implements the buffered block class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "geodesic/buffered_block.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BufferedBlock public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/28/2024
\param[in] other Object to initialize from
*/
template <int verts_per_face>
BufferedBlock<verts_per_face>::BufferedBlock(const BufferedBlock& other)
                             : StenciledBlock<verts_per_face>(static_cast<StenciledBlock<verts_per_face>>(other))
{
   if (other.side_length != -1) {
      SetDimensions(other.side_length, other.ghost_width, other.n_shells, other.ghost_height, true);
   };
};

/*!
\author Vladimir Florinski
\date 06/28/2024
*/
template <int verts_per_face>
BufferedBlock<verts_per_face>::~BufferedBlock()
{
   FreeStorage();
};

/*!
\author Vladimir Florinski
\date 05/29/2024
\param[in] width     Length of the side, without ghost cells
\param[in] wghost    Width of the ghost cell layer outside the sector
\param[in] height    Hight of the block, without ghost cells
\param[in] hghost    Number of ghost shells outside the slab
\param[in] construct Set to true when called from a constructor
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::SetDimensions(int height, int width, int hghost, int wghost, bool construct)
{
   int site, part;
   NeighborType ntype;

// Call base method
   if (!construct) StenciledBlock<verts_per_face>::SetDimensions(width, wghost, height, hghost, false);

// Free up storage (not that of the base class) because this could be a repeat call
   FreeStorage();

   ComputeRotationMatrices();

// Note: it is possible to have only one participant at a t-face exchange site that is on the simulation boundary, which in principle requires no exchange. We still count it here, and it is up to ExchangeSite to figure out how to handle this efficiently.
   ExchangeSiteCount(verts_per_face, exch_site_count);

// This is the default, the actual number of parts at some sites could be smaller. For simplicity, the face and shell maps are generated for the default configuration.
   ExchangePartCount(edges_per_vert, default_part_per_site);

// Space for the conserved variables
   cons_vars = Create2D<ConservedVariables>(n_shells_withghost, n_faces_withghost);

// Buffer sizes: TFACE
// TODO Check if "buf_width" is needed and whether it should be multiplied by "square_fill"
   buf_length[GEONBR_TFACE] = side_length - (square_fill + 1) * ghost_width;
   buf_width[GEONBR_TFACE] = buf_length[GEONBR_TFACE];
   buf_area[GEONBR_TFACE] = Sqr(buf_length[GEONBR_TFACE]);
   buf_height[GEONBR_TFACE] = ghost_height;

// Buffer sizes: RFACE
   buf_length[GEONBR_RFACE] = side_length - (3 - square_fill) * ghost_width;
   buf_width[GEONBR_RFACE] = ghost_width;
   buf_area[GEONBR_RFACE] = buf_length[GEONBR_RFACE] * buf_width[GEONBR_RFACE] - (square_fill - 1) * Sqr(ghost_width);
   buf_height[GEONBR_RFACE] = n_shells - 2 * ghost_height;

// Buffer sizes: TEDGE
   buf_length[GEONBR_TEDGE] = buf_length[GEONBR_RFACE];
   buf_width[GEONBR_TEDGE] = buf_width[GEONBR_RFACE];
   buf_area[GEONBR_TEDGE] = buf_area[GEONBR_RFACE];
   buf_height[GEONBR_TEDGE] = buf_height[GEONBR_TFACE];

// Buffer sizes: REDGE
   buf_length[GEONBR_REDGE] = ghost_width;
   buf_width[GEONBR_REDGE] = buf_length[GEONBR_REDGE];
   buf_area[GEONBR_REDGE] = Sqr(buf_length[GEONBR_REDGE]);
   buf_height[GEONBR_REDGE] = buf_height[GEONBR_RFACE];

// Buffer sizes: VERTX
   buf_area[GEONBR_VERTX] = buf_length[GEONBR_REDGE];
   buf_width[GEONBR_VERTX] = buf_width[GEONBR_REDGE];
   buf_area[GEONBR_VERTX] = buf_area[GEONBR_REDGE];
   buf_height[GEONBR_VERTX] = buf_height[GEONBR_TEDGE];

// Compute the volumes and check if some are 0 or negative
   for (ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      buf_volume[ntype] = buf_area[ntype] * buf_height[ntype];
      if (buf_volume[ntype] <= 0) {
         PrintError(__FILE__, __LINE__, "Buffer with zero or negative size encountered", true);
         return;
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Allocate storage for buffer face maps and starting shells. Because of redundancy, pointers are often used instead of allocation.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE: 2 sites with 2 PPS
// 1 unique map (1:4)
   buf_face_translation[GEONBR_TFACE] = new int**[exch_site_count[GEONBR_TFACE]];
   buf_face_translation[GEONBR_TFACE][0] = new int*[default_part_per_site[GEONBR_TFACE]];
   buf_face_translation[GEONBR_TFACE][0][0] = new int[buf_area[GEONBR_TFACE]];
   for (part = 1; part < default_part_per_site[GEONBR_TFACE]; part++) {
      buf_face_translation[GEONBR_TFACE][0][part] = buf_face_translation[GEONBR_TFACE][0][0];
   };
   for (site = 1; site < exch_site_count[GEONBR_TFACE]; site++) {
      buf_face_translation[GEONBR_TFACE][site] = buf_face_translation[GEONBR_TFACE][0];
   };

   buf_base_vert[GEONBR_TFACE] = new int*[exch_site_count[GEONBR_TFACE]];
   buf_base_vert[GEONBR_TFACE][0] = new int[default_part_per_site[GEONBR_TFACE]];
   buf_rotation[GEONBR_TFACE] = new int*[exch_site_count[GEONBR_TFACE]];
   buf_rotation[GEONBR_TFACE][0] = new int[default_part_per_site[GEONBR_TFACE]];
   for (site = 1; site < exch_site_count[GEONBR_TFACE]; site++) {
      buf_base_vert[GEONBR_TFACE][site] = buf_base_vert[GEONBR_TFACE][0];
      buf_rotation[GEONBR_TFACE][site] = buf_rotation[GEONBR_TFACE][0];
   };

// 4 unique starting shells (1:1)
   buf_shell_start[GEONBR_TFACE] = new int*[exch_site_count[GEONBR_TFACE]];
   for (site = 0; site < exch_site_count[GEONBR_TFACE]; site++) {
      buf_shell_start[GEONBR_TFACE][site] = new int[default_part_per_site[GEONBR_TFACE]];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE: 3/4 sites with 2 PPS
// 6/8 unique maps (1:1)
   buf_face_translation[GEONBR_RFACE] = new int**[exch_site_count[GEONBR_RFACE]];
   for (site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_face_translation[GEONBR_RFACE][site] = new int*[default_part_per_site[GEONBR_RFACE]];
      for (part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         buf_face_translation[GEONBR_RFACE][site][part] = new int[buf_area[GEONBR_RFACE]];
      };
   };

   buf_base_vert[GEONBR_RFACE] = new int*[exch_site_count[GEONBR_RFACE]];
   buf_rotation[GEONBR_RFACE] = new int*[exch_site_count[GEONBR_RFACE]];
   for (site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_base_vert[GEONBR_RFACE][site] = new int[default_part_per_site[GEONBR_RFACE]];
      buf_rotation[GEONBR_RFACE][site] = new int[default_part_per_site[GEONBR_RFACE]];
   };

// 1 unique starting shell (1:6/1:4)
   buf_shell_start[GEONBR_RFACE] = new int*[exch_site_count[GEONBR_RFACE]];
   buf_shell_start[GEONBR_RFACE][0] = new int[default_part_per_site[GEONBR_RFACE]];
   for (site = 1; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_shell_start[GEONBR_RFACE][site] = buf_shell_start[GEONBR_RFACE][0];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TEDGE: 6/8 sites with 4 PPS
// 6/8 unique maps (1:4)
   buf_face_translation[GEONBR_TEDGE] = new int**[exch_site_count[GEONBR_TEDGE]];
   for (site = 0; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      buf_face_translation[GEONBR_TEDGE][site] = new int*[default_part_per_site[GEONBR_TEDGE]];
      buf_face_translation[GEONBR_TEDGE][site + exch_site_count[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_TEDGE][site];
      for (part = 0; part < exch_part_count[GEONBR_TEDGE] / 2; part++) {
         buf_face_translation[GEONBR_TEDGE][site][part] = new int[buf_area[GEONBR_TEDGE]];
         buf_face_translation[GEONBR_TEDGE][site][part + exch_part_count[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_TEDGE][site][part];
      };
   };



   buf_base_vert[GEONBR_TEDGE] = new int*[exch_site_count[GEONBR_TEDGE]];
   buf_rotationq[GEONBR_TEDGE] = new int*[exch_site_count[GEONBR_TEDGE]];
   for (site = 0; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      buf_face_translation[GEONBR_TEDGE][site] = new int*[default_part_per_site[GEONBR_TEDGE]];
      buf_face_translation[GEONBR_TEDGE][site + exch_site_count[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_TEDGE][site];


      for (part = 0; part < exch_part_count[GEONBR_TEDGE] / 2; part++) {
         buf_face_translation[GEONBR_TEDGE][site][part] = new int[buf_area[GEONBR_TEDGE]];
         buf_face_translation[GEONBR_TEDGE][site][part + exch_part_count[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_TEDGE][site][part];
      };
   };





// 4 unique starting shells (1:6/1:4)
   buf_shell_start[GEONBR_TEDGE] = new int*[exch_site_count[GEONBR_TEDGE]];
   buf_shell_start[GEONBR_TEDGE][0] = new int[default_part_per_site[GEONBR_TEDGE]];
   buf_shell_start[GEONBR_TEDGE][exch_site_count[GEONBR_TEDGE] / 2] = new int[default_part_per_site[GEONBR_TEDGE]];
   for (site = 1; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      buf_shell_start[GEONBR_TEDGE][site] = buf_shell_start[GEONBR_TEDGE][0];
      buf_shell_start[GEONBR_TEDGE][site + exch_site_count[GEONBR_TEDGE] / 2] = buf_shell_start[GEONBR_TEDGE][exch_site_count[GEONBR_TEDGE] / 2];
   };
   
//----------------------------------------------------------------------------------------------------------------------------------------------------

// REDGE: 3/4 sites with 6/4 PPS
// 18/16 unique maps (1:1)
   buf_face_translation[GEONBR_REDGE] = new int**[exch_site_count[GEONBR_REDGE]];
   for (site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
      buf_face_translation[GEONBR_REDGE][site] = new int*[default_part_per_site[GEONBR_REDGE]];
      for (part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
         buf_face_translation[GEONBR_REDGE][site][part] = new int[buf_area[GEONBR_REDGE]];
      };
   };

// 1 unique starting shell (1:18/1:16)
   buf_shell_start[GEONBR_REDGE] = new int*[exch_site_count[GEONBR_REDGE]];
   buf_shell_start[GEONBR_REDGE][0] = new int[default_part_per_site[GEONBR_REDGE]];
   for (site = 1; site < exch_site_count[GEONBR_REDGE]; site++) {
      buf_shell_start[GEONBR_REDGE][site] = buf_shell_start[GEONBR_REDGE][0];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// VERTX: 6/8 sites with 12/8 PPS
// 18/16 unique maps (1:4)
   buf_face_translation[GEONBR_VERTX] = new int**[exch_site_count[GEONBR_VERTX]];
   for (site = 0; site < exch_site_count[GEONBR_VERTX] / 2; site++) {
      buf_face_translation[GEONBR_VERTX][site] = new int*[default_part_per_site[GEONBR_VERTX]];
      buf_face_translation[GEONBR_VERTX][site + exch_site_count[GEONBR_VERTX] / 2] = buf_face_translation[GEONBR_VERTX][site];
      for (part = 0; part < exch_part_count[GEONBR_VERTX] / 2; part++) {
         buf_face_translation[GEONBR_VERTX][site][part] = new int[buf_area[GEONBR_VERTX]];
         buf_face_translation[GEONBR_VERTX][site][part + exch_part_count[GEONBR_VERTX] / 2] = buf_face_translation[GEONBR_VERTX][site][part];
      };
   };

// 4 unique starting shells (1:18/1:16)
   buf_shell_start[GEONBR_VERTX] = new int*[exch_site_count[GEONBR_VERTX]];
   buf_shell_start[GEONBR_VERTX][0] = new int[default_part_per_site[GEONBR_VERTX]];
   buf_shell_start[GEONBR_VERTX][exch_site_count[GEONBR_VERTX] / 2] = new int[default_part_per_site[GEONBR_VERTX]];
   for(site = 1; site < exch_site_count[GEONBR_VERTX] / 2; site++) {
      buf_shell_start[GEONBR_VERTX][site] = buf_shell_start[GEONBR_VERTX][0];
      buf_shell_start[GEONBR_VERTX][site + exch_site_count[GEONBR_VERTX] / 2] = buf_shell_start[GEONBR_VERTX][exch_site_count[GEONBR_VERTX] / 2];
   };
};

/*!
\author Vladimir Florinski
\date 06/18/2024
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::FreeStorage(void)
{
   int site, part;
   NeighborType ntype;

// Free up storage for the variables
   Delete2D(cons_vars);

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE
   delete[] buf_face_translation[GEONBR_TFACE][0][0];
   delete[] buf_face_translation[GEONBR_TFACE][0];
   delete[] buf_face_translation[GEONBR_TFACE];

   for (site = 0; site < exch_site_count[GEONBR_TFACE]; site++) {
      delete[] buf_shell_start[GEONBR_TFACE][site];
   };
   delete[] buf_shell_start[GEONBR_TFACE];

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE
   for (site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      for (part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         delete[] buf_face_translation[GEONBR_RFACE][site][part];
      };
      delete[] buf_face_translation[GEONBR_RFACE][site];
   };
   delete[] buf_face_translation[GEONBR_RFACE];

   delete[] buf_shell_start[GEONBR_RFACE][0];
   delete[] buf_shell_start[GEONBR_RFACE];

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TEDGE
   for (site = 0; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      for (part = 0; part < exch_part_count[GEONBR_TEDGE] / 2; part++) {
         delete[] buf_face_translation[GEONBR_TEDGE][site][part];
      };
      delete[] buf_face_translation[GEONBR_TEDGE][site];
   };
   delete[] buf_face_translation[GEONBR_TEDGE];

   delete[] buf_shell_start[GEONBR_TEDGE][0];
   delete[] buf_shell_start[GEONBR_TEDGE][exch_site_count[GEONBR_TEDGE] / 2];
   delete[] buf_shell_start[GEONBR_TEDGE];

//----------------------------------------------------------------------------------------------------------------------------------------------------

// REDGE
   for (site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
      for (part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
         delete[] buf_face_translation[GEONBR_REDGE][site][part];
      };
      delete[] buf_face_translation[GEONBR_REDGE][site];
   };
   delete[] buf_face_translation[GEONBR_REDGE];

   delete[] buf_shell_start[GEONBR_REDGE][0];
   delete[] buf_shell_start[GEONBR_REDGE];

//----------------------------------------------------------------------------------------------------------------------------------------------------

// VERTX
   for (site = 0; site < exch_site_count[GEONBR_VERTX] / 2; site++) {
      for (part = 0; part < exch_part_count[GEONBR_VERTX] / 2; part++) {
         delete[] buf_face_translation[GEONBR_VERTX][site][part];
      };
      delete[] buf_face_translation[GEONBR_VERTX][site];
   };
   delete[] buf_face_translation[GEONBR_VERTX];

   delete[] buf_shell_start[GEONBR_VERTX][0];
   delete[] buf_shell_start[GEONBR_VERTX][exch_site_count[GEONBR_VERTX] / 2];
   delete[] buf_shell_start[GEONBR_VERTX];

// Call the base class memory deallocator
   StenciledBlock<verts_per_face>::FreeStorage();
};

/*!
\author Vladimir Florinski
\date 06/26/2024
\param[in] ntype         Neighbor type
\param[in] exch_sites_in Array of pointers to exchange sites of this type
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::ImportExchangeSites(NeighborType ntype, std::vector<ExchangeSite<ConservedVariables>*> exch_sites_in)
{
   exch_sites[ntype] = exch_sites_in;
};

/*!
\author Vladimir Florinski
\date 06/27/2024
\param[in] ntype Neighbor type
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::PackBuffers(NeighborType ntype)
{
   int site, my_part, shell, face_buf, face;
   size_t bufidx;
   ConservedVariables* buffer;

   for (site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->part_lookup[block_index];
      buffer = exch_sites[ntype][site]->buffer_entry[my_part];
      bufidx = 0;

// TODO the code below is incorrect
/*
      for (shell = buf_shell_start[ntype][site][my_part]; shell < buf_shell_start[ntype][site][my_part] + buf_height[ntype]; shell++) {
         for (face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
            face = buf_face_translation[ntype][site][my_part][face_buf];
            buffer[bufidx++] = cons_vars[shell][face];
         };
      };
*/

   };
};

/*!
\author Vladimir Florinski
\date 06/27/2024
\param[in] ntype Neighbor type
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::UnPackBuffers(NeighborType ntype)
{
   int site, part, my_part, shell, face_buf, face;
   long int bufidx;
   ConservedVariables* buffer;

   for(site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->part_lookup[block_index];
      for(part = 0; part < default_part_per_site[ntype]; part++) {

// The part equal to ours will be skipped
         if(part == my_part) continue;
         buffer = exch_sites[ntype][site]->buffer_entry[part];
         bufidx = 0;

// TODO the code below is incorrect

/*
// Copy each zone in the buffer to its proper place in the block
         for(shell = buf_shell_start[ntype][site][part]; shell < buf_shell_start[ntype][site][part] + buf_height[ntype]; shell++) {
            for(face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
               face = buf_face_translation[ntype][site][part][face_buf];
               cons_vars[shell][face] = buffers[ntype][site][part][bufidx++];
            };
         };
*/

      };
   };
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 06/18/2024
\param[in] ntype Neighbor type
\param[in] site  Site index of this type
\param[in] part  Participant index in this site
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::PrintBufferMap(NeighborType ntype, int site, int part) const
{
   for(face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
      std::cerr << std::setw(6) << face_buf << std::setw(6) << fbuf_face_translation[ntype][site][part][face_buf] << std::endl;
   };
};

// TODO introduce an exchange testing routine

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BufferedBlock protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2024
*/
template <>
void BufferedBlock<3>::ComputeRotationMatrices(void)
{
// Rotation 0 - yellow
   face_rotation_matrix[0][0][0] =  0; face_rotation_matrix[0][1][0] =  0;
   face_rotation_matrix[0][0][1] =  1; face_rotation_matrix[0][1][1] =  0;
   face_rotation_matrix[0][0][2] =  0; face_rotation_matrix[0][1][2] =  1;
   face_rotation_matrix[0][0][3] =  0; face_rotation_matrix[0][1][3] =  0;
   face_rotation_matrix[0][0][4] =  0; face_rotation_matrix[0][1][4] =  0;

// Rotation 1 - red
   face_rotation_matrix[1][0][0] =  0; face_rotation_matrix[1][1][0] =  1;
   face_rotation_matrix[1][0][1] =  1; face_rotation_matrix[1][1][1] =  2;
   face_rotation_matrix[1][0][2] =  0; face_rotation_matrix[1][1][2] =  0;
   face_rotation_matrix[1][0][3] = -1; face_rotation_matrix[1][1][3] =  0;
   face_rotation_matrix[1][0][4] =  0; face_rotation_matrix[1][1][4] = -1;

// Rotation 2 - green
   face_rotation_matrix[2][0][0] = -1; face_rotation_matrix[2][1][0] =  0;
   face_rotation_matrix[2][0][1] =  0; face_rotation_matrix[2][1][1] =  2;
   face_rotation_matrix[2][0][2] =  0; face_rotation_matrix[2][1][2] = -1;
   face_rotation_matrix[2][0][3] = -1; face_rotation_matrix[2][1][3] =  0;
   face_rotation_matrix[2][0][4] =  0; face_rotation_matrix[2][1][4] =  0;

// Rotation 3 - purple
   face_rotation_matrix[3][0][0] = -1; face_rotation_matrix[3][1][0] = -1;
   face_rotation_matrix[3][0][1] = -1; face_rotation_matrix[3][1][1] =  0;
   face_rotation_matrix[3][0][2] =  0; face_rotation_matrix[3][1][2] = -1;
   face_rotation_matrix[3][0][3] =  0; face_rotation_matrix[3][1][3] =  0;
   face_rotation_matrix[3][0][4] =  0; face_rotation_matrix[3][1][4] =  0;

// Rotation 4 - orange
   face_rotation_matrix[4][0][0] = -1; face_rotation_matrix[4][1][0] = -2;
   face_rotation_matrix[4][0][1] = -1; face_rotation_matrix[4][1][1] = -2;
   face_rotation_matrix[4][0][2] =  0; face_rotation_matrix[4][1][2] =  0;
   face_rotation_matrix[4][0][3] =  1; face_rotation_matrix[4][1][3] =  0;
   face_rotation_matrix[4][0][4] =  0; face_rotation_matrix[4][1][4] =  1;

// Rotation 5 - blue
   face_rotation_matrix[5][0][0] =  0; face_rotation_matrix[5][1][0] = -1;
   face_rotation_matrix[5][0][1] =  0; face_rotation_matrix[5][1][1] = -2;
   face_rotation_matrix[5][0][2] =  0; face_rotation_matrix[5][1][2] =  1;
   face_rotation_matrix[5][0][3] =  1; face_rotation_matrix[5][1][3] =  0;
   face_rotation_matrix[5][0][4] =  0; face_rotation_matrix[5][1][4] =  0;
};

/*!
\author Vladimir Florinski
\date 06/14/2024
*/
template <>
void BufferedBlock<4>::ComputeRotationMatrices(void)
{
// Rotation 0 - yellow
   face_rotation_matrix[0][0][0] =  0; face_rotation_matrix[0][1][0] =  0;
   face_rotation_matrix[0][0][1] =  1; face_rotation_matrix[0][1][1] =  0;
   face_rotation_matrix[0][0][2] =  0; face_rotation_matrix[0][1][2] =  1;
   face_rotation_matrix[0][0][3] =  0; face_rotation_matrix[0][1][3] =  0; // must be 0
   face_rotation_matrix[0][0][4] =  0; face_rotation_matrix[0][1][4] =  0; // must be 0

// Rotation 1 - green
   face_rotation_matrix[1][0][0] = -1; face_rotation_matrix[1][1][0] =  0;
   face_rotation_matrix[1][0][1] =  0; face_rotation_matrix[1][1][1] =  0;
   face_rotation_matrix[1][0][2] =  0; face_rotation_matrix[1][1][2] =  1;
   face_rotation_matrix[1][0][3] =  0; face_rotation_matrix[1][1][3] =  0; // must be 0
   face_rotation_matrix[1][0][4] =  0; face_rotation_matrix[1][1][4] =  0; // must be 0

// Rotation 2 - purple
   face_rotation_matrix[2][0][0] = -1; face_rotation_matrix[2][1][0] = -1;
   face_rotation_matrix[2][0][1] = -1; face_rotation_matrix[2][1][1] =  0;
   face_rotation_matrix[2][0][2] =  0; face_rotation_matrix[2][1][2] =  0;
   face_rotation_matrix[2][0][3] =  0; face_rotation_matrix[2][1][3] =  0; // must be 0
   face_rotation_matrix[2][0][4] =  0; face_rotation_matrix[2][1][4] =  0; // must be 0

// Rotation 3 - blue
   face_rotation_matrix[3][0][0] =  0; face_rotation_matrix[3][1][0] = -1;
   face_rotation_matrix[3][0][1] =  0; face_rotation_matrix[3][1][1] =  0;
   face_rotation_matrix[3][0][2] =  0; face_rotation_matrix[3][1][2] = -1;
   face_rotation_matrix[3][0][3] =  0; face_rotation_matrix[3][1][3] =  0; // must be 0
   face_rotation_matrix[3][0][4] =  0; face_rotation_matrix[3][1][4] =  0; // must be 0
};




template <>
void BufferedBlock<3>::ComputeBaseVertices(void)
{
   






};



/*!
\author Vladimir Florinski
\date 06/18/2024
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::ComputeBufferTranslations(void)
{
   int ibuf, jbuf, site, part, face, bufidx;
   NeighborType ntype;
   std::pair<int, int> ij;
   std::pair<int, int> base_vertex;

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE: 2 sites with 2 PPS
// 1 unique map (1:4)
   bufidx = 0;
   base_vertex = std::make_pair(2 * square_fill * ghost_width, 2 * ghost_width);
   for (ibuf = 0; ibuf < buf_length[GEONBR_TFACE]; ibuf++) {
      for (jbuf = 0; jbuf <= MaxFaceJ(buf_width[GEONBR_TFACE], buf_length[GEONBR_TFACE]); jbuf++) {
         ij = RotateFace(base_vertex, 0, ibuf, jbuf);
         buf_face_translation[GEONBR_TFACE][0][0][bufidx++] = face_index_sector[ij.first][ij.second];
      };
   };   

// 4 unique starting shells (1:1)
   buf_shell_start[GEONBR_TFACE][0][0] = 0;
   buf_shell_start[GEONBR_TFACE][0][1] = ghost_height;
   buf_shell_start[GEONBR_TFACE][1][0] = n_shells_withghost - 2 * ghost_height;
   buf_shell_start[GEONBR_TFACE][1][1] = n_shells_withghost - ghost_height;

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE neighbor face translation - 6/8 maps


   site = 0;
//   part = 
   bufidx = 0;
   
   rot = 




   rot = edges_per_vert / 2;
   for (site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      my_part = exch_sites[GEONBR_RFACE][site]->part_lookup[block_index];
      
      base_vertex = std::make_pair(total_length - (3 - square_fill) * ghost_width, 2 * ghost_width);   
   
   
      for (part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         bufidx = 0;
//         base_vertex = ???;
         for(ibuf = 0; ibuf < buf_length[GEONBR_RFACE]; ibuf++) {
            for(jbuf = 0; jbuf <= MaxFaceJ(buf_length[GEONBR_RFACE], buf_width[GEONBR_RFACE]); jbuf++) {
//               ij = RotateFace(base_vertex, ???, ibuf, jbuf);
               buf_face_translation[GEONBR_RFACE][site][part][bufidx++] = face_index_sector[ij.first][ij.second];
            };
         };   
      };
   };
  

// TODO Write the rest of this function


};

template class BufferedBlock<3>;
template class BufferedBlock<4>;

};
