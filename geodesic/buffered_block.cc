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
   if(!construct) StenciledBlock<verts_per_face>::SetDimensions(width, wghost, height, hghost, false);
   ComputeRotationMatrices();

// Note: it is possible to have only one participant at a t-face exchange site that is on the simulation boundary, which in principle requires no exchange. We still count it here, and it is up to ExchangeSite to figure out how to handle this efficiently.
   ExchangeSiteCount(verts_per_face, exch_site_count);

// This is the default, the actual number of parts at some sites could be smaller. For simplicity, the face and shell maps are generated for the default configuration.
   ExchangePartCount(edges_per_vert, default_part_per_site);

// Space for the conserved variables
   cons_vars = Create2D<ConservedVariables>(n_shells_withghost, n_faces_withghost);

// Buffer sizes: TFACE
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
   for(ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      buf_volume[ntype] = buf_area[ntype] * buf_height[ntype];
      if(buf_volume[ntype] <= 0) {
         PrintError(__FILE__, __LINE__, "Buffer with zero or negative size encountered", true);
         return;
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Allocate storage for buffer face maps. There is much redundancy in these maps, so only the minimum amount of storage is used; the remaining indices are just pointers.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE: 2 sites with 2 PPS, all 4 maps are identical
   buf_face_translation[GEONBR_TFACE] = new int**[exch_site_count[GEONBR_TFACE]];
   buf_face_translation[GEONBR_TFACE][0] = new int*[default_part_per_site[GEONBR_TFACE]];
   buf_face_translation[GEONBR_TFACE][0][0] = new int[buf_area[GEONBR_TFACE]];
   for(part = 1; part < default_part_per_site[GEONBR_TFACE]; part++) {
      buf_face_translation[GEONBR_TFACE][0][part] = buf_face_translation[GEONBR_TFACE][0][0];
   };
   for(site = 1; site < exch_site_count[GEONBR_TFACE]; site++) {
      buf_face_translation[GEONBR_TFACE][site] = buf_face_translation[GEONBR_TFACE][0];
   };

// All 4 starting shells are different
   buf_shell_start[GEONBR_TFACE] = new int*[exch_site_count[GEONBR_TFACE]];
   for(site = 0; site < exch_site_count[GEONBR_TFACE]; site++) {
      buf_shell_start[GEONBR_TFACE][site] = new int[default_part_per_site[GEONBR_TFACE]];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE: 3/4 sites with 2 PPS, all 6/8 maps are different
   buf_face_translation[GEONBR_RFACE] = new int**[exch_site_count[GEONBR_RFACE]];
   for(site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_face_translation[GEONBR_RFACE][site] = new int*[default_part_per_site[GEONBR_RFACE]];
      for(part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         buf_face_translation[GEONBR_RFACE][site][part] = new int[buf_area[GEONBR_RFACE]];
      };
   };

// All 6/8 starting shells are identical; the code computing the shells must assign identical values to three sets
   buf_shell_start[GEONBR_RFACE] = new int*[exch_site_count[GEONBR_RFACE]];
   buf_shell_start[GEONBR_RFACE][0] = new int[default_part_per_site[GEONBR_RFACE]];
   for(site = 1; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_shell_start[GEONBR_RFACE][site] = buf_shell_start[GEONBR_RFACE][0];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TEDGE: 6/8 sites with 4 PPS, the maps are identical to those of RFACE
   buf_face_translation[GEONBR_TEDGE] = new int**[exch_site_count[GEONBR_TEDGE]];
   for(site = 0; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      buf_face_translation[GEONBR_TEDGE][site] = buf_face_translation[GEONBR_RFACE][site];
      buf_face_translation[GEONBR_TEDGE][site + exch_site_count[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_RFACE][site];
   };

// The starting shells are identical to those of TFACE, but the arrays are longer, so some storage must be allocated here
   buf_shell_start[GEONBR_TEDGE] = new int*[exch_site_count[GEONBR_TEDGE]];
   buf_shell_start[GEONBR_TEDGE][0] = new int[default_part_per_site[GEONBR_TEDGE]];
   buf_shell_start[GEONBR_TEDGE][exch_site_count[GEONBR_TEDGE] / 2] = new int[default_part_per_site[GEONBR_TEDGE]];
   for(site = 1; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      buf_shell_start[GEONBR_TEDGE][site] = buf_shell_start[GEONBR_TEDGE][0];
      buf_shell_start[GEONBR_TEDGE][site + exch_site_count[GEONBR_TEDGE] / 2] = buf_shell_start[GEONBR_TEDGE][exch_site_count[GEONBR_TEDGE] / 2];
   };
   
//----------------------------------------------------------------------------------------------------------------------------------------------------

// REDGE: 3/4 sites with 6/4 PPS, all 18/16 maps are different
   buf_face_translation[GEONBR_REDGE] = new int**[exch_site_count[GEONBR_REDGE]];
   for(site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
      buf_face_translation[GEONBR_REDGE][site] = new int*[default_part_per_site[GEONBR_REDGE]];
      for(part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
         buf_face_translation[GEONBR_REDGE][site][part] = new int[buf_area[GEONBR_REDGE]];
      };
   };

// The starting shells are identical to those of RFACE, but the arrays are longer, so some storage must be allocated here
   buf_shell_start[GEONBR_REDGE] = new int*[exch_site_count[GEONBR_REDGE]];
   buf_shell_start[GEONBR_REDGE][0] = new int[default_part_per_site[GEONBR_REDGE]];
   for(site = 1; site < exch_site_count[GEONBR_REDGE]; site++) {
      buf_shell_start[GEONBR_REDGE][site] = buf_shell_start[GEONBR_REDGE][0];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// VERTX: 6/8 sites with 12/8 PPS, the maps are identical to those of REDGE
   buf_face_translation[GEONBR_VERTX] = new int**[exch_site_count[GEONBR_VERTX]];
   for(site = 0; site < exch_site_count[GEONBR_VERTX] / 4; site++) {
      buf_face_translation[GEONBR_VERTX][site] = buf_face_translation[GEONBR_REDGE][site];
      buf_face_translation[GEONBR_VERTX][site + exch_site_count[GEONBR_VERTX] / 4] = buf_face_translation[GEONBR_REDGE][site];
      buf_face_translation[GEONBR_VERTX][site + exch_site_count[GEONBR_VERTX] / 2] = buf_face_translation[GEONBR_REDGE][site];
      buf_face_translation[GEONBR_VERTX][site + 3 * exch_site_count[GEONBR_VERTX] / 4] = buf_face_translation[GEONBR_REDGE][site];
   };

// The starting shells are identical to those of TFACE and TEDGE, but the arrays are longer, so some storage must be allocated here
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
// Deallocate storage for buffer face maps.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE
   delete[] buf_face_translation[GEONBR_TFACE][0][0];
   delete[] buf_face_translation[GEONBR_TFACE][0];
   delete[] buf_face_translation[GEONBR_TFACE];
   for(site = 0; site < exch_site_count[GEONBR_TFACE]; site++) {
      delete[] buf_shell_start[GEONBR_TFACE][site];
   };
   delete[] buf_shell_start[GEONBR_TFACE];

// RFACE
   for(site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      for(part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         delete[] buf_face_translation[GEONBR_RFACE][site][part];
      };
      delete[] buf_face_translation[GEONBR_RFACE][site];
   };
   delete[] buf_face_translation[GEONBR_RFACE];
   delete[] buf_shell_start[GEONBR_RFACE][0];
   delete[] buf_shell_start[GEONBR_RFACE];

// TEDGE
   delete[] buf_face_translation[GEONBR_TEDGE];
   delete[] buf_shell_start[GEONBR_TEDGE][0];
   delete[] buf_shell_start[GEONBR_TEDGE][exch_site_count[GEONBR_TEDGE] / 2];
   delete[] buf_shell_start[GEONBR_TEDGE];

// REDGE
   for(site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
      for(part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
         delete[] buf_face_translation[GEONBR_REDGE][site][part];
      };
      delete[] buf_face_translation[GEONBR_REDGE][site];
   };
   delete[] buf_face_translation[GEONBR_REDGE];
   delete[] buf_shell_start[GEONBR_REDGE][0];
   delete[] buf_shell_start[GEONBR_REDGE];

// VERTX
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
\param ntype  Neighbor type
\param comms  A global array of communicators
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::ImportExchangeSites(NeighborType ntype, std::vector<std::shared_ptr<ExchangeSite<ConservedVariables>>> exch_sites_in)
{
   exch_sites[ntype] = exch_sites_in;
};

/*!
\author Vladimir Florinski
\date 06/27/2024
\param ntype Neighbor type
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::PackBuffers(NeighborType ntype)
{
   int site, my_part, shell, face_buf, face;
   long int bufidx;
   ConservedVariables* buffer;

   for(site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->part_lookup[block_index];
      buffer = exch_sites[ntype][site]->buffer_entry[my_part];
      bufidx = 0;

// TODO the code below is incorrect
/*
      for(shell = buf_shell_start[ntype][site][my_part]; shell < buf_shell_start[ntype][site][my_part] + buf_height[ntype]; shell++) {
         for(face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
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
\param ntype Neighbor type
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

// TFACE neighbor face translation (1)
   bufidx = 0;
   base_vertex = std::make_pair(2 * square_fill * ghost_width, 2 * ghost_width);
   for(ibuf = 0; ibuf < buf_length[GEONBR_TFACE]; ibuf++) {
      for(jbuf = 0; jbuf <= MaxFaceJ(buf_length[GEONBR_TFACE], buf_width[GEONBR_TFACE]); jbuf++) {
         ij = RotateFace(base_vertex, 0, ibuf, jbuf);
         buf_face_translation[GEONBR_TFACE][0][0][bufidx++] = face_index_sector[ij.first][ij.second];
      };
   };   

// RFACE neighbor face translation (6)
   for(site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      for(part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
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

};
