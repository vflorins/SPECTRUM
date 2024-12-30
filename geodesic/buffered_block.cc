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
\date 12/23/2024
\param[in] width     Length of the side, without ghost cells
\param[in] wghost    Width of the ghost cell layer outside the sector
\param[in] height    Hight of the block, without ghost cells
\param[in] hghost    Number of ghost shells outside the slab
\param[in] construct Set to true when called from a constructor
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::SetDimensions(int width, int wghost, int height, int hghost, bool construct)
{
// Call base method
   if (!construct) StenciledBlock<verts_per_face>::SetDimensions(width, wghost, height, hghost, false);

// Free up storage (not that of the base class) because this could be a repeat call
   FreeStorage();

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
   buf_length[GEONBR_VERTX] = buf_length[GEONBR_REDGE];
   buf_width[GEONBR_VERTX] = buf_width[GEONBR_REDGE];
   buf_area[GEONBR_VERTX] = buf_area[GEONBR_REDGE];
   buf_height[GEONBR_VERTX] = buf_height[GEONBR_TEDGE];

// Compute the volumes and check if some are 0 or negative. Note that the area is _not_ always equal to buf_length*buf_width because the width is always measured in vertices, not faces.
   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
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
   for (auto site = 0; site < exch_site_count[GEONBR_TFACE] / 2; site++) {
      buf_face_translation[GEONBR_TFACE][site] = new int*[default_part_per_site[GEONBR_TFACE]];
      buf_face_translation[GEONBR_TFACE][site + exch_site_count[GEONBR_TFACE] / 2] = buf_face_translation[GEONBR_TFACE][site];
      for (auto part = 0; part < default_part_per_site[GEONBR_TFACE] / 2; part++) {
         buf_face_translation[GEONBR_TFACE][site][part] = new int[buf_area[GEONBR_TFACE]];
         buf_face_translation[GEONBR_TFACE][site][part + default_part_per_site[GEONBR_TFACE] / 2] = buf_face_translation[GEONBR_TFACE][site][part];
      };
   };

// 4 unique starting shells (1:1)
   buf_shell_start[GEONBR_TFACE] = new int**[exch_site_count[GEONBR_TFACE]];
   for (auto site = 0; site < exch_site_count[GEONBR_TFACE]; site++) {
      buf_shell_start[GEONBR_TFACE][site] = new int*[default_part_per_site[GEONBR_TFACE]];
      for (auto part = 0; part < default_part_per_site[GEONBR_TFACE]; part++) {
         buf_shell_start[GEONBR_TFACE][site][part] = buf_shell_tab + site * exch_site_count[GEONBR_TFACE] + part;
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE: 3/4 sites with 2 PPS
// 6/8 unique maps (1:1)
   buf_face_translation[GEONBR_RFACE] = new int**[exch_site_count[GEONBR_RFACE]];
   for (auto site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_face_translation[GEONBR_RFACE][site] = new int*[default_part_per_site[GEONBR_RFACE]];
      for (auto part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         buf_face_translation[GEONBR_RFACE][site][part] = new int[buf_area[GEONBR_RFACE]];
      };
   };

// 1 unique starting shell (1:6/1:8)
   buf_shell_start[GEONBR_RFACE] = new int**[exch_site_count[GEONBR_RFACE]];
   buf_shell_start[GEONBR_RFACE][0] = new int*[default_part_per_site[GEONBR_RFACE]];
   for (auto part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
      buf_shell_start[GEONBR_RFACE][0][part] = buf_shell_tab + 1;
   };
   for (auto site = 1; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_shell_start[GEONBR_RFACE][site] = buf_shell_start[GEONBR_RFACE][0];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TEDGE: 6/8 sites with 4 PPS
// 6/8 unique maps (1:4), storage points to RFACE
   buf_face_translation[GEONBR_TEDGE] = new int**[exch_site_count[GEONBR_TEDGE]];
   for (auto site = 0; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
      buf_face_translation[GEONBR_TEDGE][site] = new int*[default_part_per_site[GEONBR_TEDGE]];
      buf_face_translation[GEONBR_TEDGE][site + exch_site_count[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_TEDGE][site];
      for (auto part = 0; part < default_part_per_site[GEONBR_TEDGE] / 2; part++) {
         buf_face_translation[GEONBR_TEDGE][site][part] = buf_face_translation[GEONBR_RFACE][site][part];
         buf_face_translation[GEONBR_TEDGE][site][part + default_part_per_site[GEONBR_TEDGE] / 2] = buf_face_translation[GEONBR_TEDGE][site][part];
      };
   };

// 4 unique starting shells (1:6/1:8)
   buf_shell_start[GEONBR_TEDGE] = new int**[exch_site_count[GEONBR_TEDGE]];
   for (auto site = 0; site < exch_site_count[GEONBR_TEDGE] / verts_per_face; site++) {
      buf_shell_start[GEONBR_TEDGE][site * verts_per_face] = new int*[default_part_per_site[GEONBR_TEDGE]];
      for (auto part = 0; part < 2; part++) {
         for (auto part1 = 0; part1 < default_part_per_site[GEONBR_TEDGE] / 2; part1++) {
            buf_shell_start[GEONBR_TEDGE][site * verts_per_face][2 * part + part1] = buf_shell_tab + 2 * site + part;
         };
      };
      for (auto site1 = 1; site1 < exch_site_count[GEONBR_TEDGE] / 2; site1++) {
         buf_shell_start[GEONBR_TEDGE][site * verts_per_face + site1] = buf_shell_start[GEONBR_TEDGE][site * verts_per_face];
      };         
   };
   
//----------------------------------------------------------------------------------------------------------------------------------------------------

// REDGE: 3/4 sites with 6/4 PPS
// 18/16 unique maps (1:1)
   buf_face_translation[GEONBR_REDGE] = new int**[exch_site_count[GEONBR_REDGE]];
   for (auto site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
      buf_face_translation[GEONBR_REDGE][site] = new int*[default_part_per_site[GEONBR_REDGE]];
      for (auto part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
         buf_face_translation[GEONBR_REDGE][site][part] = new int[buf_area[GEONBR_REDGE]];
      };
   };

// 1 unique starting shell (1:18/1:16)
   buf_shell_start[GEONBR_REDGE] = new int**[exch_site_count[GEONBR_REDGE]];
   buf_shell_start[GEONBR_REDGE][0] = new int*[default_part_per_site[GEONBR_REDGE]];
   for (auto part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
      buf_shell_start[GEONBR_REDGE][0][part] = buf_shell_tab + 1;
   };
   for (auto site = 1; site < exch_site_count[GEONBR_RFACE]; site++) {
      buf_shell_start[GEONBR_REDGE][site] = buf_shell_start[GEONBR_REDGE][0];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// VERTX: 6/8 sites with 12/8 PPS
// 18/16 unique maps (1:4), storage points to REDGE
   buf_face_translation[GEONBR_VERTX] = new int**[exch_site_count[GEONBR_VERTX]];
   for (auto site = 0; site < exch_site_count[GEONBR_VERTX] / 2; site++) {
      buf_face_translation[GEONBR_VERTX][site] = new int*[default_part_per_site[GEONBR_VERTX]];
      buf_face_translation[GEONBR_VERTX][site + exch_site_count[GEONBR_VERTX] / 2] = buf_face_translation[GEONBR_VERTX][site];
      for (auto part = 0; part < default_part_per_site[GEONBR_VERTX] / 2; part++) {
         buf_face_translation[GEONBR_VERTX][site][part] = buf_face_translation[GEONBR_REDGE][site][part];
         buf_face_translation[GEONBR_VERTX][site][part + default_part_per_site[GEONBR_VERTX] / 2] = buf_face_translation[GEONBR_VERTX][site][part];
      };
   };

// 4 unique starting shells (1:18/1:16)
   buf_shell_start[GEONBR_VERTX] = new int**[exch_site_count[GEONBR_VERTX]];
   for (auto site = 0; site < exch_site_count[GEONBR_VERTX] / verts_per_face; site++) {
      buf_shell_start[GEONBR_VERTX][site * verts_per_face] = new int*[default_part_per_site[GEONBR_VERTX]];
      for (auto part = 0; part < 2; part++) {
         for (auto part1 = 0; part1 < default_part_per_site[GEONBR_VERTX] / 2; part1++) {
            buf_shell_start[GEONBR_VERTX][site * verts_per_face][2 * part + part1] = buf_shell_tab + 2 * site + part;
         };
      };
      for (auto site1 = 1; site1 < exch_site_count[GEONBR_VERTX] / 2; site1++) {
         buf_shell_start[GEONBR_VERTX][site * verts_per_face + site1] = buf_shell_start[GEONBR_VERTX][site * verts_per_face];
      };         
   };
};

/*!
\author Vladimir Florinski
\date 12/23/2024
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::FreeStorage(void)
{
// Free up storage for the variables
   Delete2D(cons_vars);

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE
   if(buf_face_translation[GEONBR_TFACE] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_TFACE] / 2; site++) {
         for (auto part = 0; part < default_part_per_site[GEONBR_TFACE] / 2; part++) {
            delete[] buf_face_translation[GEONBR_TFACE][site][part];
         };
         delete[] buf_face_translation[GEONBR_TFACE][site];
      };
      delete[] buf_face_translation[GEONBR_TFACE];
   };

   if(buf_shell_start[GEONBR_TFACE] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_TFACE]; site++) {
         delete[] buf_shell_start[GEONBR_TFACE][site];
      };
      delete[] buf_shell_start[GEONBR_TFACE];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE
   if(buf_face_translation[GEONBR_RFACE] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
         for (auto part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
            delete[] buf_face_translation[GEONBR_RFACE][site][part];
         };
         delete[] buf_face_translation[GEONBR_RFACE][site];
      };
      delete[] buf_face_translation[GEONBR_RFACE];
   };

   if(buf_shell_start[GEONBR_RFACE] != nullptr) {
      delete[] buf_shell_start[GEONBR_RFACE][0];
      delete[] buf_shell_start[GEONBR_RFACE];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TEDGE
   if(buf_face_translation[GEONBR_TEDGE] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_TEDGE] / 2; site++) {
         delete[] buf_face_translation[GEONBR_TEDGE][site];
      };
      delete[] buf_face_translation[GEONBR_TEDGE];
   };

   if(buf_shell_start[GEONBR_TEDGE] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_TEDGE] / verts_per_face; site++) {
         delete[] buf_shell_start[GEONBR_TEDGE][site * verts_per_face];
      };
      delete[] buf_shell_start[GEONBR_TEDGE];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// REDGE
   if(buf_face_translation[GEONBR_REDGE] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
         for (auto part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
            delete[] buf_face_translation[GEONBR_REDGE][site][part];
         };
         delete[] buf_face_translation[GEONBR_REDGE][site];
      };
      delete[] buf_face_translation[GEONBR_REDGE];
   };

   if(buf_shell_start[GEONBR_REDGE] != nullptr) {
      delete[] buf_shell_start[GEONBR_REDGE][0];
      delete[] buf_shell_start[GEONBR_REDGE];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// VERTX
   if(buf_face_translation[GEONBR_VERTX] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_VERTX] / 2; site++) {
         delete[] buf_face_translation[GEONBR_VERTX][site];
      };
      delete[] buf_face_translation[GEONBR_VERTX];
   };

   if(buf_shell_start[GEONBR_VERTX] != nullptr) {
      for (auto site = 0; site < exch_site_count[GEONBR_VERTX] / verts_per_face; site++) {
         delete[] buf_shell_start[GEONBR_VERTX][site * verts_per_face];
      };
      delete[] buf_shell_start[GEONBR_VERTX];
   };
};

/*!
\author Vladimir Florinski
\date 06/26/2024
\param[in] ntype         Neighbor type
\param[in] exch_sites_in Array of pointers to exchange sites of this type
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::ImportExchangeSites(NeighborType ntype, const std::vector<std::shared_ptr<ExchangeSite<ConservedVariables>>>& exch_sites_in)
{
   exch_sites[ntype] = exch_sites_in;
/*
   exch_sites[ntype].clear();
   for (auto site : exch_sites_in) {
      exch_sites[ntype].emplace_back();
      exch_sites[ntype].back() = site;
   };
*/
};

/*!
\author Vladimir Florinski
\date 06/27/2024
\param[in] ntype Neighbor type
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::PackBuffers(NeighborType ntype) const
{
   int my_part, face, starting_shell;
   size_t bufidx;
   ConservedVariables* buffer;

   for (auto site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->part_lookup[block_index];
      buffer = exch_sites[ntype][site]->buffer_entry[my_part];
      starting_shell = *buf_shell_start[ntype][site][my_part];
      bufidx = 0;

      for (auto shell = starting_shell; shell < starting_shell + buf_height[ntype]; shell++) {
         for (auto face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
            face = buf_face_translation[ntype][site][my_part][face_buf];
            buffer[bufidx++] = cons_vars[shell][face];
         };
      };
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
   int my_part, face, starting_shell;
   size_t bufidx;
   ConservedVariables* buffer;

   for (auto site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->part_lookup[block_index];

// TODO check if this works at singular corners
      for (auto part = 0; part < default_part_per_site[ntype]; part++) {

// The part equal to ours will be skipped
         if (part == my_part) continue;

         buffer = exch_sites[ntype][site]->buffer_entry[part];
         starting_shell = *buf_shell_start[ntype][site][part];
         bufidx = 0;

// Copy each zone in the buffer to its proper place in the block
         for (auto shell = starting_shell; shell < starting_shell + buf_height[ntype]; shell++) {
            for (auto face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
               face = buf_face_translation[ntype][site][part][face_buf];
               cons_vars[shell][face] = buffer[bufidx++];
            };
         };
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
   for(auto face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
      std::cerr << std::setw(6) << face_buf << std::setw(6) << buf_face_translation[ntype][site][part][face_buf] << std::endl;
   };
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BufferedBlock protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2024
*/
template <int verts_per_face>
void BufferedBlock<verts_per_face>::ComputeBufferTranslations(void)
{
   int bufidx, i_origin, j_origin, vert_origin, i0, j0, i, j, face_origin, rot, my_part, foreign;
   std::pair<int, int> origin_arr[exch_site_count[GEONBR_RFACE] * default_part_per_site[GEONBR_RFACE]];

//----------------------------------------------------------------------------------------------------------------------------------------------------

// TFACE: 1 unique map
   rot = 0;
   bufidx = 0;

// Vertex origin for the sub-block
   i_origin = 2 * square_fill * ghost_width;
   j_origin = 2 * ghost_width;
   vert_origin = vert_index_sector[i_origin][j_origin];

// Face origin for the sub-block
   face_origin = vf_local[vert_origin][3];
   i0 = face_index_i[face_origin];
   j0 = face_index_j[face_origin];

   for (auto i_buf = 0; i_buf < buf_length[GEONBR_TFACE]; i_buf++) {
      for (auto j_buf = 0; j_buf <= MaxFaceJ(buf_length[GEONBR_TFACE], buf_width[GEONBR_TFACE], i_buf); j_buf++) {
         i = i0 + i_buf * rotated_faces[rot][0][0] + j_buf * rotated_faces[rot][0][1 + j_buf % square_fill];
         j = j0 + i_buf * rotated_faces[rot][1][0] + j_buf * rotated_faces[rot][1][1 + j_buf % square_fill];
         buf_face_translation[GEONBR_TFACE][0][0][bufidx++] = face_index_sector[i][j];
      };
   };   

//----------------------------------------------------------------------------------------------------------------------------------------------------

// RFACE/TEDGE: 6/8 unique maps
   origin_arr[0].first = total_length - (3 - square_fill) * ghost_width;
   origin_arr[0].second = 2 * ghost_width;

   origin_arr[0 + verts_per_face].first = 2 * ghost_width;
   origin_arr[0 + verts_per_face].second = 0;

   origin_arr[1].first = total_length - 2 * ghost_width;
   origin_arr[1].second = total_length - (square_fill + 1) * ghost_width;

   origin_arr[1 + verts_per_face].first = total_length;
   origin_arr[1 + verts_per_face].second = 2 * ghost_width;

   origin_arr[2].first = 2 * ghost_width;
   origin_arr[2].second = (verts_per_face == 3 ? ghost_width : total_length - 2 * ghost_width);

   origin_arr[2 + verts_per_face].first = total_length - 2 * ghost_width;
   origin_arr[2 + verts_per_face].second = total_length - (square_fill - 1) * ghost_width;

   if (verts_per_face == 4) {
      origin_arr[3].first = 2 * ghost_width;
      origin_arr[3].second = 2 * ghost_width;

      origin_arr[3 + verts_per_face].first = 0;
      origin_arr[3 + verts_per_face].second = total_length - 2 * ghost_width;
   };

   for (auto site = 0; site < exch_site_count[GEONBR_RFACE]; site++) {
      my_part = exch_sites[GEONBR_RFACE][site]->part_lookup[block_index];
      for (auto part = 0; part < default_part_per_site[GEONBR_RFACE]; part++) {
         foreign = (part == my_part ? 0 : 1);

// The base vertices are on the opposite sides of two sub-blocks, and "foreign" selects which side it is.
         i_origin = origin_arr[site + foreign * verts_per_face].first;
         j_origin = origin_arr[site + foreign * verts_per_face].second;
         vert_origin = vert_index_sector[i_origin][j_origin];

         face_origin = vf_local[vert_origin][(site * square_fill + foreign * cardinal_directions) % edges_per_vert];
         i0 = face_index_i[face_origin];
         j0 = face_index_j[face_origin];

         rot = (corner_rotation[site] + foreign * cardinal_directions) % edges_per_vert;

         bufidx = 0;
         for (auto i_buf = 0; i_buf < buf_length[GEONBR_RFACE]; i_buf++) {
            for(auto j_buf = 0; j_buf <= MaxFaceJ(buf_length[GEONBR_RFACE], buf_width[GEONBR_RFACE], i_buf); j_buf++) {
               i = i0 + i_buf * rotated_faces[rot][0][0] + j_buf * rotated_faces[rot][0][1 + j_buf % square_fill];
               j = j0 + i_buf * rotated_faces[rot][1][0] + j_buf * rotated_faces[rot][1][1 + j_buf % square_fill];
               buf_face_translation[GEONBR_RFACE][site][part][bufidx++] = face_index_sector[i][j];
            };
         };   
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// REDGE/VERTX: 18/16 unique maps
   origin_arr[0].first = square_fill * ghost_width;
   origin_arr[0].second = ghost_width;

   origin_arr[1].first = total_length - ghost_width;
   origin_arr[1].second = ghost_width;

   origin_arr[2].first = total_length - ghost_width;
   origin_arr[2].second = total_length - ghost_width;

   if (verts_per_face == 4) {
      origin_arr[3].first = square_fill * ghost_width;
      origin_arr[3].second = total_length - ghost_width;
   };

   for (auto site = 0; site < exch_site_count[GEONBR_REDGE]; site++) {
      i_origin = origin_arr[site].first;
      j_origin = origin_arr[site].second;
      vert_origin = vert_index_sector[i_origin][j_origin];

      my_part = exch_sites[GEONBR_REDGE][site]->part_lookup[block_index];
      for (auto part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
         face_origin = vf_local[vert_origin][(part + 3 * cardinal_directions - my_part) % edges_per_vert];
         i0 = face_index_i[face_origin];
         j0 = face_index_j[face_origin];

         rot = (part + edges_per_vert - my_part) % edges_per_vert;

         bufidx = 0;
         for (auto i_buf = 0; i_buf < buf_length[GEONBR_REDGE]; i_buf++) {
            for(auto j_buf = 0; j_buf <= MaxFaceJ(buf_length[GEONBR_REDGE], buf_width[GEONBR_REDGE], i_buf); j_buf++) {
               i = i0 + i_buf * rotated_faces[rot][0][0] + j_buf * rotated_faces[rot][0][1 + j_buf % square_fill];
               j = j0 + i_buf * rotated_faces[rot][1][0] + j_buf * rotated_faces[rot][1][1 + j_buf % square_fill];
               buf_face_translation[GEONBR_REDGE][site][part][bufidx++] = face_index_sector[i][j];
            };
         };   
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Starting shells
   buf_shell_tab[0] = 0;
   buf_shell_tab[1] = ghost_height;
   buf_shell_tab[2] = n_shells_withghost - 2 * ghost_height;
   buf_shell_tab[3] = n_shells_withghost - ghost_height;
};

template class BufferedBlock<3>;
template class BufferedBlock<4>;

};
