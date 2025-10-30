/*!
\file buffered_block.cc
\brief Implements the buffered block class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <cstring>
#include <utility>

#include "geodesic/buffered_block.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BufferedBlock public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/08/2025
*/
template <int verts_per_face, typename _datatype>
BufferedBlock<verts_per_face, _datatype>::BufferedBlock(void)
                                       : StenciledBlock<verts_per_face>()
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Default constructing a BufferedBlock\n";
#endif
#endif
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] other Object to initialize from
\note The copy constructor for "StenciledBlock" calls its "SetDimensions()" and "AssociateMesh()" methods
*/
template <int verts_per_face, typename _datatype>
BufferedBlock<verts_per_face, _datatype>::BufferedBlock(const BufferedBlock& other)
                                       : StenciledBlock<verts_per_face>(static_cast<const StenciledBlock<verts_per_face>&>(other))
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Copy constructing a BufferedBlock\n";
#endif
#endif

   if (other.side_length == -1) return;

   SetDimensions(other.side_length, other.ghost_width, other.n_shells, other.ghost_height, true);
   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      exch_sites[ntype] = other.exch_sites[ntype];
   };
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] other Object to move into this
*/
template <int verts_per_face, typename _datatype>
BufferedBlock<verts_per_face, _datatype>::BufferedBlock(BufferedBlock&& other) noexcept
                                       : StenciledBlock<verts_per_face>(std::move(static_cast<StenciledBlock<verts_per_face>&&>(other)))
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Move constructing a BufferedBlock (moving the content)\n";
#endif
#endif

   if (other.side_length == -1) {
      PrintMessage(__FILE__, __LINE__, "Move constructor called, but the dimension of the moved object was not set", true);
      return;
   };

// Move the variables
   zone_cons = other.zone_cons;
   other.zone_cons = nullptr;

// Move the exchange sites
   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      exch_sites[ntype] = std::move(other.exch_sites[ntype]);
   };
   
// Static arrays must be copied
   std::memcpy(buf_length, other.buf_length, N_NBRTYPES * sizeof(int));
   std::memcpy(buf_width, other.buf_width, N_NBRTYPES * sizeof(int));
   std::memcpy(buf_area, other.buf_area, N_NBRTYPES * sizeof(int));
   std::memcpy(buf_height, other.buf_height, N_NBRTYPES * sizeof(int));

// Move the buffer maps
   std::memcpy(buf_face_translation, other.buf_face_translation, N_NBRTYPES * sizeof(int***));
   std::memcpy(buf_shell_start, other.buf_shell_start, N_NBRTYPES * sizeof(int**));
   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      other.buf_face_translation[ntype] = nullptr;
      other.buf_shell_start[ntype] = nullptr;
   };
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] width  Length of the side, without ghost cells
\param[in] wgohst Width of the ghost cell layer outside the sector
\param[in] height Hight of the block, without ghost cells
\param[in] hghost Number of ghost shells outside the slab
*/
template <int verts_per_face, typename _datatype>
BufferedBlock<verts_per_face, _datatype>::BufferedBlock(int width, int wghost, int height, int hghost)
                                       : StenciledBlock<verts_per_face>(width, wghost, height, hghost)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Argument constructing a BufferedBlock\n";
#endif
#endif

   SetDimensions(width, wghost, height, hghost, true);
};

/*!
\author Vladimir Florinski
\date 01/08/2025
*/
template <int verts_per_face, typename _datatype>
BufferedBlock<verts_per_face, _datatype>::~BufferedBlock()
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Destructing a BufferedBlock\n";
#endif
#endif

   FreeStorage();
};

/*!
\author Vladimir Florinski
\date 01/09/2025
\param[in] width     Length of the side, without ghost cells
\param[in] wghost    Width of the ghost cell layer outside the sector
\param[in] height    Hight of the block, without ghost cells
\param[in] hghost    Number of ghost shells outside the slab
\param[in] construct Set to true when called from a constructor
\note It is possible to have only one participant at a t-face exchange site that is on the simulation boundary, which in principle requires no exchange. We still count it here, and it is up to ExchangeSite to figure out how to handle this efficiently.
*/
template <int verts_per_face, typename _datatype>
void BufferedBlock<verts_per_face, _datatype>::SetDimensions(int width, int wghost, int height, int hghost, bool construct)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Setting dimensions " << width << " by " << height << " for a BufferedBlock\n";
#endif
#endif

// Call base method
   if (!construct) StenciledBlock<verts_per_face>::SetDimensions(width, wghost, height, hghost, false);

// Free up storage (not that of the base class) because this could be a repeat call
   FreeStorage();

// Space for the conserved variables
   zone_cons = Create2D<_datatype>(n_shells_withghost, n_faces_withghost);

// Buffer sizes: TFACE
   buf_length[GEONBR_TFACE] = side_length - (square_fill + 1) * ghost_width;
   buf_width[GEONBR_TFACE] = buf_length[GEONBR_TFACE];
   buf_area[GEONBR_TFACE] = Sqr(buf_length[GEONBR_TFACE]);
   buf_height[GEONBR_TFACE] = ghost_height;

// Buffer sizes: RFACE
   buf_length[GEONBR_RFACE] = side_length - (3 - square_fill) * ghost_width;
   buf_width[GEONBR_RFACE] = ghost_width;
   buf_area[GEONBR_RFACE] = buf_length[GEONBR_RFACE] * square_fill * buf_width[GEONBR_RFACE] - (square_fill - 1) * Sqr(ghost_width);
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

// Check the buffer sized and print an information message if some buffers are zero length (for very small blocks)
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3

   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      if (buf_area[ntype] * buf_height[ntype] <= 0) {
         PrintMessage(__FILE__, __LINE__, "Block " + std::to_string(block_index) + " buffer of type " + neighbor_text[ntype] + " has zero or negative size", true);
      };
   };

#endif
#endif

// Allocate storage for the _three_ unique buffer face maps
   for (auto ntype : {GEONBR_TFACE, GEONBR_RFACE, GEONBR_REDGE}) {

// Don't allocate memory if the buffer area is zero (for very small blocks)
      if (buf_area[ntype] <= 0) continue;
   
      buf_face_translation[ntype] = new int**[exch_site_count[ntype] / n_slab_parts[ntype]];
      for (auto site = 0; site < exch_site_count[ntype] / n_slab_parts[ntype]; site++) {
         buf_face_translation[ntype][site] = new int*[default_part_per_site[ntype] / n_slab_parts[ntype]];
         for (auto part = 0; part < default_part_per_site[ntype] / n_slab_parts[ntype]; part++) {
            buf_face_translation[ntype][site][part] = new int[buf_area[ntype]];
         };
      };
   };

// Use pointers for the remaining two face maps 
   buf_face_translation[GEONBR_TEDGE] = buf_face_translation[GEONBR_RFACE];
   buf_face_translation[GEONBR_VERTX] = buf_face_translation[GEONBR_REDGE];

// Starting unique shells for multi-slab sites
   buf_shell_start[GEONBR_TFACE] = new int*[2];
   for (auto site = 0; site < 2; site++) buf_shell_start[GEONBR_TFACE][site] = new int[2];
   buf_shell_start[GEONBR_TFACE][0][0] = 0;
   buf_shell_start[GEONBR_TFACE][0][1] = ghost_height;
   buf_shell_start[GEONBR_TFACE][1][0] = n_shells_withghost - 2 * ghost_height;
   buf_shell_start[GEONBR_TFACE][1][1] = n_shells_withghost - ghost_height;
   
// Use pointers for the remaining two multi-slab sites
   buf_shell_start[GEONBR_TEDGE] = buf_shell_start[GEONBR_TFACE];
   buf_shell_start[GEONBR_VERTX] = buf_shell_start[GEONBR_TFACE];

// Starting shell for single slab sites
   buf_shell_start[GEONBR_RFACE] = new int*[1];
   buf_shell_start[GEONBR_RFACE][0] = new int[1];
   buf_shell_start[GEONBR_RFACE][0][0] = 2 * ghost_height;

// Use a pointer for the remaining single slab site
   buf_shell_start[GEONBR_REDGE] = buf_shell_start[GEONBR_RFACE];
};

/*!
\author Vladimir Florinski
\date 12/23/2024
*/
template <int verts_per_face, typename _datatype>
void BufferedBlock<verts_per_face, _datatype>::FreeStorage(void)
{
// Free up storage for the variables
   Delete2D(zone_cons);

// Free up storage for buffer face maps
   for (auto ntype : {GEONBR_TFACE, GEONBR_RFACE, GEONBR_REDGE}) {
      if(buf_face_translation[ntype] != nullptr) {
         for (auto site = 0; site < exch_site_count[ntype] / n_slab_parts[ntype]; site++) {
            for (auto part = 0; part < default_part_per_site[ntype] / n_slab_parts[ntype]; part++) {
               delete[] buf_face_translation[ntype][site][part];
            }
            delete[] buf_face_translation[ntype][site];
         };
         delete[] buf_face_translation[ntype];
      };
   };

// Free up storage for starting shells
   if(buf_shell_start[GEONBR_TFACE] != nullptr) {
      for (auto site = 0; site < 2; site++) delete[] buf_shell_start[GEONBR_TFACE][site];
      delete[] buf_shell_start[GEONBR_TFACE];
   };
   if(buf_shell_start[GEONBR_RFACE] != nullptr) {
      delete[] buf_shell_start[GEONBR_RFACE][0];
      delete[] buf_shell_start[GEONBR_RFACE];
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Starting shells have fixed values
   buf_shell_tab[0] = 0;
   buf_shell_tab[1] = ghost_height;
   buf_shell_tab[2] = n_shells_withghost - 2 * ghost_height;
   buf_shell_tab[3] = n_shells_withghost - ghost_height;
   buf_shell_tab[4] = 2 * ghost_height;
};

/*!
\author Vladimir Florinski
\date 01/09/2025
\param[in] ntype         Neighbor type
\param[in] exch_sites_in Array of pointers to exchange sites of this type
*/
template <int verts_per_face, typename _datatype>
void BufferedBlock<verts_per_face, _datatype>::ImportExchangeSites(NeighborType ntype, const std::vector<std::shared_ptr<ExchangeSite<_datatype>>>& exch_sites_in)
{
   exch_sites[ntype] = exch_sites_in;

// Compute the buffer translations
   int bufidx, i_origin, j_origin, vert_origin, i0, j0, i, j, face_origin, rot, my_part, missing_part, old_part, foreign;
   int** buf_face_translation_tmp;
   std::pair<int, int> origin_arr[exch_site_count[GEONBR_RFACE] * default_part_per_site[GEONBR_RFACE]];
//   int site_stride = exch_site_count[ntype] / n_slab_parts[ntype];
   int part_stride = default_part_per_site[ntype] / n_slab_parts[ntype];

   switch(ntype) {
   
// TFACE: 1 unique map
   case GEONBR_TFACE:
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
         i = i0 + i_buf * rotated_faces[rot][0][0];
         j = j0 + i_buf * rotated_faces[rot][1][0];
         for (auto j_buf = 0; j_buf <= MaxFaceJ(buf_length[GEONBR_TFACE], buf_width[GEONBR_TFACE], i_buf); j_buf++) {
            buf_face_translation[GEONBR_TFACE][0][0][bufidx++] = face_index_sector[i][j];
            i += rotated_faces[rot][0][1 + j_buf % square_fill];
            j += rotated_faces[rot][1][1 + j_buf % square_fill];
         };
      };

// Correction for the bottom block where the part numbers are shifted (2/3 are now 0/1, etc.)
      if (border_type[0]) buf_shell_start[GEONBR_TFACE][0][0] = ghost_height;

      break;

// RFACE/TEDGE: 6/8 unique maps. GEONBR_TEDGE sites are used because GEONBR_RFACE sites can be left un-initialized (for small size blocks).
   case GEONBR_TEDGE:
      origin_arr[0].first = total_length - (3 - square_fill) * ghost_width;
      origin_arr[0].second = 2 * ghost_width;

      origin_arr[0 + verts_per_face].first = 2 * ghost_width;
      origin_arr[0 + verts_per_face].second = 0;

      origin_arr[1].first = total_length - 2 * ghost_width;
      origin_arr[1].second = total_length - (square_fill + 1) * ghost_width;

      origin_arr[1 + verts_per_face].first = total_length;
      origin_arr[1 + verts_per_face].second = 2 * ghost_width;

      origin_arr[2].first = (1 + square_fill) * ghost_width;
      origin_arr[2].second = (verts_per_face == 3 ? ghost_width : total_length - 2 * ghost_width);

      origin_arr[2 + verts_per_face].first = total_length - 2 * ghost_width;
      origin_arr[2 + verts_per_face].second = total_length - 2 * (square_fill - 1) * ghost_width;

      if (verts_per_face == 4) {
         origin_arr[3].first = 2 * ghost_width;
         origin_arr[3].second = 2 * ghost_width;

         origin_arr[3 + verts_per_face].first = 0;
         origin_arr[3 + verts_per_face].second = total_length - 2 * ghost_width;
      };

      for (auto site = 0; site < exch_site_count[ntype] / n_slab_parts[ntype]; site++) {
         my_part = exch_sites[ntype][site]->PartOfLabel(block_index) % part_stride;
         for (auto part = 0; part < default_part_per_site[ntype] / n_slab_parts[ntype]; part++) {
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
               i = i0 + i_buf * rotated_faces[rot][0][0];
               j = j0 + i_buf * rotated_faces[rot][1][0];
               for (auto j_buf = 0; j_buf <= MaxFaceJ(buf_length[GEONBR_RFACE], buf_width[GEONBR_RFACE], i_buf); j_buf++) {
                  buf_face_translation[GEONBR_RFACE][site][part][bufidx++] = face_index_sector[i][j];
                  i += rotated_faces[rot][0][1 + j_buf % square_fill];
                  j += rotated_faces[rot][1][1 + j_buf % square_fill];
               };
            };   
         };
      };

      break;

// REDGE/VERTX: 18(15)/16(12) unique maps. GEONBR_VERTX sites are used because GEONBR_REDGE sites can be left un-initialized (for small size blocks).
   case GEONBR_VERTX:
      origin_arr[0].first = square_fill * ghost_width;
      origin_arr[0].second = ghost_width;

      origin_arr[1].first = total_length - ghost_width;
      origin_arr[1].second = ghost_width;

      origin_arr[2].first = total_length - ghost_width;
      origin_arr[2].second = MaxVertJ(total_length, origin_arr[2].first) - ghost_width;

      if (verts_per_face == 4) {
         origin_arr[3].first = square_fill * ghost_width;
         origin_arr[3].second = total_length - ghost_width;
      };

      for (auto site = 0; site < exch_site_count[ntype] / n_slab_parts[ntype]; site++) {
         i_origin = origin_arr[site].first;
         j_origin = origin_arr[site].second;
         vert_origin = vert_index_sector[i_origin][j_origin];

         my_part = exch_sites[ntype][site]->PartOfLabel(block_index) % part_stride;
         missing_part = -1;

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
         std::cerr << "Site " << site << " index " << exch_sites[GEONBR_REDGE][site]->GetIndex()
                   << " my part is " << exch_sites[GEONBR_REDGE][site]->PartOfLabel(block_index) << std::endl;
#endif
#endif

         for (auto part = 0; part < default_part_per_site[ntype] / n_slab_parts[ntype]; part++) {
            rot = (corner_rotation[site] + part + cardinal_directions - my_part + edges_per_vert) % edges_per_vert;
            face_origin = vf_local[vert_origin][(rot + cardinal_directions) % edges_per_vert];

// Missing part at a singular corner is flagged
            if (face_origin == -1) {
               missing_part = part;
               continue;
            };

            i0 = face_index_i[face_origin];
            j0 = face_index_j[face_origin];

            bufidx = 0;
            for (auto i_buf = 0; i_buf < buf_length[GEONBR_REDGE]; i_buf++) {
               i = i0 + i_buf * rotated_faces[rot][0][0];
               j = j0 + i_buf * rotated_faces[rot][1][0];
               for (auto j_buf = 0; j_buf <= MaxFaceJ(buf_length[GEONBR_REDGE], buf_width[GEONBR_REDGE], i_buf); j_buf++) {
                  buf_face_translation[GEONBR_REDGE][site][part][bufidx++] = face_index_sector[i][j];
                  i += rotated_faces[rot][0][1 + j_buf % square_fill];
                  j += rotated_faces[rot][1][1 + j_buf % square_fill];
               };
            };   
         };

// Compress the translation arrays to remove the hole at the "missing_part". We use the fact that "my_part" is never a missing part.
         if (missing_part != -1) {
            buf_face_translation_tmp = new int*[default_part_per_site[GEONBR_REDGE]];
            buf_face_translation_tmp[my_part] = buf_face_translation[GEONBR_REDGE][site][my_part];
            buf_face_translation_tmp[default_part_per_site[GEONBR_REDGE] - 1] = buf_face_translation[GEONBR_REDGE][site][missing_part];

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
            std::cerr << "Part " << my_part << " is mapped to " << my_part << std::endl;
            std::cerr << "Part " << default_part_per_site[GEONBR_REDGE] - 1 << " is mapped to " << missing_part << std::endl;
#endif
#endif

// Half loop toward decreasing "part"
            old_part = my_part + default_part_per_site[GEONBR_REDGE];
            for (auto part = my_part - 1; part >= 0; part--) {
               old_part--;
               if(part == missing_part) old_part--;
               buf_face_translation_tmp[part] = buf_face_translation[GEONBR_REDGE][site][old_part % default_part_per_site[GEONBR_REDGE]];

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
               std::cerr << "Part " << part << " is mapped to " << old_part % default_part_per_site[GEONBR_REDGE] << std::endl;
#endif
#endif
            };

// Half loop toward increasing "part"
            old_part = my_part + default_part_per_site[GEONBR_REDGE];
            for (auto part = my_part + 1; part < default_part_per_site[GEONBR_REDGE] - 1; part++) {
               old_part++;
               if(part == missing_part) old_part++;
               buf_face_translation_tmp[part] = buf_face_translation[GEONBR_REDGE][site][old_part % default_part_per_site[GEONBR_REDGE]];

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
               std::cerr << "Part " << part << " is mapped to " << old_part % default_part_per_site[GEONBR_REDGE] << std::endl;
#endif
#endif
            };
            
// Replace the map
            for (auto part = 0; part < default_part_per_site[GEONBR_REDGE]; part++) {
               buf_face_translation[GEONBR_REDGE][site][part] = buf_face_translation_tmp[part];
            };
            delete[] buf_face_translation_tmp;
         };
      };

      break;

   default:
      break;

   };
};

/*!
\author Vladimir Florinski
\date 01/10/2025
\param[in] ntype Neighbor type
*/
template <int verts_per_face, typename _datatype>
void BufferedBlock<verts_per_face, _datatype>::PackBuffers(NeighborType ntype, int test_block) const
{
// Skip the operation if the buffer size is zero (for small size blocks)
   if ((buf_area[ntype] <= 0) || (buf_height[ntype] <= 0)) return;

   int my_part, face, starting_shell;
   size_t bufidx;
   int site_stride = exch_site_count[ntype] / n_slab_parts[ntype];
   int part_stride = default_part_per_site[ntype] / n_slab_parts[ntype];

   _datatype* buffer;

   for (auto site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->PartOfLabel(block_index);

// Find the starting position in the correct buffer
      buffer = exch_sites[ntype][site]->BufferAddress(my_part);
      starting_shell = buf_shell_start[ntype][site / site_shel_mult[ntype]][my_part / n_sect_parts[ntype]];
      bufidx = 0;

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 2
      if(block_index == test_block) {
         std::cerr << "Packing block " << block_index << " site " << site << " type " << neighbor_text[ntype] << " part " << my_part << std::endl;
      };
#endif
#endif

// Copy each zone into the buffer
      for (auto shell = starting_shell; shell < starting_shell + buf_height[ntype]; shell++) {
         for (auto face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
            face = buf_face_translation[ntype][site % site_stride][my_part % part_stride][face_buf];
            buffer[bufidx++] = zone_cons[shell][face];
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 01/10/2025
\param[in] ntype Neighbor type
*/
template <int verts_per_face, typename _datatype>
void BufferedBlock<verts_per_face, _datatype>::UnPackBuffers(NeighborType ntype, int test_block)
{
// Skip the operation if the buffer size is zero (for small size blocks)
   if ((buf_area[ntype] <= 0) || (buf_height[ntype] <= 0)) return;

   int my_part, face, starting_shell;
   size_t bufidx;
   int site_stride = exch_site_count[ntype] / n_slab_parts[ntype];
   int part_stride = default_part_per_site[ntype] / n_slab_parts[ntype];

   _datatype* buffer;

   for (auto site = 0; site < exch_site_count[ntype]; site++) {
      my_part = exch_sites[ntype][site]->PartOfLabel(block_index);

      for (auto part = 0; part < exch_sites[ntype][site]->GetPartCount(); part++) {

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 2
         if(block_index == test_block) {
            std::cerr << "Unpacking block " << block_index << " site " << site << " type " << neighbor_text[ntype] << " part " << part << (part == my_part ? " *" : "") << std::endl;
         };
#endif
#endif

// The part equal to ours is skipped
         if (part == my_part) continue;

// Find the starting position in the correct buffer
         buffer = exch_sites[ntype][site]->BufferAddress(part);
         starting_shell = buf_shell_start[ntype][site / site_shel_mult[ntype]][part / n_sect_parts[ntype]]; 
         bufidx = 0;

// Copy each zone from the buffer
         for (auto shell = starting_shell; shell < starting_shell + buf_height[ntype]; shell++) {
            for (auto face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
               face = buf_face_translation[ntype][site % site_stride][part % part_stride][face_buf];
               zone_cons[shell][face] = buffer[bufidx++];
            };
         };
         actual_part++;
      };
   };
};

#ifdef USE_SILO

/*!
\author Vladimir Florinski
\date 03/19/2025
\param[in] silofile   SILO database file object
\param[in] phys_units Use physical units for output
\return Error code from SILO file I/O (zero if no error)
*/
template <int verts_per_face, typename _datatype>
template <typename datatype, std::enable_if_t<std::is_arithmetic<datatype>::value, bool>>
int BufferedBlock<verts_per_face, _datatype>::WriteSiloData(DBfile* silofile, bool phys_units) const
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 2
   std::cerr << "Writing data for block " << block_index << "\n";
#endif
#endif

   int err, k, k1, k2, face, face_silo, idx1, idx2, n_zones, n_intzones;

   n_zones = this->n_faces_silo * n_shells;
   if (!border_type[0]) n_zones += this->n_faces_silo;
   if (!border_type[1]) n_zones += this->n_faces_silo;
   n_intzones = n_faces * n_shells;

// A 1D array containing all data in the radial-index-first format.
   double* linear_data = new double[n_zones];

// The first index tracks the interior zones, while the second enumerates the ghost zones.
   idx1 = 0;
   idx2 = n_intzones;
   k1 = (border_type[0] ? ghost_height : ghost_height - 1);
   k2 = (border_type[1] ? n_shells_withghost - ghost_height - 1 : n_shells_withghost - ghost_height);

// Copy the variable into the 1D "linear_data" array for SILO consumption.
   for (k = k1; k <= k2; k++) {
      for (face_silo = 0; face_silo < this->n_faces_silo; face_silo++) {      
         face = this->silo_to_face[face_silo];

// Add to the appropriate place in the array
         if (IsInteriorShellOfSlab(k) && IsInteriorFaceOfSector(face)) linear_data[idx1++] = zone_cons[k][face];
         else linear_data[idx2++] = zone_cons[k][face];
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Unique indices for the variables
   char blck_num[blck_length + 1];
   sprintf(blck_num, blck_format.c_str(), block_index);
   std::string meshname = um_base + "_" + blck_num;
   std::string varname  = uv_base + "_" + blck_num;

// Write the unstructured var object to the SILO database
   err = DBPutUcdvar1(silofile, varname.c_str(), meshname.c_str(), linear_data, n_zones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

// Clean up
   delete[] linear_data;
   return err;
};

/*!
\author Vladimir Florinski
\date 03/19/2025
\param[in] silofile   SILO database file object
\param[in] phys_units Use physical units for output
\return Error code from file I/O (zero if no error)
*/
template <int verts_per_face, typename _datatype>
template <typename datatype, std::enable_if_t<std::is_arithmetic<datatype>::value, bool>>
int BufferedBlock<verts_per_face, _datatype>::WriteSilo(DBfile* silofile, bool phys_units) const
{
   int err1, err2;
   
   err1 = StenciledBlock<verts_per_face>::WriteSilo(silofile, phys_units);
   err2 = WriteSiloData(silofile, phys_units);
   return err1 | err2;
};

#endif

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 06/18/2024
\param[in] ntype Neighbor type
\param[in] site  Site index of this type
\param[in] part  Participant index in this site
*/
template <int verts_per_face, typename _datatype>
void BufferedBlock<verts_per_face, _datatype>::PrintBufferMaps(NeighborType ntype) const
{
   int site_stride = exch_site_count[ntype] / n_slab_parts[ntype];
   int part_stride = default_part_per_site[ntype] / n_slab_parts[ntype];

   std::cerr << std::endl;
   std::cerr << "Printing buffer face map of type " << neighbor_text[ntype] << " for block " << block_index << ":\n";
   for (auto site = 0; site < exch_site_count[ntype]; site++) {
      std::cerr << "   Site " << site << std::endl;
      for (auto part = 0; part < default_part_per_site[ntype]; part++) {
         std::cerr << "      Part " << part << "\n      ";
         std::cerr << std::setw(5) << buf_shell_start[ntype][site / site_shel_mult[ntype]][part / n_sect_parts[ntype]] << " | ";
         for (auto face_buf = 0; face_buf < buf_area[ntype]; face_buf++) {
            std::cerr << std::setw(5) << buf_face_translation[ntype][site % site_stride][part % part_stride][face_buf];
         };
         std::cerr << std::endl;
      };
   };
};

/*!
\author Vladimir Florinski
\date 03/19/2025
\param[in] val Value to fill with (-1 for face index)
*/
template <int verts_per_face, typename _datatype>
template <typename datatype, std::enable_if_t<std::is_arithmetic<datatype>::value, bool>>
void BufferedBlock<verts_per_face, _datatype>::FillUniform(_datatype val)
{
   for (auto shell = 0; shell < n_shells_withghost; shell++) {
      for (auto face = 0; face < n_faces_withghost; face++) {

// Assign only the internal zones (one could use shell information as well).
         if((shell >= ghost_height) && (shell < n_shells + ghost_height) && BITS_RAISED(face_mask[face], GEOELM_INTR)) {
            zone_cons[shell][face] = (val == -1 ? face : val);
         }
         else zone_cons[shell][face] = -1;
      };
   };
};

/*!
\author Vladimir Florinski
\date 01/10/2025
*/
template <>
void BufferedBlock<3, int>::PrintContents(void) const
{
   int face;
   std::cerr << std::endl;
   std::cerr << "Printing the contents of block " << block_index << std::endl;
   for (auto shell = n_shells_withghost - 1; shell >= 0; shell--) {
      for (auto twoj = total_length - 1; twoj >= 0; twoj--) {
         for (auto k = 0; k < total_length - twoj; k++) std::cerr << "    ";

         for (auto i = 0; i < total_length; i++) {
            if(i) {
               face = face_index_sector[i][2 * twoj + 1];
               if(face == -1) std::cerr << "    ";
               else std::cerr << std::setw(4) << zone_cons[shell][face];
            };
            face = face_index_sector[i][2 * twoj];
            if(face == -1) std::cerr << "    ";
            else std::cerr << std::setw(4) << zone_cons[shell][face];
         };
         std::cerr << std::endl;
      };
      std::cerr << std::endl;
   };
};

/*!
\author Vladimir Florinski
\date 01/10/2025
*/
template <>
void BufferedBlock<4, int>::PrintContents(void) const
{
   std::cerr << std::endl;
   std::cerr << "Printing the contents of block " << block_index << std::endl;
   for (auto shell = n_shells_withghost - 1; shell >= 0; shell--) {
      for (auto i = 0; i < total_length; i++) {
         std::cerr << std::endl;
         for (auto j = 0; j <= MaxFaceJ(total_length, i); j++) {
            std::cerr << std::setw(4) << zone_cons[shell][face_index_sector[i][j]];
         };
      };
      std::cerr << std::endl;
   };
};

#endif

template class BufferedBlock<3, int>;

#ifdef USE_SILO
template int BufferedBlock<3, int>::WriteSilo(DBfile* silofile, bool phys_units) const;
#endif

#ifdef GEO_DEBUG
template void BufferedBlock<3, int>::FillUniform(int val);
#endif


//template class BufferedBlock<4, int>;

};
