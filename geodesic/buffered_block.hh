/*!
\file buffered_block.hh
\brief Extends the "StenciledBlock" class with conserved variables and connections to neighbors
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BUFFERED_BLOCK_HH
#define SPECTRUM_BUFFERED_BLOCK_HH

#include <common/exchange_site.hh>
#include <geodesic/stenciled_block.hh>
#include <geodesic/neighbors.hh>

namespace Spectrum {

#ifdef USE_SILO
//! Pattern for _test_ variable naming
const std::string uv_base = "geouv";
#endif

/*!
\brief A class describing buffered and stenciled grid blocks
\author Vladimir Florinski

Extends the StencilBlock class with (a) storage for conserved variables and (b) buffers for variable exchange between neighboring blocks
*/
template <int verts_per_face, typename _datatype>
class BufferedBlock : public StenciledBlock<verts_per_face>
{
protected:

   using SphericalSlab::n_shells;
   using SphericalSlab::ghost_height;
   using SphericalSlab::n_shells_withghost;
   using SphericalSlab::IsInteriorShellOfSlab;
   using PolygonalAddressing<verts_per_face>::edges_per_vert;
   using PolygonalAddressing<verts_per_face>::square_fill;
   using PolygonalAddressing<verts_per_face>::cardinal_directions;
   using PolygonalAddressing<verts_per_face>::FaceCount;
   using GeodesicSector<verts_per_face>::n_faces_withghost;
   using GeodesicSector<verts_per_face>::side_length;
   using GeodesicSector<verts_per_face>::total_length;
   using GeodesicSector<verts_per_face>::ghost_width;
   using GeodesicSector<verts_per_face>::n_faces;
   using GeodesicSector<verts_per_face>::vert_index_sector;
   using GeodesicSector<verts_per_face>::face_index_sector;
   using GeodesicSector<verts_per_face>::face_index_i;
   using GeodesicSector<verts_per_face>::face_index_j;
   using GeodesicSector<verts_per_face>::vf_local;
   using GeodesicSector<verts_per_face>::MaxVertJ;
   using GridBlock<verts_per_face>::block_index;
   using GridBlock<verts_per_face>::corner_rotation;
   using GridBlock<verts_per_face>::rotated_faces;
   using GridBlock<verts_per_face>::face_mask;
   using GridBlock<verts_per_face>::border_type;
   using GridBlock<verts_per_face>::IsInteriorFaceOfSector;
   using StenciledBlock<verts_per_face>::MaxFaceJ;

//! Conserved variables storage
   _datatype** zone_cons = nullptr;

//! Number of exchange sites
   static constexpr int exch_site_count[N_NBRTYPES] = {2, verts_per_face, 2 * verts_per_face, verts_per_face, 2 * verts_per_face};

//! Default (not actual) number of parts per site
   static constexpr int default_part_per_site[N_NBRTYPES] = {2, 2, 4, edges_per_vert, 2 * edges_per_vert};

//! Number of exchange sites with identical starting shells
   static constexpr int site_shel_mult[N_NBRTYPES] = {1, verts_per_face, verts_per_face, verts_per_face, verts_per_face};

//! Number of participating sectors per site
   static constexpr int n_sect_parts[N_NBRTYPES] = {1, 2, 2, edges_per_vert, edges_per_vert};

//! Pointers of exchange sites
   std::vector<std::shared_ptr<ExchangeSite<_datatype>>> exch_sites[N_NBRTYPES];

//! Buffer side length
   int buf_length[N_NBRTYPES];

//! Buffer width (measured in vertex units - multiply by "square_fill" to get the area)
   int buf_width[N_NBRTYPES];

//! Buffer area
   int buf_area[N_NBRTYPES];

//! Buffer height
   int buf_height[N_NBRTYPES];

//! Translation tables for each buffer (site, part, face)
   int*** buf_face_translation[N_NBRTYPES] = {nullptr};

//! Starting shell indices for each buffer (site, part), accessed through a pointer
   int** buf_shell_start[N_NBRTYPES] = {nullptr};

//! Write the data description to a SILO database
#ifdef USE_SILO
   template <typename datatype = _datatype, std::enable_if_t<std::is_arithmetic<datatype>::value, bool> = true>
   int WriteSiloData(DBfile* silofile, bool phys_units) const;
#endif

public:

   using datatype = _datatype;

//! Default constructor
   BufferedBlock(void);

//! Copy constructor
   BufferedBlock(const BufferedBlock& other);

//! Move constructor
   BufferedBlock(BufferedBlock&& other) noexcept;

//! Constructor with arguments
   BufferedBlock(int width, int wghost, int height, int hghost);

//! Destructor
   ~BufferedBlock();

//! Allocate memory
   void SetDimensions(int width, int wghost, int height, int hghost, bool construct);

//! Free all dynamically allocated memory
   void FreeStorage(void);

//! Assign exchange site objects and generate address maps for the buffers
   void ImportExchangeSites(NeighborType ntype, const std::vector<std::shared_ptr<ExchangeSite<_datatype>>>& exch_sites_in);

//! Return the number of exchange sites
   static constexpr int GetExchangeSiteCount(NeighborType ntype) {return exch_site_count[ntype];};

//! Return the number of exchange sites
   static constexpr int GetPartCount(NeighborType ntype) {return default_part_per_site[ntype];};

//! Return the size of the buffer (in units of "_datatype")
   int GetBufferSize(NeighborType ntype) const {return buf_area[ntype] * buf_height[ntype];};

//! Pack ghost regions for exchange
   void PackBuffers(NeighborType ntype, int test_block = -1) const;

//! Assign ghost values after exchange
   void UnPackBuffers(NeighborType ntype, int test_block = -1);

#ifdef USE_SILO
//! Write the entire block to a SILO database
   template <typename datatype = _datatype, std::enable_if_t<std::is_arithmetic<datatype>::value, bool> = true>
   int WriteSilo(DBfile* silofile, bool phys_units) const;
#endif

#ifdef GEO_DEBUG

//! Print a buffer map
   void PrintBufferMaps(NeighborType ntype) const;

//! Fill the data with a uniform value or index data
   template <typename datatype = _datatype, std::enable_if_t<std::is_arithmetic<datatype>::value, bool> = true>
   void FillUniform(_datatype val = -1);

//! Print the contents of the entire block
   void PrintContents(void) const;



// TODO
   bool use_storage;
   void SetStorageType(bool);




#endif

};

};

#endif