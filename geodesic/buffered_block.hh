/*!
\file buffered_block.hh
\brief Extends the "StenciledBlock" class with conserved variables and connections to neighbors
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BUFFERED_BLOCK_HH
#define SPECTRUM_BUFFERED_BLOCK_HH

#include "common/exchange_site.hh"
#include "geodesic/stenciled_block.hh"

namespace Spectrum {

//! Number of neighbor types
#define N_NBRTYPES 5

//! Neighbor types
enum NeighborType {
   GEONBR_TFACE,
   GEONBR_RFACE,
   GEONBR_TEDGE,
   GEONBR_REDGE,
   GEONBR_VERTX
};

/*!
\brief Compute the number of exchange sites
\author Vladimir Florinski
\date 06/24/2024
\param[in]  p          Vertices per face
\param[out] site_count Array of site counts of each type
*/
inline void ExchangeSiteCount(int p, int* site_count)
{
   site_count[GEONBR_TFACE] = 2;
   site_count[GEONBR_RFACE] = p;
   site_count[GEONBR_TEDGE] = 2 * p;
   site_count[GEONBR_REDGE] = p;
   site_count[GEONBR_VERTX] = 2 * p;
};

/*!
\brief Compute the number of participants of each exchange site
\author Vladimir Florinski
\date 06/24/2024
\param[in]  q          Edges meeting at each vertex
\param[out] part_count Array of participant counts of each type
*/
inline void ExchangePartCount(int q, int* part_count)
{
   part_count[GEONBR_TFACE] = 2;
   part_count[GEONBR_RFACE] = 2;
   part_count[GEONBR_TEDGE] = 4;
   part_count[GEONBR_REDGE] = q;
   part_count[GEONBR_VERTX] = 2 * q;
};

/*!
\brief A class describing buffered and stenciled grid blocks
\author Vladimir Florinski

Extends the StencilBlock class with (a) storage for conserved variables and (b) buffers for variable exchange between neighboring blocks
*/
template <int verts_per_face, typename _datatype>
class BufferedBlock : public StenciledBlock<verts_per_face>
{
public:

   using datatype = _datatype;

protected:

   using SphericalSlab::n_shells;
   using SphericalSlab::ghost_height;
   using SphericalSlab::n_shells_withghost;
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
   using StenciledBlock<verts_per_face>::MaxFaceJ;

//! Number of exchange sites
   int exch_site_count[N_NBRTYPES];

//! Default (not actual) number of parts per site
   int default_part_per_site[N_NBRTYPES];

//! Pointers of exchange sites
   std::vector<std::shared_ptr<ExchangeSite<datatype>>> exch_sites[N_NBRTYPES];

//! Buffer side length
   int buf_length[N_NBRTYPES];

//! Buffer width
   int buf_width[N_NBRTYPES];

//! Buffer area
   int buf_area[N_NBRTYPES];

//! Buffer height
   int buf_height[N_NBRTYPES];

//! Buffer volume
   int buf_volume[N_NBRTYPES];

//! Conserved variables storage
   datatype** zone_cons = nullptr;

//! Translation tables for each buffer (site, part, face)
   int*** buf_face_translation[N_NBRTYPES] = {nullptr};

//! Starting shell table
   int buf_shell_tab[5];

//! Starting shell indices for each buffer (site, part), accessed through a pointer
   int*** buf_shell_start[N_NBRTYPES] = {nullptr};

public:

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
   void ImportExchangeSites(NeighborType ntype, const std::vector<std::shared_ptr<ExchangeSite<datatype>>>& exch_sites_in);

//! Return the number of exchange sites
   int GetExchangeSiteCount(NeighborType ntype) const {return exch_site_count[ntype];};

//! Return the number of exchange sites
   int GetPartCount(NeighborType ntype) const {return default_part_per_site[ntype];};

//! Return the size of the buffer (in units of datatype)
   int GetBufferSize(NeighborType ntype) const {return buf_volume[ntype];};

//! Pack ghost regions for exchange
   void PackBuffers(NeighborType ntype) const;

//! Assign ghost values after exchange
   void UnPackBuffers(NeighborType ntype);

#ifdef GEO_DEBUG

//! Print a buffer map
   void PrintBufferMap(NeighborType ntype, int site, int part) const;

//! Fill the data with index values
   void FillWithIndexData(void);

//! Print the contents of the entire block
   void PrintContents(void);

#endif

};

};

#endif
