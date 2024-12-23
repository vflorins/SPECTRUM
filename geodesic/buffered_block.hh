/*!
\file buffered_block.hh
\brief Extends the "StenciledBlock" class with conserved variables and connections to neighbors
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BUFFERED_BLOCK_HH
#define SPECTRUM_BUFFERED_BLOCK_HH

#include "geodesic/stenciled_block.hh"
#include "fluid/variables.hh"
#include "common/exchange_site.hh"
#include "common/mpi_config.hh"

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
template <int verts_per_face>
class BufferedBlock : public StenciledBlock<verts_per_face>
{
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
   using GridBlock<verts_per_face>::block_index;
   using GridBlock<verts_per_face>::corner_rotation;
   using GridBlock<verts_per_face>::rotated_faces;

//! Number of exchange sites
   int exch_site_count[N_NBRTYPES];

//! Default (not actual) number of parts per site
   int default_part_per_site[N_NBRTYPES];

//! Pointers of exchange sites
   std::vector<std::shared_ptr<ExchangeSite<ConservedVariables>>> exch_sites[N_NBRTYPES];

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
   ConservedVariables** cons_vars = nullptr;

//! Translation tables for each buffer (site, part, face)
   int*** buf_face_translation[N_NBRTYPES] = {nullptr};

//! Starting shell table
   int buf_shell_tab[4];

//! Starting shell indices for each buffer (site, part), accessed through a pointer
   int*** buf_shell_start[N_NBRTYPES] = {nullptr};

//! Maximum value of the vertex j-index for the given vertex i-index for a truncated block
   SPECTRUM_DEVICE_FUNC int MaxVertJ(int len, int height, int i) const;

//! Maximum value of the face j-index for the given face i-index for a truncated block
   SPECTRUM_DEVICE_FUNC int MaxFaceJ(int len, int height, int i) const;

//! Generate address maps for the buffers
   void ComputeBufferTranslations(void);

public:

//! Default constructor
   BufferedBlock(void) = default;

//! Copy constructor
   BufferedBlock(const BufferedBlock& other);

//! Constructor with arguments
   BufferedBlock(int width, int wghost, int height, int hghost);

//! Destructor
   ~BufferedBlock();

//! Allocate memory
   void SetDimensions(int width, int wghost, int height, int hghost, bool construct);

//! Free all dynamically allocated memory
   void FreeStorage(void);

//! Assign exchange site objects
   void ImportExchangeSites(NeighborType ntype, const std::vector<std::shared_ptr<ExchangeSite<ConservedVariables>>>& exch_sites_in);

//! Return the size of the buffer (in units of ConservedVariables)
   int GetBufferSize(NeighborType ntype) const {return buf_volume[ntype];};

//! Pack ghost regions for exchange
   void PackBuffers(NeighborType ntype) const;

//! Assign ghost values after exchange
   void UnPackBuffers(NeighborType ntype);

#ifdef GEO_DEBUG
//! Print a buffer map
   void PrintBufferMap(NeighborType ntype, int site, int part) const;
#endif

};

/*!
\author Vladimir Florinski
\date 06/14/2024
\param[in] len    Side length
\param[in] height Height of the block in _vertices_
\param[in] i      i-index
\return Largest j-index for this i
*/
template <>
SPECTRUM_DEVICE_FUNC inline int BufferedBlock<3>::MaxVertJ(int len, int height, int i) const
{
   return std::min(i, height);
};

/*!
\author Vladimir Florinski
\date 06/14/2024
\param[in] len    Side length
\param[in] height Height of the block in _vertices_
\param[in] i      i-index
\return Largest j-index for this i
*/
template <>
SPECTRUM_DEVICE_FUNC inline int BufferedBlock<4>::MaxVertJ(int len, int height, int i) const
{
   return height;
};

/*!
\author Vladimir Florinski
\date 06/14/2024
\param[in] len    Side length
\param[in] height Height of the block in _vertices_
\param[in] i      i-index
\return Largest j-index for this i
*/
template <>
SPECTRUM_DEVICE_FUNC inline int BufferedBlock<3>::MaxFaceJ(int len, int height, int i) const
{
   return (i < height ? square_fill * i : square_fill * height - 1);
};

/*!
\author Vladimir Florinski
\date 06/14/2024
\param[in] len    Side length
\param[in] height Height of the block in _vertices_
\param[in] i      i-index
\return Largest j-index for this i
*/
template <>
SPECTRUM_DEVICE_FUNC inline int BufferedBlock<4>::MaxFaceJ(int len, int height, int i) const
{
   return height - 1;
};

};

#endif
