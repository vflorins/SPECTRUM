/*!
\file stenciled_block.hh
\brief Extends the "GridBlock" class by adding stencils
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_STENCILED_BLOCK_HH
#define SPECTRUM_STENCILED_BLOCK_HH

#include <eigen3/Eigen/Dense>

#include <geodesic/grid_block.hh>

namespace Spectrum {

//! Stencil element (F)
#define GEOELM_STEN 0x0020

#ifdef GEO_DEBUG
//! Stencil name prefix string
const std::string stfile_pref = "stencil_";
#endif

/*!
\brief A class describing grid blocks with stencils
\author Vladimir Florinski

This class extends the basec GridBlock functionality by adding moments, stencils, and geometry matrices. Because of the use of exponential coordinates the zone aspect ratio is the same for all shells, so a single set of matrices is required per face. Stencil generation is based on a template, an idealized mesh made up of identical elements. The actual stencil is not always identical to the template because of singular corners.
*/
template <int verts_per_face>
class StenciledBlock : public GridBlock<verts_per_face>
{
protected:

   using PolygonalAddressing<verts_per_face>::square_fill;
   using SphericalSlab::ghost_height;
   using SphericalSlab::n_shells_withghost;
   using GeodesicSector<verts_per_face>::total_length;
   using GeodesicSector<verts_per_face>::ghost_width;
   using GeodesicSector<verts_per_face>::n_faces_withghost;
   using GeodesicSector<verts_per_face>::n_edges_withghost;
   using GeodesicSector<verts_per_face>::face_index_sector;
   using GeodesicSector<verts_per_face>::ev_local;
   using GeodesicSector<verts_per_face>::fv_local;
   using GeodesicSector<verts_per_face>::ff_local;
   using GeodesicSector<verts_per_face>::edge_mask;
   using GeodesicSector<verts_per_face>::face_mask;
   using GeodesicSector<verts_per_face>::MaxVertJ;
   using GeodesicSector<verts_per_face>::MaxFaceJ;
   using GridBlock<verts_per_face>::drp_ratio;
   using GridBlock<verts_per_face>::block_vert_cart;

//! A local Polynomial object
//   Polynomial poly;

//! The total number of stencils
   static constexpr int n_stencils = 1 + 2 * verts_per_face;

//! Number of stenciled faces
//   int n_faces_stenciled;

//! Number of zones in each stencil, _excluding_ the principal
   int zones_per_stencil[n_stencils];

//! Face areas on the US
   double* face_area = nullptr;

//! Face centers of mass for the US (do _not_ lie on the US)
   GeoVector* face_cmass = nullptr;

//! Edge lengths on the US
   double* edge_length = nullptr;

//! Edge centers of mass for the US (do _not_ lie on the US)
   GeoVector* edge_cmass = nullptr;

//! Stencil zonelists
   int*** stencil_zonelist = nullptr;

//! Transposed geometry matrices
   Eigen::MatrixXd** geom_matr_At = nullptr;

//! LU decomposition of At*A
   Eigen::PartialPivLU<Eigen::MatrixXd>** geom_matr_LU = {nullptr};

//! Build a list of stenciled zones and compute the stencil mask
   void MarkStenciledArea(void);

//! Maximum value of the vertex j-index for the given vertex i-index for a truncated block
   int MaxVertJ(int len, int height, int i) const;

//! Maximum value of the face j-index for the given face i-index for a truncated block
   int MaxFaceJ(int len, int height, int i) const;

//! Determine whether a zone is in the block's interior - ijk version
//   bool IsInteriorZoneOfBlock(int i, int j, int k) const;

//! Determine whether a zone is in the block's interior - face+k version
//   bool IsInteriorZoneOfBlock(int face, int k) const;

//! Compute the geometry matrix for one stencil of a single face
   void ComputeOneMatrix(int pface, int stencil);

//! Compute the zone lists for each stencil for each zone
   void BuildAllStencils(void);

//! Compute the volumetric moments
   void ComputeMoments(void);

public:

//! Default constructor
   StenciledBlock(void);

//! Copy constructor
   StenciledBlock(const StenciledBlock& other);

//! Move constructor
   StenciledBlock(StenciledBlock&& other) noexcept;

//! Constructor with arguments
   StenciledBlock(int width, int wghost, int height, int hghost);

//! Destructor
   ~StenciledBlock();

//! Allocate memory
   void SetDimensions(int width, int wghost, int height, int hghost, bool construct);

//! Free all dynamically allocated memory
   void FreeStorage(void);

//! Set up the dimensions and geometry of the mesh
   void AssociateMesh(double ximin, double ximax, const bool* corners, const bool* borders,
                      const GeoVector* vcart, std::shared_ptr<DistanceBase> dist_map_in, bool construct);

#ifdef GEO_DEBUG

//! Print the properties of a stencil
   void PrintStencilProps(int pface, int stencil) const;

//! Draw stencils for a given shell and face
   void DrawStencil(int pshell, int pface, int stencil, double rot_z, double rot_x) const;

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
inline int StenciledBlock<3>::MaxVertJ(int len, int height, int i) const
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
inline int StenciledBlock<4>::MaxVertJ(int len, int height, int i) const
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
inline int StenciledBlock<3>::MaxFaceJ(int len, int height, int i) const
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
inline int StenciledBlock<4>::MaxFaceJ(int len, int height, int i) const
{
   return height - 1;
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] i First index
\param[in] j Second index
\param[in] k The r-shell index
\return True if the zone is interior to the block
*/
/*
template <int verts_per_face>
inline bool StenciledBlock<verts_per_face>::IsInteriorZoneOfBlock(int i, int j, int k) const
{
   return SphericalSlab::IsInteriorShellOfSlab(k) && GeodesicSector<verts_per_face>::IsInteriorFaceOfSector(i, j);
};
*/
/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] face The face index
\param[in] k    The r-shell index
\return True if the zone is interior to the block
*/
/*
template <int verts_per_face>
inline bool StenciledBlock<verts_per_face>::IsInteriorZoneOfBlock(int face, int k) const
{
   return SphericalSlab::IsInteriorShellOfSlab(k) && GridBlock<verts_per_face>::IsInteriorFaceOfSector(face);
};
*/
};

#endif

