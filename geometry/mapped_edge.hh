/*!
\file mapped_edge.hh
\brief Class representing a curve geometric element
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MAPPED_EDGE
#define SPECTRUM_MAPPED_EDGE

#ifdef GEO_DEBUG
#include <vector>
#include <iostream>
#include <iomanip>
#endif
#include "geometry/mapped_vertex.hh"

namespace Spectrum {

//          0-----------------------1     (order 1)
//          0-----------2-----------1     (order 2)
//          0-------2-------3-------1     (order 3)

//! Number of anchor nodes per edge
constexpr int n_edge_anchors[MAX_ELEMENT_ORDER] = {2, 3, 4};

//! Reference coordinates of the edge anchors (1 order)
constexpr double edge_anchors_1[n_edge_anchors[0]] = {0.0, 1.0};

//! Reference coordinates of the edge anchors (2 order)
constexpr double edge_anchors_2[n_edge_anchors[1]] = {0.0, 1.0, 0.5};

//! Reference coordinates of the edge anchors (3 order)
constexpr double edge_anchors_3[n_edge_anchors[2]] = {0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0};

//! Vertex anchor list (1 order)
constexpr int edge_anchors_list_1[2] = {0, 1};

//! Vertex anchor list (2 order)
constexpr int edge_anchors_list_2[2] = {0, 1};

//! Vertex anchor list (3 order)
constexpr int edge_anchors_list_3[2] = {0, 1};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MappedEdge class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A class describing a mapped edge on a curved line in 3D space
\author Vladimir Florinski
*/
template <int order>
class MappedEdge : public MappedElement<n_edge_anchors[order - 1]>
{
protected:

   using MappedElement<n_edge_anchors[order - 1]>::anchors;
   using MappedElement<n_edge_anchors[order - 1]>::ref_anchors;

//! Vertex anchor list (use pointer to save on memory)
   const int* vert_anchors_list;

//! Cosine of the central angle
   double c;

//! Compute the set of basis functions for the given point
   SPECTRUM_DEVICE_FUNC void BasisFunctions(double del, double* psi) const;

//! Compute the set of basis function derivatives for the given point
   SPECTRUM_DEVICE_FUNC void BasisDerivatives(double del, double* dpdd) const;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC MappedEdge(void);

//! Constructor with arguments
   SPECTRUM_DEVICE_FUNC MappedEdge(const GeoVector* new_anc);

//! Initialize the anchor set
   SPECTRUM_DEVICE_FUNC void SetAnchors(const GeoVector* new_anc);

//! Extract the edge anchors
   SPECTRUM_DEVICE_FUNC void GetVertAnchors(int iv, GeoVector* edge_anc) const;

//! Compute the position vector on the edge
   SPECTRUM_DEVICE_FUNC GeoVector Position(double del) const;

//! Compute the tangent vector with physical arc length on the edge
   SPECTRUM_DEVICE_FUNC GeoVector Tangent(double del) const;

#ifdef GEO_DEBUG

//! Tests if positions of the anchors are correct
   void TestAnchors(void) const;

//! Tests the coordinate mappings using exact values
   void TestCoordMapping(const std::vector<double>& del, const std::vector<GeoVector>& pos, const std::vector<GeoVector>& tan) const;

#endif

};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedEdge<1>::MappedEdge(void)
{
   ref_anchors = edge_anchors_1;
   vert_anchors_list = edge_anchors_list_1;
};

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedEdge<2>::MappedEdge(void)
{
   ref_anchors = edge_anchors_2;
   vert_anchors_list = edge_anchors_list_2;
};

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedEdge<3>::MappedEdge(void)
{
   ref_anchors = edge_anchors_3;
   vert_anchors_list = edge_anchors_list_3;
};

/*!
\author Vladimir Florinski
\date 01/06/2020
\param[in] new_anc New physical coordinates of the anchors
*/
template <int order>
SPECTRUM_DEVICE_FUNC inline MappedEdge<order>::MappedEdge(const GeoVector* new_anc)
                                             : MappedEdge<order>()
{
   SetAnchors(new_anc);
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in] new_anc New physical coordinates of the anchors
*/
template <int order>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<order>::SetAnchors(const GeoVector* new_anc)
{
   MappedElement<n_edge_anchors[order - 1]>::SetAnchors(new_anc);
   c = anchors[0] * anchors[1];
};

/*!
\author Vladimir Florinski
\date 04/03/2024
\param[in]  iv       Vertex index
\param[out] vert_anc Anchor sub-set for this vertex
*/
template <int order>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<order>::GetVertAnchors(int iv, GeoVector* vert_anc) const
{
   vert_anc[0] = anchors[vert_anchors_list[iv]];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Keyvan Ghanbari
\date 07/10/2019
\param[in]  del Reference coordinate along the edge
\param[out] psi All edge basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<1>::BasisFunctions(double del, double* psi) const
{
   psi[0] = 1.0 - del;
   psi[1] = del;
};

/*!
\author Vladimir Florinski
\author Keyvan Ghanbari
\date 07/10/2019
\param[in]  del Reference coordinate along the edge
\param[out] psi All edge basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<2>::BasisFunctions(double del, double* psi) const
{
   double omd = 1.0 - del;

   psi[0] = omd * (2.0 * omd - 1.0);
   psi[1] = del * (2.0 * del - 1.0);
   psi[2] = 4.0 * omd * del;
};

/*!
\author Vladimir Florinski
\author Keyvan Ghanbari
\date 07/10/2019
\param[in]  del Reference coordinate along the edge
\param[out] psi All edge basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<3>::BasisFunctions(double del, double* psi) const
{
   double omd = 1.0 - del;
   double tomdm1 = 3.0 * omd - 1.0;
   double tdelm1 = 3.0 * del - 1.0;
   
   psi[0] = -0.5 * omd * tomdm1 * tdelm1;
   psi[1] = -0.5 * del * tdelm1 * tomdm1;
   psi[2] =  4.5 * omd * del * tomdm1;
   psi[3] =  4.5 * del * omd * tdelm1;
};

/*!
\author Vladimir Florinski
\date 07/10/2019
\param[in]  del  Reference coordinate along the edge
\param[out] dpdd All edge basis function delta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<1>::BasisDerivatives(double del, double* dpdd) const
{
   dpdd[0] = -1.0;
   dpdd[1] =  1.0;
};

/*!
\author Vladimir Florinski
\date 07/10/2019
\param[in]  del  Reference coordinate along the edge
\param[out] dpdd All edge basis function delta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<2>::BasisDerivatives(double del, double* dpdd) const
{
   double omd = 1.0 - del;

   dpdd[0] = -4.0 * omd + 1.0;
   dpdd[1] =  4.0 * del - 1.0;
   dpdd[2] =  4.0 * (omd - del);
};

/*!
\author Vladimir Florinski
\date 07/10/2019
\param[in]  del  Reference coordinate along the edge
\param[out] dpdd All edge basis function delta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedEdge<3>::BasisDerivatives(double del, double* dpdd) const
{
   double omd = 1.0 - del;
   double tomdm1 = 3.0 * omd - 1.0;
   double tdelm1 = 3.0 * del - 1.0;

   dpdd[0] =  4.5 * omd * tdelm1 - 1.0;
   dpdd[1] = -4.5 * del * tomdm1 + 1.0;
   dpdd[2] =  4.5 * (omd * (1.0 - 9.0 * del) + 1.0);
   dpdd[3] = -4.5 * (del * (1.0 - 9.0 * omd) + 1.0);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] del Reference coordinate along the edge
\return Position vector
*/
template <int order>
SPECTRUM_DEVICE_FUNC inline GeoVector MappedEdge<order>::Position(double del) const
{
   double psi[n_edge_anchors[order - 1]];
   GeoVector pos = gv_zeros;

   BasisFunctions(del, psi);
   for(auto apt = 0; apt < n_edge_anchors[order - 1]; apt++) pos += psi[apt] * anchors[apt];
   return pos;
};

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] del Reference coordinate along the edge
\return Tangent vector with physical arc length
*/
template <int order>
SPECTRUM_DEVICE_FUNC inline GeoVector MappedEdge<order>::Tangent(double del) const
{
   double dpdd[n_edge_anchors[order - 1]];
   GeoVector tan = gv_zeros;

   BasisDerivatives(del, dpdd);
   for(auto apt = 0; apt < n_edge_anchors[order - 1]; apt++) tan += dpdd[apt] * anchors[apt];
   return tan;
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
template <int order>
inline void MappedEdge<order>::TestAnchors(void) const
{
   GeoVector error;

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Edge anchor position test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Anchor |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(auto apt = 0; apt < n_edge_anchors[order - 1]; apt++) {
      error = Position(ref_anchors[apt]) - anchors[apt];
      std::cerr << std::setw(5) << apt << "   |   "
                << std::setprecision(15)
                << std::fixed << std::setw(19) << anchors[apt][0]
                << std::fixed << std::setw(19) << anchors[apt][1]
                << std::fixed << std::setw(19) << anchors[apt][2] << "   |   "
                << std::fixed << std::setw(9) << std::setprecision(5) << error.Norm() / anchors[apt].Norm() * 100.0 << " %"
                << std::endl;
   };
   std::cerr << "-------------------------------------------------------------------------------------------\n";
};

/*!
\author Vladimir Florinski
\date 04/02/2024
\param[in] del Reference coordinate delta
\param[in] pos Exact positions of test points
\param[in] tan Exact tangent vectors _times length_ at test points
*/
template <int order>
inline void MappedEdge<order>::TestCoordMapping(const std::vector<double>& del, const std::vector<GeoVector>& pos,
                                                                                const std::vector<GeoVector>& tan) const
{
   int tpt;
   int n_test_points = del.size();
   GeoVector error;

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Edge position test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Point  |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(tpt = 0; tpt < n_test_points; tpt++) {
      error = Position(del[tpt]) - pos[tpt];
      std::cerr << std::setw(5) << tpt << "   |   "
                << std::setprecision(15)
                << std::fixed << std::setw(19) << pos[tpt][0]
                << std::fixed << std::setw(19) << pos[tpt][1]
                << std::fixed << std::setw(19) << pos[tpt][2] << "   |   "
                << std::fixed << std::setw(9) << std::setprecision(5) << error.Norm() / pos[tpt].Norm() * 100.0 << " %"
                << std::endl;
   };

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Edge tangent length test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Point  |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(tpt = 0; tpt < n_test_points; tpt++) {
      error = Tangent(del[tpt]) - tan[tpt];
      std::cerr << std::setw(5) << tpt << "   |   "
                << std::setprecision(15)
                << std::fixed << std::setw(19) << tan[tpt][0]
                << std::fixed << std::setw(19) << tan[tpt][1]
                << std::fixed << std::setw(19) << tan[tpt][2] << "   |   "
                << std::fixed << std::setw(9) << std::setprecision(5) << error.Norm() / tan[tpt].Norm() * 100.0 << " %"
                << std::endl;
   };
   std::cerr << "-------------------------------------------------------------------------------------------\n";
};

#endif

};

#endif
