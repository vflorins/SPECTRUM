/*!
\file mapped_face.hh
\brief Class representing a surface geometric element
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MAPPED_FACE
#define SPECTRUM_MAPPED_FACE

#ifdef GEO_DEBUG
#include <vector>
#include <iostream>
#include <iomanip>
#endif
#include "geometry/mapped_edge.hh"

namespace Spectrum {

//! Area of the reference (plane) element
constexpr double ref_face_area[2] = {sqrtthr / 4.0, 1.0};

//                      2                    3-----------------------2
//                     / \                   |                       |
//                    /   \                  |          E2           |
//                   /     \                 |                       |
//                  /       \                |                       |
//                 /         \               |                       |
//                /  E2   E1  \              |   E3             E1   |     (order 1)
//               /             \             |                       |
//              /               \            |                       |
//             /                 \           |                       |
//            /        E0         \          |          E0           |
//           /                     \         |                       |
//          0-----------------------1        0-----------------------1
//          0-----------------------1
//
//                      2                    3-----------6-----------2
//                     / \                   |                       |
//                    /   \                  |          E2           |
//                   /     \                 |                       |
//                  /       \                |                       |
//                 /         \               |                       |
//                5  E2   E1  4              7   E3             E1   5     (order 2)
//               /             \             |                       |
//              /               \            |                       |
//             /                 \           |                       |
//            /        E0         \          |          E0           |
//           /                     \         |                       |
//          0-----------3-----------1        0-----------4-----------1
//          0-----------2-----------1
//
//                      2                    3-------9-------8-------2
//                     / \                   |                       |
//                    /   \                  |          E2           |
//                   /     \                 |                       |
//                  7       6               10                       7
//                 /         \               |                       |
//                /  E2   E1  \              |   E3             E1   |     (order 3)
//               /             \             |                       |
//              8       9       5           11                       6
//             /                 \           |                       |
//            /        E0         \          |          E0           |
//           /                     \         |                       |
//          0-------3-------4-------1        0-------4-------5-------1
//          0-------2-------3-------1

//! Number of anchor nodes per face
constexpr int n_face_anchors[2][MAX_ELEMENT_ORDER] = {{3, 6, 10}, {4, 8, 12}};

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Reference coordinates of the face anchors (triangles, 1 order)
constexpr double face_anchors_tria_1[n_face_anchors[0][0]][2] = {{0.0,       0.0          },
                                                                 {1.0,       0.0          },
                                                                 {1.0 / 2.0, sqrtthr / 2.0}};

//! Reference coordinates of the face anchors (triangles, 2 order)
constexpr double face_anchors_tria_2[n_face_anchors[0][1]][2] = {{0.0,       0.0          },
                                                                 {1.0,       0.0          },
                                                                 {1.0 / 2.0, sqrtthr / 2.0},
                                                                 {1.0 / 2.0, 0.0          },
                                                                 {3.0 / 4.0, sqrtthr / 4.0},
                                                                 {1.0 / 4.0, sqrtthr / 4.0}};

//! Reference coordinates of the face anchors (triangles, 3 order)
constexpr double face_anchors_tria_3[n_face_anchors[0][2]][2] = {{0.0,       0.0                  },
                                                                 {1.0,       0.0                  },
                                                                 {1.0 / 2.0, sqrtthr / 2.0        },
                                                                 {1.0 / 3.0, 0.0                  },
                                                                 {2.0 / 3.0, 0.0                  },
                                                                 {5.0 / 6.0, 1.0 / (2.0 * sqrtthr)},
                                                                 {2.0 / 3.0, 1.0 / sqrtthr        },
                                                                 {1.0 / 3.0, 1.0 / sqrtthr        },
                                                                 {1.0 / 6.0, 1.0 / (2.0 * sqrtthr)},
                                                                 {1.0 / 2.0, 1.0 / (2.0 * sqrtthr)}};

//! Reference coordinates of the face anchors (quads, 1 order)
constexpr double face_anchors_quad_1[n_face_anchors[1][0]][2] = {{0.0, 0.0},
                                                                 {1.0, 0.0},
                                                                 {1.0, 1.0},
                                                                 {0.0, 1.0}};

//! Reference coordinates of the face anchors (quads, 2 order)
constexpr double face_anchors_quad_2[n_face_anchors[1][1]][2] = {{0.0,       0.0      },
                                                                 {1.0,       0.0      },
                                                                 {1.0,       1.0      },
                                                                 {0.0,       1.0      },
                                                                 {1.0 / 2.0, 0.0      },
                                                                 {1.0,       1.0 / 2.0},
                                                                 {1.0 / 2.0, 1.0      },
                                                                 {0.0,       1.0 / 2.0}};

//! Reference coordinates of the face anchors (quads, 3 order)
constexpr double face_anchors_quad_3[n_face_anchors[1][2]][2] = {{0.0,       0.0      },
                                                                 {1.0,       0.0      },
                                                                 {1.0,       1.0      },
                                                                 {0.0,       1.0      },
                                                                 {1.0 / 3.0, 0.0      },
                                                                 {2.0 / 3.0, 0.0      },
                                                                 {1.0,       1.0 / 3.0},
                                                                 {1.0,       2.0 / 3.0},
                                                                 {2.0 / 3.0, 1.0      },
                                                                 {1.0 / 3.0, 1.0      },
                                                                 {0.0,       2.0 / 3.0},
                                                                 {0.0,       1.0 / 3.0}};

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Edge anchor list (triangles, 1 order)
constexpr int edge_anchors_list_tria_1[3][n_edge_anchors[0]] = {{0, 1},
                                                                {1, 2},
                                                                {2, 0}};

//! Edge anchor list (triangles, 2 order)
constexpr int edge_anchors_list_tria_2[3][n_edge_anchors[1]] = {{0, 1, 3},
                                                                {1, 2, 4},
                                                                {2, 0, 5}};

//! Edge anchor list (triangles, 3 order)
constexpr int edge_anchors_list_tria_3[3][n_edge_anchors[2]] = {{0, 1, 3, 4},
                                                                {1, 2, 5, 6},
                                                                {2, 0, 7, 8}};

//! Edge anchor list (quads, 1 order)
constexpr int edge_anchors_list_quad_1[4][n_edge_anchors[0]] = {{0, 1},
                                                                {1, 2},
                                                                {2, 3},
                                                                {3, 0}};

//! Edge anchor list (quads, 2 order)
constexpr int edge_anchors_list_quad_2[4][n_edge_anchors[1]] = {{0, 1, 4},
                                                                {1, 2, 5},
                                                                {2, 3, 6},
                                                                {3, 0, 7}};

//! Edge anchor list (quads, 3 order)
constexpr int edge_anchors_list_quad_3[4][n_edge_anchors[2]] = {{0, 1,  4,  5},
                                                                {1, 2,  6,  7},
                                                                {2, 3,  8,  9},
                                                                {3, 0, 10, 11}};

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Edge coordinates of the edge anchors (triangles, 1 order)
constexpr double edge_anchors_edge_tria_1[3][n_edge_anchors[0]] = {{0.0, 1.0},
                                                                   {0.0, 1.0},
                                                                   {0.0, 1.0}};

//! Edge coordinates of the edge anchors (triangles, 2 order)
constexpr double edge_anchors_edge_tria_2[3][n_edge_anchors[1]] = {{0.0, 1.0, 1.0 / 2.0},
                                                                   {0.0, 1.0, 1.0 / 2.0},
                                                                   {0.0, 1.0, 1.0 / 2.0}};

//! Edge coordinates of the edge anchors (triangles, 3 order)
constexpr double edge_anchors_edge_tria_3[3][n_edge_anchors[2]] = {{0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0},
                                                                   {0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0},
                                                                   {0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0}};

//! Edge coordinates of the edge anchors (quads, 1 order)
constexpr double edge_anchors_edge_quad_1[4][n_edge_anchors[0]] = {{0.0, 1.0},
                                                                   {0.0, 1.0},
                                                                   {0.0, 1.0},
                                                                   {0.0, 1.0}};

//! Edge coordinates of the edge anchors (quads, 2 order)
constexpr double edge_anchors_edge_quad_2[4][n_edge_anchors[1]] = {{0.0, 1.0, 1.0 / 2.0},
                                                                   {0.0, 1.0, 1.0 / 2.0},
                                                                   {0.0, 1.0, 1.0 / 2.0},
                                                                   {0.0, 1.0, 1.0 / 2.0}};

//! Edge coordinates of the edge anchors (quads, 3 order)
constexpr double edge_anchors_edge_quad_3[4][n_edge_anchors[2]] = {{0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0},
                                                                   {0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0},
                                                                   {0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0},
                                                                   {0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0}};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MappedFace class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A class describing a mapped face on a curved 2D surface in 3D space
\author Vladimir Florinski
*/
template <int n_verts, int order>
class MappedFace : public MappedElement<n_face_anchors[n_verts - 3][order - 1]>
{
protected:

   using MappedElement<n_face_anchors[n_verts - 3][order - 1]>::anchors;
   using MappedElement<n_face_anchors[n_verts - 3][order - 1]>::ref_anchors;

//! Edge anchor list (use pointer to save on memory)
   const int* edge_anchors_list;

//! Edge coordinates of the edge anchors (use pointer to save on memory)
   const double* edge_anchors_edge;

//! Cosines of the central angles
   double c[n_verts];

//! Compute the set of basis functions for the given point
   SPECTRUM_DEVICE_FUNC void BasisFunctions(double alp, double bet, double* psi) const;

//! Compute the complete set of basis function derivatives for the given point
   SPECTRUM_DEVICE_FUNC void BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC MappedFace(void);

//! Constructor with arguments
   SPECTRUM_DEVICE_FUNC MappedFace(const GeoVector* new_anc);

//! Initialize the anchor set
   SPECTRUM_DEVICE_FUNC void SetAnchors(const GeoVector* new_anc);

//! Extract the edge anchors
   SPECTRUM_DEVICE_FUNC void GetEdgeAnchors(int ie, GeoVector* edge_anc) const;

//! Compute the position vector on the face
   SPECTRUM_DEVICE_FUNC GeoVector Position(double alp, double bet) const;

//! Compute the first tangent vector with physical length factor on the face
   SPECTRUM_DEVICE_FUNC GeoVector TangentA(double alp, double bet) const;

//! Compute the second tangent vector with physical length factor on the face
   SPECTRUM_DEVICE_FUNC GeoVector TangentB(double alp, double bet) const;

//! Compute the normal vector with physical area factor on the face
   SPECTRUM_DEVICE_FUNC GeoVector Normal(double alp, double bet) const;

#ifdef GEO_DEBUG

//! Tests if positions of the anchors are correct
   void TestAnchors(void) const;

//! Tests the coordinate mappings using exact values
   void TestCoordMapping(const std::vector<double>& alp, const std::vector<double>& bet,
                         const std::vector<GeoVector>& pos, const std::vector<GeoVector>& nrm) const;

#endif

};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedFace<3, 1>::MappedFace(void)
{
   ref_anchors = face_anchors_tria_1[0];
   edge_anchors_list = edge_anchors_list_tria_1[0];
   edge_anchors_edge = edge_anchors_edge_tria_1[0];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedFace<3, 2>::MappedFace(void)
{
   ref_anchors = face_anchors_tria_2[0];
   edge_anchors_list = edge_anchors_list_tria_2[0];
   edge_anchors_edge = edge_anchors_edge_tria_2[0];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedFace<3, 3>::MappedFace(void)
{
   ref_anchors = &face_anchors_tria_3[0][0];
   edge_anchors_list = &edge_anchors_list_tria_3[0][0];
   edge_anchors_edge = &edge_anchors_edge_tria_3[0][0];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedFace<4, 1>::MappedFace(void)
{
   ref_anchors = face_anchors_quad_1[0];
   edge_anchors_list = edge_anchors_list_quad_1[0];
   edge_anchors_edge = edge_anchors_edge_quad_1[0];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedFace<4, 2>::MappedFace(void)
{
   ref_anchors = face_anchors_quad_2[0];
   edge_anchors_list = edge_anchors_list_quad_2[0];
   edge_anchors_edge = edge_anchors_edge_quad_2[0];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline MappedFace<4, 3>::MappedFace(void)
{
   ref_anchors = face_anchors_quad_3[0];
   edge_anchors_list = edge_anchors_list_quad_3[0];
   edge_anchors_edge = edge_anchors_edge_quad_3[0];
};

/*!
\author Vladimir Florinski
\date 01/06/2020
\param[in] new_anc New physical coordinates of the anchors
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline MappedFace<n_verts, order>::MappedFace(const GeoVector* new_anc)
                                                      : MappedFace<n_verts, order>()
{
   SetAnchors(new_anc);
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in] new_anc New physical coordinates of the anchors
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline void MappedFace<n_verts, order>::SetAnchors(const GeoVector* new_anc)
{
   MappedElement<n_face_anchors[n_verts - 3][order - 1]>::SetAnchors(new_anc);
   for(auto ie = 0; ie < n_verts; ie++) {
      c[ie] = anchors[edge_anchors_list[ie * n_edge_anchors[order - 1]]] * anchors[edge_anchors_list[ie * n_edge_anchors[order - 1] + 1]];
   };
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  ie       Edge index
\param[out] edge_anc Anchor sub-set for this edge
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline void MappedFace<n_verts, order>::GetEdgeAnchors(int ie, GeoVector* edge_anc) const
{
   for(auto apt = 0; apt < n_edge_anchors[order - 1]; apt++) {
      edge_anc[apt] = anchors[edge_anchors_list[ie * n_edge_anchors[order - 1] + apt]];
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp First reference coordinate
\param[in]  bet Second reference coordinate
\param[out] psi All face basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<3, 1>::BasisFunctions(double alp, double bet, double* psi) const
{
   psi[0] = 1.0 - alp - oo_sqrtthr * bet;
   psi[1] = alp - oo_sqrtthr * bet;
   psi[2] = 2.0 * oo_sqrtthr * bet;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp First reference coordinate
\param[in]  bet Second reference coordinate
\param[out] psi All face basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<3, 2>::BasisFunctions(double alp, double bet, double* psi) const
{
// Natural coordinates
   double L[3];

   L[0] = 1.0 - alp - oo_sqrtthr * bet;
   L[1] = alp - oo_sqrtthr * bet;
   L[2] = 2.0 * oo_sqrtthr * bet;

   psi[0] = L[0] * (2.0 * L[0] - 1.0);
   psi[1] = L[1] * (2.0 * L[1] - 1.0);
   psi[2] = L[2] * (2.0 * L[2] - 1.0);

   psi[3] = 4.0 * L[0] * L[1];
   psi[4] = 4.0 * L[1] * L[2];
   psi[5] = 4.0 * L[2] * L[0];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp First reference coordinate
\param[in]  bet Second reference coordinate
\param[out] psi All face basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<3, 3>::BasisFunctions(double alp, double bet, double* psi) const
{
// Natural coordinates
   double L[3], tl1[3];

   L[0] = 1.0 - alp - oo_sqrtthr * bet;
   L[1] = alp - oo_sqrtthr * bet;
   L[2] = 2.0 * oo_sqrtthr * bet;

   tl1[0] = 3.0 * L[0] - 1.0;
   tl1[1] = 3.0 * L[1] - 1.0;
   tl1[2] = 3.0 * L[2] - 1.0;

   psi[0] = 0.5 * L[0] * tl1[0] * (tl1[0] - 1.0);
   psi[1] = 0.5 * L[1] * tl1[1] * (tl1[1] - 1.0);
   psi[2] = 0.5 * L[2] * tl1[2] * (tl1[2] - 1.0);

   psi[3] = 4.5 * L[0] * L[1] * tl1[0];
   psi[4] = 4.5 * L[0] * L[1] * tl1[1];
   psi[5] = 4.5 * L[1] * L[2] * tl1[1];
   psi[6] = 4.5 * L[1] * L[2] * tl1[2];
   psi[7] = 4.5 * L[2] * L[0] * tl1[2];
   psi[8] = 4.5 * L[2] * L[0] * tl1[0];

   psi[9] = 27.0 * L[0] * L[1] * L[2];
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp First reference coordinate
\param[in]  bet Second reference coordinate
\param[out] psi All face basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<4, 1>::BasisFunctions(double alp, double bet, double* psi) const
{
   double oma = 1.0 - alp;
   double omb = 1.0 - bet;

   psi[0] = oma * omb;
   psi[1] = alp * omb;
   psi[2] = alp * bet;
   psi[3] = bet * oma;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp First reference coordinate
\param[in]  bet Second reference coordinate
\param[out] psi All face basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<4, 2>::BasisFunctions(double alp, double bet, double* psi) const
{
   double oma = 1.0 - alp;
   double omb = 1.0 - bet;

   psi[0] =  oma * omb * (1.0 - 2.0 * (alp + bet));
   psi[1] = -alp * omb * (1.0 - 2.0 * (alp - bet));
   psi[2] = -alp * bet * (3.0 - 2.0 * (alp + bet));
   psi[3] = -bet * oma * (1.0 - 2.0 * (bet - alp));

   psi[4] = 4.0 * alp * oma * omb;
   psi[5] = 4.0 * alp * bet * omb;
   psi[6] = 4.0 * bet * alp * oma;
   psi[7] = 4.0 * bet * omb * oma;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp First reference coordinate
\param[in]  bet Second reference coordinate
\param[out] psi All face basis function values
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<4, 3>::BasisFunctions(double alp, double bet, double* psi) const
{
   double oma = 1.0 - alp;
   double omb = 1.0 - bet;
   double tab = 2.0 - 9.0 * alp * oma - 9.0 * bet * omb;
   double ota = 1.0 - 3.0 * alp;
   double otb = 1.0 - 3.0 * bet;
   double tta = 2.0 - 3.0 * alp;
   double ttb = 2.0 - 3.0 * bet;

   psi[0] = 0.5 * oma * omb * tab;
   psi[1] = 0.5 * alp * omb * tab;
   psi[2] = 0.5 * alp * bet * tab;
   psi[3] = 0.5 * bet * oma * tab;

   psi[4]  =  4.5 * alp * oma * omb * tta;
   psi[5]  = -4.5 * alp * oma * omb * ota;
   psi[6]  =  4.5 * alp * bet * omb * ttb;
   psi[7]  = -4.5 * alp * bet * omb * otb;
   psi[8]  = -4.5 * bet * alp * oma * ota;
   psi[9]  =  4.5 * bet * alp * oma * tta;
   psi[10] = -4.5 * bet * omb * oma * otb;
   psi[11] =  4.5 * bet * omb * oma * ttb;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp  First reference coordinate
\param[in]  bet  Second reference coordinate
\param[out] dpda All face basis function alpha derivatives
\param[out] dpdb All face basis function beta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<3, 1>::BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const
{
   dpda[0] = -1.0;
   dpda[1] =  1.0;
   dpda[2] =  0.0;

   dpdb[0] = -oo_sqrtthr;
   dpdb[1] = -oo_sqrtthr;
   dpdb[2] = 2.0 * oo_sqrtthr;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp  First reference coordinate
\param[in]  bet  Second reference coordinate
\param[out] dpda All face basis function alpha derivatives
\param[out] dpdb All face basis function beta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<3, 2>::BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const
{
// Natural coordinates
   double L[3];

   L[0] = 1.0 - alp - oo_sqrtthr * bet;
   L[1] = alp - oo_sqrtthr * bet;
   L[2] = 2.0 * oo_sqrtthr * bet;

   dpda[0] = -4.0 * L[0] + 1.0;
   dpda[1] =  4.0 * L[1] - 1.0;
   dpda[2] =  0.0;
   dpda[3] =  4.0 * (L[0] - L[1]);
   dpda[4] =  4.0 * L[2];
   dpda[5] = -4.0 * L[2];

   dpdb[0] =        oo_sqrtthr * dpda[0];
   dpdb[1] =       -oo_sqrtthr * dpda[1];
   dpdb[2] =  2.0 * oo_sqrtthr * (4.0 * L[2] - 1.0);
   dpdb[3] = -4.0 * oo_sqrtthr * (L[0] + L[1]);
   dpdb[4] =  4.0 * oo_sqrtthr * (2.0 * L[1] - L[2]);
   dpdb[5] =  4.0 * oo_sqrtthr * (2.0 * L[0] - L[2]);
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp  First reference coordinate
\param[in]  bet  Second reference coordinate
\param[out] dpda All face basis function alpha derivatives
\param[out] dpdb All face basis function beta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<3, 3>::BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const
{
// Natural coordinates
   double L[3], tl1[3], sl1[3];

   L[0] = 1.0 - alp - oo_sqrtthr * bet;
   L[1] = alp - oo_sqrtthr * bet;
   L[2] = 2.0 * oo_sqrtthr * bet;

   tl1[0] = 3.0 * L[0] - 1.0;
   tl1[1] = 3.0 * L[1] - 1.0;
   tl1[2] = 3.0 * L[2] - 1.0;

   sl1[0] = 6.0 * L[0] - 1.0;
   sl1[1] = 6.0 * L[1] - 1.0;
   sl1[2] = 6.0 * L[2] - 1.0;

   dpda[0] = -9.0 * L[0] * (1.5 * L[0] - 1.0) - 1.0;
   dpda[1] =  9.0 * L[1] * (1.5 * L[1] - 1.0) + 1.0;
   dpda[2] =  0.0;
   dpda[3] =  4.5 * (L[0] * tl1[0] - L[1] * sl1[0]);
   dpda[4] = -4.5 * (L[1] * tl1[1] - L[0] * sl1[1]);
   dpda[5] =  4.5 * L[2] * sl1[1];
   dpda[6] =  4.5 * L[2] * tl1[2];
   dpda[7] = -4.5 * L[2] * tl1[2];
   dpda[8] = -4.5 * L[2] * sl1[0];
   dpda[9] = 27.0 * L[2] * (L[0] - L[1]);

   dpdb[0] =        oo_sqrtthr * dpda[0];
   dpdb[1] =       -oo_sqrtthr * dpda[1];
   dpdb[2] =  2.0 * oo_sqrtthr * (9.0 * L[2] * (1.5 * L[2] - 1.0) + 1.0);
   dpdb[3] = -4.5 * oo_sqrtthr * (L[0] * tl1[0] + L[1] * sl1[0]);
   dpdb[4] = -4.5 * oo_sqrtthr * (L[1] * tl1[1] + L[0] * sl1[1]);
   dpdb[5] =  4.5 * oo_sqrtthr * (2.0 * L[1] * tl1[1] - L[2] * sl1[1]);
   dpdb[6] = -4.5 * oo_sqrtthr * (L[2] * tl1[2] - 2.0 * L[1] * sl1[2]);
   dpdb[7] = -4.5 * oo_sqrtthr * (L[2] * tl1[2] - 2.0 * L[0] * sl1[2]);   
   dpdb[8] =  4.5 * oo_sqrtthr * (2.0 * L[0] * tl1[0] - L[2] * sl1[0]);
   dpdb[9] = 27.0 * oo_sqrtthr * (2.0 * L[0] * L[1] - L[2] * (L[0] + L[1]));
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp  First reference coordinate
\param[in]  bet  Second reference coordinate
\param[out] dpda All face basis function alpha derivatives
\param[out] dpdb All face basis function beta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<4, 1>::BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const
{
   double oma = 1.0 - alp;
   double omb = 1.0 - bet;

   dpda[0] = -omb;
   dpda[1] =  omb;
   dpda[2] =  bet;
   dpda[3] = -bet;

   dpdb[0] = -oma;
   dpdb[1] = -alp;
   dpdb[2] =  alp;
   dpdb[3] =  oma;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp  First reference coordinate
\param[in]  bet  Second reference coordinate
\param[out] dpda All face basis function alpha derivatives
\param[out] dpdb All face basis function beta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<4, 2>::BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const
{
   double oma = 1.0 - alp;
   double omb = 1.0 - bet;

   dpda[0] = -omb * (3.0 - 4.0 * alp - 2.0 * bet);
   dpda[1] = -omb * (1.0 - 4.0 * alp + 2.0 * bet);
   dpda[2] = -bet * (3.0 - 4.0 * alp - 2.0 * bet);
   dpda[3] = -bet * (1.0 - 4.0 * alp + 2.0 * bet);
   dpda[4] =  4.0 * omb * (1.0 - 2.0 * alp);
   dpda[5] =  4.0 * bet * omb;
   dpda[6] =  4.0 * bet * (1.0 - 2.0 * alp);
   dpda[7] = -4.0 * bet * omb;

   dpdb[0] = -oma * (3.0 - 4.0 * bet - 2.0 * alp);
   dpdb[1] = -alp * (1.0 - 4.0 * bet + 2.0 * alp);
   dpdb[2] = -alp * (3.0 - 4.0 * bet - 2.0 * alp);
   dpdb[3] = -oma * (1.0 - 4.0 * bet + 2.0 * alp);
   dpdb[4] = -4.0 * alp * oma;
   dpdb[5] =  4.0 * alp * (1.0 - 2.0 * bet);
   dpdb[6] =  4.0 * alp * oma;
   dpdb[7] =  4.0 * oma * (1.0 - 2.0 * bet);
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in]  alp  First reference coordinate
\param[in]  bet  Second reference coordinate
\param[out] dpda All face basis function alpha derivatives
\param[out] dpdb All face basis function beta derivatives
*/
template <>
SPECTRUM_DEVICE_FUNC inline void MappedFace<4, 3>::BasisDerivatives(double alp, double bet, double* dpda, double* dpdb) const
{
   double oma = 1.0 - alp;
   double omb = 1.0 - bet;
   double ota = 1.0 - 3.0 * alp;
   double otb = 1.0 - 3.0 * bet;
   double tta = 2.0 - 3.0 * alp;
   double ttb = 2.0 - 3.0 * bet;

   dpda[0]  = -0.5 * omb * (2.0 + 9.0 * oma * ota - 9.0 * bet * omb);
   dpda[1]  =  0.5 * omb * (2.0 - 9.0 * alp * tta - 9.0 * bet * omb);
   dpda[2]  =  0.5 * bet * (2.0 - 9.0 * alp * tta - 9.0 * bet * omb);
   dpda[3]  = -0.5 * bet * (2.0 + 9.0 * oma * ota - 9.0 * bet * omb);
   dpda[4]  =  4.5 * omb * (2.0 - 9.0 * alp * oma - alp);
   dpda[5]  = -4.5 * omb * (1.0 - 9.0 * alp * oma + alp);
   dpda[6]  =  4.5 * omb * bet * ttb;
   dpda[7]  = -4.5 * omb * bet * otb;
   dpda[8]  = -4.5 * bet * (1.0 - 9.0 * alp * oma + alp);
   dpda[9]  =  4.5 * bet * (2.0 - 9.0 * alp * oma - alp);
   dpda[10] =  4.5 * bet * omb * otb;
   dpda[11] = -4.5 * bet * omb * ttb;

   dpdb[0]  = -0.5 * oma * (2.0 + 9.0 * omb * otb - 9.0 * alp * oma);
   dpdb[1]  = -0.5 * alp * (2.0 + 9.0 * omb * otb - 9.0 * alp * oma);
   dpdb[2]  =  0.5 * alp * (2.0 - 9.0 * bet * ttb - 9.0 * alp * oma);
   dpdb[3]  =  0.5 * oma * (2.0 - 9.0 * bet * ttb - 9.0 * alp * oma);
   dpdb[4]  = -4.5 * alp * oma * tta;
   dpdb[5]  =  4.5 * alp * oma * ota;
   dpdb[6]  =  4.5 * alp * (2.0 - 9.0 * bet * omb - bet);
   dpdb[7]  = -4.5 * alp * (1.0 - 9.0 * bet * omb + bet);
   dpdb[8]  = -4.5 * oma * alp * ota;
   dpdb[9]  =  4.5 * oma * alp * tta;
   dpdb[10] = -4.5 * oma * (1.0 - 9.0 * bet * omb + bet);
   dpdb[11] =  4.5 * oma * (2.0 - 9.0 * bet * omb - bet);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] alp First reference coordinate
\param[in] bet Second reference coordinate
\return Position vector
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline GeoVector MappedFace<n_verts, order>::Position(double alp, double bet) const
{
   double psi[n_face_anchors[n_verts - 3][order - 1]];
   GeoVector pos = gv_zeros;

   BasisFunctions(alp, bet, psi);
   for(auto apt = 0; apt < n_face_anchors[n_verts - 3][order - 1]; apt++) pos += psi[apt] * anchors[apt];
   return pos;
};

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] alp First reference coordinate
\param[in] bet Second reference coordinate
\return Tangent vector in the first coordinate with physical length factor
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline GeoVector MappedFace<n_verts, order>::TangentA(double alp, double bet) const
{
   double dpda[n_face_anchors[n_verts - 3][order - 1]], dpdb[n_face_anchors[n_verts - 3][order - 1]];
   GeoVector tana = gv_zeros;

   BasisDerivatives(alp, bet, dpda, dpdb);
   for(auto apt = 0; apt < n_face_anchors[n_verts - 3][order - 1]; apt++) tana += dpda[apt] * anchors[apt];

// To ensure proper normalization the physical area (i.e., the Jacobian) must include the area of the reference element.
   return sqrt(ref_face_area[n_verts - 3]) * tana;
};

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] alp First reference coordinate
\param[in] bet Second reference coordinate
\return Tangent vector in the second coordinate with physical length factor
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline GeoVector MappedFace<n_verts, order>::TangentB(double alp, double bet) const
{
   double dpda[n_face_anchors[n_verts - 3][order - 1]], dpdb[n_face_anchors[n_verts - 3][order - 1]];
   GeoVector tanb = gv_zeros;

   BasisDerivatives(alp, bet, dpda, dpdb);
   for(auto apt = 0; apt < n_face_anchors[n_verts - 3][order - 1]; apt++) tanb += dpdb[apt] * anchors[apt];

// To ensure proper normalization the physical area (i.e., the Jacobian) must include the area of the reference element.
   return sqrt(ref_face_area[n_verts - 3]) * tanb;
};

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] alp First reference coordinate
\param[in] bet Second reference coordinate
\return Normal vector with physical area factor
*/
template <int n_verts, int order>
SPECTRUM_DEVICE_FUNC inline GeoVector MappedFace<n_verts, order>::Normal(double alp, double bet) const
{
   double dpda[n_face_anchors[n_verts - 3][order - 1]], dpdb[n_face_anchors[n_verts - 3][order - 1]];
   GeoVector tana = gv_zeros, tanb = gv_zeros;

   BasisDerivatives(alp, bet, dpda, dpdb);
   for(auto apt = 0; apt < n_face_anchors[n_verts - 3][order - 1]; apt++) {
      tana += dpda[apt] * anchors[apt];
      tanb += dpdb[apt] * anchors[apt];
   };

// To ensure proper normalization the physical area (i.e., the Jacobian) must include the area of the reference element.
   return ref_face_area[n_verts - 3] * (tana ^ tanb);
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
template <int n_verts, int order>
inline void MappedFace<n_verts, order>::TestAnchors(void) const
{
   GeoVector error;

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Face anchor position test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Anchor |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(auto apt = 0; apt < n_face_anchors[n_verts - 3][order - 1]; apt++) {
      error = Position(ref_anchors[2 * apt], ref_anchors[2 * apt + 1]) - anchors[apt];
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
\param[in] face_ref Reference coordinates of test points (alpha, beta)
\param[in] face_pos Exact positions of facial test points
\param[in] face_nrm Exact normal vectors _times area_ at facial test points
*/
template <int n_verts, int order>
inline void MappedFace<n_verts, order>::TestCoordMapping(const std::vector<double>& alp, const std::vector<double>& bet,
                                                         const std::vector<GeoVector>& pos, const std::vector<GeoVector>& nrm) const
{
   int tpt;
   int n_test_points = alp.size();
   GeoVector error;

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Face position test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Point  |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(tpt = 0; tpt < n_test_points; tpt++) {
      error = Position(alp[tpt], bet[tpt]) - pos[tpt];
      std::cerr << std::setw(5) << tpt << "   |   "
                << std::setprecision(15)
                << std::fixed << std::setw(19) << pos[tpt][0]
                << std::fixed << std::setw(19) << pos[tpt][1]
                << std::fixed << std::setw(19) << pos[tpt][2] << "   |   "
                << std::fixed << std::setw(9) << std::setprecision(5) << error.Norm() / pos[tpt].Norm() * 100.0 << " %"
                << std::endl;
   };

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Face normal area test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Point  |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(tpt = 0; tpt < n_test_points; tpt++) {
      error = Normal(alp[tpt], bet[tpt]) - nrm[tpt];
      std::cerr << std::setw(5) << tpt << "   |   "
                << std::setprecision(15)
                << std::fixed << std::setw(19) << nrm[tpt][0]
                << std::fixed << std::setw(19) << nrm[tpt][1]
                << std::fixed << std::setw(19) << nrm[tpt][2] << "   |   "
                << std::fixed << std::setw(9) << std::setprecision(5) << error.Norm() / nrm[tpt].Norm() * 100.0 << " %"
                << std::endl;
   };
   std::cerr << "-------------------------------------------------------------------------------------------\n";
};

#endif

};

#endif
