/*!
\file polyhedron.hh
\brief Properties of regular convex polyhedra (Platonic solids)
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_POLYHEDRON_HH
#define SPECTRUM_POLYHEDRON_HH

#include <string>
#include <array>

#include <common/definitions.hh>

namespace Spectrum {

//! Number of polyhedron types
#define N_POLYTYPES 8

/*!
\brief A class enumerating common polyhedra used as templates for mesh objects
\author Vladimir Florinski

The code supports all five regular polyhedra (tetrahedron, hexahedron, octahedron, dodecahedron, and icosahedron) and three types of prisms (with 3, 5, and 6 sides). The hexahedron (i.e., cube) is both a regular polyhedron and a 4-sided prism, so it has two names. Regular polyhedra are used for spherical tesselations and prisms are used to define block topologies.
*/
enum class PolyType {
   POLY_TETRAHEDRON,
   POLY_THREE_SIDED_PRISM,
   POLY_HEXAHEDRON,
   POLY_FOUR_SIDED_PRISM = POLY_HEXAHEDRON,
   POLY_OCTAHEDRON,
   POLY_FIVE_SIDED_PRISM,
   POLY_DODECAHEDRON,
   POLY_SIX_SIDED_PRISM,
   POLY_ICOSAHEDRON
};

//! Polyhedron names
inline constexpr std::array<std::string_view, N_POLYTYPES> PolyNames = {"tetrahedron"     , "three sided prism", "hexahedron"     , "octahedron" ,
                                                                        "five sided prism", "dodecahedron"     , "six sided prism", "icosahedron"};

/*!
\brief Schlafli symbol {p,q}
\author Vladimir Florinski
\author Lucius Schoenbaum
*/
template <int p, int q, bool prism = false>
struct Schlafli
{
//! Number of vertices
   SPECTRUM_DEVICE_FUNC static constexpr int NVerts(void) {return (prism ? 2 * p : 4 * p / (4 - (p - 2) * (q - 2)));};

//! Number of edges
   SPECTRUM_DEVICE_FUNC static constexpr int NEdges(void) {return (prism ? 3 * p : 2 * p * q / (4 - (p - 2) * (q - 2)));};

//! Number of faces
   SPECTRUM_DEVICE_FUNC static constexpr int NFaces(void) {return (prism ? 2 + p : 4 * q / (4 - (p - 2) * (q - 2)));};
};

//! Number of vertices/edges of each face (for prisms, this is for top/bottom faces, the sides all have p=4)
inline constexpr int poly_p[] = {3, 3, 4, 3, 5, 5, 6, 3};

//! Number of vertex/edge neighbors of each vertex
inline constexpr int poly_q[] = {3, 3, 3, 4, 3, 3, 3, 5};

//! Flag telling whether the polyhedron is a prism (the hexahedron is considered a prism even though it is also a regular polyhedron).
inline constexpr bool is_prism[] = {false, true, true, false, true, false, true, false};

/*!
\brief A class describing a convex polyhedron
\author Vladimir Florinski
\author Lucius Schoenbaum

This class provides a geometrical description of a convex polyhedron, with all of its edges equal in length, inscribed in a unit sphere. It requires no memory because all its members are static. It provides VEF counts, polar coordinates of the vertices on the unit sphere, and VV and FV connectivity tables. It does not have a public interface, so the class it intended as a base for mesh objects, such as blocks and tesselations. The unit-sphere coordinates are typically not useful for blocks, except to visualize the elements of the block.
*/
template <PolyType poly_type>
class Polyhedron
{
protected:

//! Number of vertices
   static constexpr int Nv = Schlafli<poly_p[static_cast<size_t>(poly_type)], poly_q[static_cast<size_t>(poly_type)],
                                    is_prism[static_cast<size_t>(poly_type)]>::NVerts(); 

//! Number of edges
   static constexpr int Ne = Schlafli<poly_p[static_cast<size_t>(poly_type)], poly_q[static_cast<size_t>(poly_type)],
                                    is_prism[static_cast<size_t>(poly_type)]>::NEdges();

//! Bumber of faces
   static constexpr int Nf = Schlafli<poly_p[static_cast<size_t>(poly_type)], poly_q[static_cast<size_t>(poly_type)],
                                    is_prism[static_cast<size_t>(poly_type)]>::NFaces();

//! Vertex latitudes
   static inline double vlat[Nv];

//! Vertex longitudes
   static inline double vlon[Nv];

//! Vertex-vertex connectivity array (ordered counter-clockwise)
   static inline int vert_vert[Nv][poly_q[static_cast<size_t>(poly_type)]];

//! Face-vertex connectivity array (ordered counter-clockwise)
   static inline int face_vert[Nf][poly_p[static_cast<size_t>(poly_type)]];

//! Coordinates of the vertices
   SPECTRUM_DEVICE_FUNC static constexpr void VertexCoords(void);

//! VV and FV connectivity
   SPECTRUM_DEVICE_FUNC static constexpr void SetConnectivity(void);

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC constexpr Polyhedron(void);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 03/17/2025
*/
template <PolyType poly_type>
SPECTRUM_DEVICE_FUNC inline constexpr Polyhedron<poly_type>::Polyhedron(void)
{
   VertexCoords();
   SetConnectivity();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Tetrahedron
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                      0
                    . | .
                  .   |   .
                .     |1    .
              .       |       .
            3 . . . . | . . . . 2
              .   2   |   0   .
                .     |3    .
                  .   |   .
                    . | . 
                      1 
*/

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_TETRAHEDRON>::VertexCoords(void)
{
   double poly_ang = asin(1.0 / 3.0);

   vlat[0] = 0.5 * M_PI;
   vlat[1] = -poly_ang;  vlat[2] = -poly_ang;        vlat[3] = -poly_ang;

   vlon[0] = 0.0 * M_PI;
   vlon[1] = 0.0 * M_PI; vlon[2] = 2.0 / 3.0 * M_PI; vlon[3] = 4.0 / 3.0 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_TETRAHEDRON>::SetConnectivity(void)
{
   vert_vert[0][0] = 1; vert_vert[0][1] = 2; vert_vert[0][2] = 3;
   vert_vert[1][0] = 0; vert_vert[1][1] = 3; vert_vert[1][2] = 2;
   vert_vert[2][0] = 0; vert_vert[2][1] = 1; vert_vert[2][2] = 3;
   vert_vert[3][0] = 0; vert_vert[3][1] = 2; vert_vert[3][2] = 1;

   face_vert[0][0] = 0; face_vert[0][1] = 1; face_vert[0][2] = 2;
   face_vert[1][0] = 0; face_vert[1][1] = 2; face_vert[1][2] = 3;
   face_vert[2][0] = 0; face_vert[2][1] = 3; face_vert[2][2] = 1;
   face_vert[3][0] = 1; face_vert[3][1] = 3; face_vert[3][2] = 2;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Three sided prism
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                      2
                    . . .
                  .   .   .
                .     0     .
              .       .       .
            0-------------------1
            |    3    .    2    |
            |         5         |
            |       .   .       |   front face = 1
            |     .       .     |
            |   .     4     .   |
            | .               . |
            3-------------------4
*/

/*!
\author Vladimir Florinski
\date 09/30/2025
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_THREE_SIDED_PRISM>::VertexCoords(void)
{
   double poly_ang = atan(0.5 * InscribedPolygonSide<poly_p[static_cast<size_t>(PolyType::POLY_THREE_SIDED_PRISM)]>());

   vlat[0] =  poly_ang;        vlat[1] =  poly_ang;        vlat[2] =  poly_ang;
   vlat[3] = -poly_ang;        vlat[4] = -poly_ang;        vlat[5] = -poly_ang;

   vlon[0] = 2.0 / 3.0 * M_PI; vlon[1] = 4.0 / 3.0 * M_PI; vlon[2] = 0.0 * M_PI;
   vlon[3] = 2.0 / 3.0 * M_PI; vlon[4] = 4.0 / 3.0 * M_PI; vlon[5] = 0.0 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_THREE_SIDED_PRISM>::SetConnectivity(void)
{
   vert_vert[0][0] = 1; vert_vert[0][1] = 2; vert_vert[0][2] = 3;
   vert_vert[1][0] = 2; vert_vert[1][1] = 0; vert_vert[1][2] = 4;
   vert_vert[2][0] = 0; vert_vert[2][1] = 1; vert_vert[2][2] = 5;
   vert_vert[3][0] = 4; vert_vert[3][1] = 0; vert_vert[3][2] = 5;
   vert_vert[4][0] = 5; vert_vert[4][1] = 1; vert_vert[4][2] = 3;
   vert_vert[5][0] = 3; vert_vert[5][1] = 2; vert_vert[5][2] = 4;

   face_vert[0][0] = 0; face_vert[0][1] = 1; face_vert[0][2] = 2; face_vert[0][3] =-1;
   face_vert[1][0] = 1; face_vert[1][1] = 0; face_vert[1][2] = 3; face_vert[1][3] = 4;
   face_vert[2][0] = 2; face_vert[2][1] = 1; face_vert[2][2] = 3; face_vert[2][3] = 5;
   face_vert[3][0] = 0; face_vert[3][1] = 2; face_vert[3][2] = 5; face_vert[3][3] = 3;
   face_vert[4][0] = 3; face_vert[4][1] = 5; face_vert[4][2] = 4; face_vert[4][3] =-1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Hexahedron
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                   3---------------2
                 . .             . |
               .   .   0       .   |
             .     .         .     |   front face = 1
            0---------------1      |   rear  face = 3
            |   4  .        |   2  |
            |      .        |      |
            |      7 . . . .|. . . 6
            |    .          |    .
            |  .       5    |  .
            |.              |.
            4---------------5
*/

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_HEXAHEDRON>::VertexCoords(void)
{
   double poly_ang = atan(0.5 * InscribedPolygonSide<poly_p[static_cast<size_t>(PolyType::POLY_HEXAHEDRON)]>());


   vlat[0] =  poly_ang;   vlat[1] =  poly_ang;   vlat[2] =  poly_ang;   vlat[3] =  poly_ang;
   vlat[4] = -poly_ang;   vlat[5] = -poly_ang;   vlat[6] = -poly_ang;   vlat[7] = -poly_ang;

   vlon[0] = 1.75 * M_PI; vlon[1] = 0.25 * M_PI; vlon[2] = 0.75 * M_PI; vlon[3] = 1.25 * M_PI;
   vlon[4] = 1.75 * M_PI; vlon[5] = 0.25 * M_PI; vlon[6] = 0.75 * M_PI; vlon[7] = 1.25 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_HEXAHEDRON>::SetConnectivity(void)
{
   vert_vert[0][0] = 1; vert_vert[0][1] = 3; vert_vert[0][2] = 4;
   vert_vert[1][0] = 2; vert_vert[1][1] = 0; vert_vert[1][2] = 5;
   vert_vert[2][0] = 3; vert_vert[2][1] = 1; vert_vert[2][2] = 6;
   vert_vert[3][0] = 0; vert_vert[3][1] = 2; vert_vert[3][2] = 7;
   vert_vert[4][0] = 5; vert_vert[4][1] = 0; vert_vert[4][2] = 7;
   vert_vert[5][0] = 6; vert_vert[5][1] = 1; vert_vert[5][2] = 4;
   vert_vert[6][0] = 7; vert_vert[6][1] = 2; vert_vert[6][2] = 5;
   vert_vert[7][0] = 4; vert_vert[7][1] = 3; vert_vert[7][2] = 6;

   face_vert[0][0] = 0; face_vert[0][1] = 1; face_vert[0][2] = 2; face_vert[0][3] = 3;
   face_vert[1][0] = 1; face_vert[1][1] = 0; face_vert[1][2] = 4; face_vert[1][3] = 5;
   face_vert[2][0] = 2; face_vert[2][1] = 1; face_vert[2][2] = 5; face_vert[2][3] = 6;
   face_vert[3][0] = 3; face_vert[3][1] = 2; face_vert[3][2] = 6; face_vert[3][3] = 7;
   face_vert[4][0] = 0; face_vert[4][1] = 3; face_vert[4][2] = 7; face_vert[4][3] = 4;
   face_vert[5][0] = 4; face_vert[5][1] = 7; face_vert[5][2] = 6; face_vert[5][3] = 5;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Octahedron
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                             0
                           . ...
                         . .  .    .
                       .  .    .       .
                     .   .     .    1      .
                   .    4 . . . . . . . . . . 3
                 .   .   .0     .          . .
               .  .      .       .      .  .
             . .          .      .   .   .        top    faces: 2, 3
            1---------------------2    .          bottom faces: 6, 7
               .           .     /  5.
                   .      4.    /  .
                       .     . / .
                           . ./.
                             5
*/

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_OCTAHEDRON>::VertexCoords(void)
{
   vlat[0] =  0.5 * M_PI;
   vlat[1] =  0.0 * M_PI; vlat[2] = 0.0 * M_PI; vlat[3] = 0.0 * M_PI; vlat[4] = 0.0 * M_PI;
   vlat[5] = -0.5 * M_PI;

   vlon[0] =  0.0 * M_PI;
   vlon[1] =  0.0 * M_PI; vlon[2] = 0.5 * M_PI; vlon[3] = 1.0 * M_PI; vlon[4] = 1.5 * M_PI;
   vlon[5] =  0.0 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_OCTAHEDRON>::SetConnectivity(void)
{
   vert_vert[0][0] = 1; vert_vert[0][1] = 2; vert_vert[0][2] = 3; vert_vert[0][3] = 4;
   vert_vert[1][0] = 0; vert_vert[1][1] = 4; vert_vert[1][2] = 5; vert_vert[1][3] = 2;
   vert_vert[2][0] = 0; vert_vert[2][1] = 1; vert_vert[2][2] = 5; vert_vert[2][3] = 3;
   vert_vert[3][0] = 0; vert_vert[3][1] = 2; vert_vert[3][2] = 5; vert_vert[3][3] = 4;
   vert_vert[4][0] = 0; vert_vert[4][1] = 3; vert_vert[4][2] = 5; vert_vert[4][3] = 1;
   vert_vert[5][0] = 1; vert_vert[5][1] = 4; vert_vert[5][2] = 3; vert_vert[5][3] = 2;

   face_vert[0][0] = 0; face_vert[0][1] = 1; face_vert[0][2] = 2;
   face_vert[1][0] = 0; face_vert[1][1] = 2; face_vert[1][2] = 3;
   face_vert[2][0] = 0; face_vert[2][1] = 3; face_vert[2][2] = 4;
   face_vert[3][0] = 0; face_vert[3][1] = 4; face_vert[3][2] = 1;
   face_vert[4][0] = 5; face_vert[4][1] = 2; face_vert[4][2] = 1;
   face_vert[5][0] = 5; face_vert[5][1] = 3; face_vert[5][2] = 2;
   face_vert[6][0] = 5; face_vert[6][1] = 4; face_vert[6][2] = 3;
   face_vert[7][0] = 5; face_vert[7][1] = 1; face_vert[7][2] = 4;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Five sided prism
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                    3
                 .  .  .
              .     .     .
           .        .        .
        4           0           2
        |.    4     .     3    .|
        | .         .         . |   front face = 1
        |  .        8        '  |
        |   0---------------1   |
        |   |  .          . |   |
        | 5.|               |.2 |
        9   |       6       |   7
         .  |               |  .
          . |               | .
           .|               |.
            5---------------6
*/

/*!
\author Vladimir Florinski
\date 10/07/2025
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_FIVE_SIDED_PRISM>::VertexCoords(void)
{
   double poly_ang = atan(0.5 * InscribedPolygonSide<poly_p[static_cast<size_t>(PolyType::POLY_FIVE_SIDED_PRISM)]>());

   vlat[0] =  poly_ang;        vlat[1] =  poly_ang;        vlat[2] =  poly_ang; vlat[3] =  poly_ang;        vlat[4] =  poly_ang;
   vlat[5] = -poly_ang;        vlat[6] = -poly_ang;        vlat[7] = -poly_ang; vlat[8] = -poly_ang;        vlat[9] = -poly_ang;

   vlon[0] = 1.0 / 5.0 * M_PI; vlon[1] = 3.0 / 5.0 * M_PI; vlon[2] = M_PI;      vlon[3] = 7.0 / 5.0 * M_PI; vlon[4] = 9.0 / 5.0 * M_PI;
   vlon[5] = 1.0 / 5.0 * M_PI; vlon[6] = 3.0 / 5.0 * M_PI; vlon[7] = M_PI;      vlon[8] = 7.0 / 5.0 * M_PI; vlon[9] = 9.0 / 5.0 * M_PI;
};

/*!
\author Vladimir Florinski
\date 10/07/2025
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_FIVE_SIDED_PRISM>::SetConnectivity(void)
{
   vert_vert[0][0] = 1; vert_vert[0][1] = 4; vert_vert[0][2] = 5;
   vert_vert[1][0] = 2; vert_vert[1][1] = 0; vert_vert[1][2] = 6;
   vert_vert[2][0] = 3; vert_vert[2][1] = 1; vert_vert[2][2] = 7;
   vert_vert[3][0] = 4; vert_vert[3][1] = 2; vert_vert[3][2] = 8;
   vert_vert[4][0] = 0; vert_vert[4][1] = 3; vert_vert[4][2] = 9;
   vert_vert[5][0] = 6; vert_vert[5][1] = 0; vert_vert[5][2] = 9;
   vert_vert[6][0] = 7; vert_vert[6][1] = 1; vert_vert[6][2] = 5;
   vert_vert[7][0] = 8; vert_vert[7][1] = 2; vert_vert[7][2] = 6;
   vert_vert[8][0] = 9; vert_vert[8][1] = 3; vert_vert[8][2] = 7;
   vert_vert[9][0] = 5; vert_vert[9][1] = 4; vert_vert[9][2] = 8;

   face_vert[0][0] = 0; face_vert[0][1] = 1; face_vert[0][2] = 2; face_vert[0][3] = 3; face_vert[0][4] = 4;
   face_vert[1][0] = 1; face_vert[1][1] = 0; face_vert[1][2] = 5; face_vert[1][3] = 6; face_vert[1][4] =-1;
   face_vert[2][0] = 2; face_vert[2][1] = 1; face_vert[2][2] = 6; face_vert[2][3] = 7; face_vert[2][4] =-1;
   face_vert[3][0] = 3; face_vert[3][1] = 2; face_vert[3][2] = 7; face_vert[3][3] = 8; face_vert[3][4] =-1;
   face_vert[4][0] = 4; face_vert[4][1] = 3; face_vert[4][2] = 8; face_vert[4][3] = 9; face_vert[4][4] =-1;
   face_vert[5][0] = 0; face_vert[5][1] = 4; face_vert[5][2] = 9; face_vert[5][3] = 5; face_vert[5][4] =-1;
   face_vert[6][0] = 5; face_vert[6][1] = 9; face_vert[6][2] = 8; face_vert[6][3] = 7; face_vert[6][4] = 6;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Dodecahedron
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
               14-----------13        13-----------12
                |             \       /             |
                |              \     /              |
                |               \   /               |
                |        5       \ /       4        |
                5                 3                11
                  .             .   .             .
                     .        .       .        .
                        .   .           .   .
                5---------4       0       2--------11       View from outside
               /           \             /           \
              /             \           /             \
             /               \         /               \
            6        1        0-------1         3      10
              .             ./         \.             .
                .        .  /           \  .        .
                  .   .    /             \    .   .
                    7     7       2       9     9
                            .           .
                              .       .
                                .   .
                                  8

                                 13
                               .     .
                             .         .
                           .             .
                   14    14      10      12    12
                  .   .    \             /    .   .
                .        .  \           /  .        .
              .             .\         /.             .
            5        6       15------19        9      11
             \               /         \               /
              \             /           \             /
               \           /     11      \           /
                6--------16              18--------10       View from inside
                        .   .           .   .
                     .        .       .        .
                  .             .   .             .
                6                17                10
                |        7       / \        8       |
                |               /   \               |
                |              /     \              |
                |             /       \             |
                7------------8         8------------9
*/

/*!
\author Vladimir Florinski
\date 03/17/2025
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_DODECAHEDRON>::VertexCoords(void)
{
   double edge_length = 4.0 / M_SQRT3 / (1.0 + M_SQRT5);
   double poly_ang1 = acos(edge_length * M_SQRT2 / sqrt(5.0 - M_SQRT5));
   double poly_ang2 = poly_ang1 - 2.0 * asin(0.5 * edge_length);

   vlat[ 0] =  poly_ang1; vlat[ 1] =  poly_ang1; vlat[ 2] =  poly_ang1; vlat[ 3] =  poly_ang1; vlat[ 4] =  poly_ang1;
   vlat[ 5] =  poly_ang2; vlat[ 6] = -poly_ang2; vlat[ 7] =  poly_ang2; vlat[ 8] = -poly_ang2; vlat[ 9] =  poly_ang2;
   vlat[10] = -poly_ang2; vlat[11] =  poly_ang2; vlat[12] = -poly_ang2; vlat[13] =  poly_ang2; vlat[14] = -poly_ang2;
   vlat[15] = -poly_ang1; vlat[16] = -poly_ang1; vlat[17] = -poly_ang1; vlat[18] = -poly_ang1; vlat[19] = -poly_ang1;

   vlon[ 0] = 0.4 * M_PI; vlon[ 1] = 0.8 * M_PI; vlon[ 2] = 1.2 * M_PI; vlon[ 3] = 1.6 * M_PI; vlon[ 4] = 0.0 * M_PI;
   vlon[ 5] = 0.0 * M_PI; vlon[ 6] = 0.2 * M_PI; vlon[ 7] = 0.4 * M_PI; vlon[ 8] = 0.6 * M_PI; vlon[ 9] = 0.8 * M_PI;
   vlon[10] = 1.0 * M_PI; vlon[11] = 1.2 * M_PI; vlon[12] = 1.4 * M_PI; vlon[13] = 1.6 * M_PI; vlon[14] = 1.8 * M_PI;
   vlon[15] = 1.8 * M_PI; vlon[16] = 0.2 * M_PI; vlon[17] = 0.6 * M_PI; vlon[18] = 1.0 * M_PI; vlon[19] = 1.4 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_DODECAHEDRON>::SetConnectivity(void)
{
   vert_vert[ 0][0] =  1; vert_vert[ 0][1] =  4; vert_vert[ 0][2] =  7;
   vert_vert[ 1][0] =  2; vert_vert[ 1][1] =  0; vert_vert[ 1][2] =  9;
   vert_vert[ 2][0] =  3; vert_vert[ 2][1] =  1; vert_vert[ 2][2] = 11;
   vert_vert[ 3][0] =  4; vert_vert[ 3][1] =  2; vert_vert[ 3][2] = 13;
   vert_vert[ 4][0] =  0; vert_vert[ 4][1] =  3; vert_vert[ 4][2] =  5;
   vert_vert[ 5][0] =  4; vert_vert[ 5][1] = 14; vert_vert[ 5][2] =  6;
   vert_vert[ 6][0] = 16; vert_vert[ 6][1] =  7; vert_vert[ 6][2] =  5;
   vert_vert[ 7][0] =  0; vert_vert[ 7][1] =  6; vert_vert[ 7][2] =  8;
   vert_vert[ 8][0] = 17; vert_vert[ 8][1] =  9; vert_vert[ 8][2] =  7;
   vert_vert[ 9][0] =  1; vert_vert[ 9][1] =  8; vert_vert[ 9][2] = 10;
   vert_vert[10][0] = 18; vert_vert[10][1] = 11; vert_vert[10][2] =  9;
   vert_vert[11][0] =  2; vert_vert[11][1] = 10; vert_vert[11][2] = 12;
   vert_vert[12][0] = 19; vert_vert[12][1] = 13; vert_vert[12][2] = 11;
   vert_vert[13][0] =  3; vert_vert[13][1] = 12; vert_vert[13][2] = 14;
   vert_vert[14][0] = 15; vert_vert[14][1] =  5; vert_vert[14][2] = 13;
   vert_vert[15][0] = 19; vert_vert[15][1] = 16; vert_vert[15][2] = 14;
   vert_vert[16][0] = 15; vert_vert[16][1] = 17; vert_vert[16][2] =  6;
   vert_vert[17][0] = 16; vert_vert[17][1] = 18; vert_vert[17][2] =  8;
   vert_vert[18][0] = 17; vert_vert[18][1] = 19; vert_vert[18][2] = 10;
   vert_vert[19][0] = 18; vert_vert[19][1] = 15; vert_vert[19][2] = 12;

   face_vert[ 0][0] =  0; face_vert[ 0][1] =  1; face_vert[ 0][2] =  2; face_vert[ 0][3] =  3; face_vert[ 0][4] =  4;
   face_vert[ 1][0] =  0; face_vert[ 1][1] =  4; face_vert[ 1][2] =  5; face_vert[ 1][3] =  6; face_vert[ 1][4] =  7;
   face_vert[ 2][0] =  1; face_vert[ 2][1] =  0; face_vert[ 2][2] =  7; face_vert[ 2][3] =  8; face_vert[ 2][4] =  9;
   face_vert[ 3][0] =  2; face_vert[ 3][1] =  1; face_vert[ 3][2] =  9; face_vert[ 3][3] = 10; face_vert[ 3][4] = 11;
   face_vert[ 4][0] =  3; face_vert[ 4][1] =  2; face_vert[ 4][2] = 11; face_vert[ 4][3] = 12; face_vert[ 4][4] = 13;
   face_vert[ 5][0] =  4; face_vert[ 5][1] =  3; face_vert[ 5][2] = 13; face_vert[ 5][3] = 14; face_vert[ 5][4] =  5;
   face_vert[ 6][0] = 15; face_vert[ 6][1] = 16; face_vert[ 6][2] =  6; face_vert[ 6][3] =  5; face_vert[ 6][4] = 14;
   face_vert[ 7][0] = 16; face_vert[ 7][1] = 17; face_vert[ 7][2] =  8; face_vert[ 7][3] =  7; face_vert[ 7][4] =  6;
   face_vert[ 8][0] = 17; face_vert[ 8][1] = 18; face_vert[ 8][2] = 10; face_vert[ 8][3] =  9; face_vert[ 8][4] =  8;
   face_vert[ 9][0] = 18; face_vert[ 9][1] = 19; face_vert[ 9][2] = 12; face_vert[ 9][3] = 11; face_vert[ 9][4] = 10;
   face_vert[10][0] = 19; face_vert[10][1] = 15; face_vert[10][2] = 14; face_vert[10][3] = 13; face_vert[10][4] = 12;
   face_vert[11][0] = 19; face_vert[11][1] = 18; face_vert[11][2] = 17; face_vert[11][3] = 16; face_vert[11][4] = 15;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Six sided prism
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                4---------------3
              . .               ..
            .   .               . .
          .     .               .  .
        5       .      0        .   2
        |.      .               . . |
        | . 5   .               . 3 |   front face = 1
        |  .   10...............9   |   rear  face = 4
        |   0---------------1    .  |
        |   |               |     . |
        | 6 |               |   2  .|
       11   |          7    |       8
         .  |               |     .
          . |               |   .
           .|               | .
            6---------------7
*/

/*!
\author Vladimir Florinski
\date 09/30/2025
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_SIX_SIDED_PRISM>::VertexCoords(void)
{
   double poly_ang = atan(0.5 * InscribedPolygonSide<poly_p[static_cast<size_t>(PolyType::POLY_SIX_SIDED_PRISM)]>());

   vlat[ 0] =  poly_ang;        vlat[ 1] =  poly_ang;  vlat[ 2] =  poly_ang;
   vlat[ 3] =  poly_ang;        vlat[ 4] =  poly_ang;  vlat[ 5] =  poly_ang;
   vlat[ 6] = -poly_ang;        vlat[ 7] = -poly_ang;  vlat[ 8] = -poly_ang;
   vlat[ 9] = -poly_ang;        vlat[10] = -poly_ang;  vlat[11] = -poly_ang;

   vlon[ 0] = 1.0 / 6.0 * M_PI; vlon[ 1] = 0.5 * M_PI; vlon[ 2] =  5.0 / 6.0 * M_PI;
   vlon[ 3] = 7.0 / 6.0 * M_PI; vlon[ 4] = 1.5 * M_PI; vlon[ 5] = 11.0 / 6.0 * M_PI;
   vlon[ 6] = 1.0 / 6.0 * M_PI; vlon[ 7] = 0.5 * M_PI; vlon[ 8] =  5.0 / 6.0 * M_PI;
   vlon[ 9] = 7.0 / 6.0 * M_PI; vlon[10] = 1.5 * M_PI; vlon[11] = 11.0 / 6.0 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_SIX_SIDED_PRISM>::SetConnectivity(void)
{
   vert_vert[ 0][0] =  1; vert_vert[ 0][1] =  5; vert_vert[ 0][2] =  6;
   vert_vert[ 1][0] =  2; vert_vert[ 1][1] =  0; vert_vert[ 1][2] =  7;
   vert_vert[ 2][0] =  3; vert_vert[ 2][1] =  1; vert_vert[ 2][2] =  8;
   vert_vert[ 3][0] =  4; vert_vert[ 3][1] =  2; vert_vert[ 3][2] =  9;
   vert_vert[ 4][0] =  5; vert_vert[ 4][1] =  3; vert_vert[ 4][2] = 10;
   vert_vert[ 5][0] =  0; vert_vert[ 5][1] =  4; vert_vert[ 5][2] = 11;
   vert_vert[ 6][0] =  7; vert_vert[ 6][1] =  0; vert_vert[ 6][2] = 11;
   vert_vert[ 7][0] =  8; vert_vert[ 7][1] =  1; vert_vert[ 7][2] =  6;
   vert_vert[ 8][0] =  9; vert_vert[ 8][1] =  2; vert_vert[ 8][2] =  7;
   vert_vert[ 9][0] = 10; vert_vert[ 9][1] =  3; vert_vert[ 9][2] =  8;
   vert_vert[10][0] = 11; vert_vert[10][1] =  4; vert_vert[10][2] =  9;
   vert_vert[11][0] =  6; vert_vert[11][1] =  5; vert_vert[11][2] = 10;

   face_vert[ 0][0] =  0; face_vert[ 0][1] =  1; face_vert[ 0][2] =  2; face_vert[ 0][3] =  3; face_vert[ 0][4] = 4; face_vert[ 0][5] = 5;
   face_vert[ 1][0] =  1; face_vert[ 1][1] =  0; face_vert[ 1][2] =  6; face_vert[ 1][3] =  7; face_vert[ 1][4] =-1; face_vert[ 1][5] =-1;
   face_vert[ 2][0] =  2; face_vert[ 2][1] =  1; face_vert[ 2][2] =  7; face_vert[ 2][3] =  8; face_vert[ 2][4] =-1; face_vert[ 2][5] =-1;
   face_vert[ 3][0] =  3; face_vert[ 3][1] =  2; face_vert[ 3][2] =  8; face_vert[ 3][3] =  9; face_vert[ 3][4] =-1; face_vert[ 3][5] =-1;
   face_vert[ 4][0] =  4; face_vert[ 4][1] =  3; face_vert[ 4][2] =  9; face_vert[ 4][3] = 10; face_vert[ 4][4] =-1; face_vert[ 4][5] =-1;
   face_vert[ 5][0] =  5; face_vert[ 5][1] =  4; face_vert[ 5][2] = 10; face_vert[ 5][3] = 11; face_vert[ 5][4] =-1; face_vert[ 5][5] =-1;
   face_vert[ 6][0] =  0; face_vert[ 6][1] =  5; face_vert[ 6][2] = 11; face_vert[ 6][3] =  6; face_vert[ 6][4] =-1; face_vert[ 6][5] =-1;
   face_vert[ 7][0] =  6; face_vert[ 7][1] = 11; face_vert[ 7][2] = 10; face_vert[ 7][3] =  9; face_vert[ 7][4] = 8; face_vert[ 7][5] = 7;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Icosahedron
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
                                      10
                                     .   .
                                   .       .
                                 .           .
                               .      13       .
                             .                   .
                        .  5-----------------------4  .
                    .      . .                   . .      .
                .              .       3       .              .
            6             .      .           .      .             9
              .      5             .       .            11      .
                .            4       .   .       2            .
                  .      .          .  0  .          .      .         View from outside
                    .           .      |      .           .
                      . .   .          |          .   . .
                        1         0    |    1         3
                        |  .           |           .  |
                        |      .       |       .      |
                        |          .   |   .          |
                        |   7       .  2  .       9   |
                        |       .             .       |
                        |   .                     .   |
                        7                             8

                        5                             4
                        |   .                     .   |
                        |       .             .       |
                        |           . 10  .           |
                        |   14     .   |   .    12    |
                        |      .       |       .      |
                        |  .           |           .  |
                        6         15   |   19         9
                      . .   .          |          .   . .
                    .           .      |      .           .
                  .      .          . 11  .          .      .         View from inside
                .           16       .   .      18            .
              .      6             .       .             10     .
            1             .      .           .      .             3
                .              .      17       .              .
                    .      . .                   . .      .
                        .  7-----------------------8  .
                             .                   .
                               .       8       .
                                 .          .
                                   .       .
                                     .   .
                                       2
*/

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_ICOSAHEDRON>::VertexCoords(void)
{
   double poly_ang = atan(0.5);

   vlat[ 0] =  0.5 * M_PI;
   vlat[ 1] =  poly_ang;   vlat[ 2] =  poly_ang;  vlat[ 3] =  poly_ang;  vlat[ 4] =  poly_ang;  vlat[ 5] =  poly_ang;
   vlat[ 6] = -poly_ang;   vlat[ 7] = -poly_ang;  vlat[ 8] = -poly_ang;  vlat[ 9] = -poly_ang;  vlat[10] = -poly_ang;
   vlat[11] = -0.5 * M_PI;

   vlon[ 0] =  0.0 * M_PI;
   vlon[ 1] =  0.2 * M_PI; vlon[ 2] = 0.6 * M_PI; vlon[ 3] = 1.0 * M_PI; vlon[ 4] = 1.4 * M_PI; vlon[ 5] = 1.8 * M_PI;
   vlon[ 6] =  0.0 * M_PI; vlon[ 7] = 0.4 * M_PI; vlon[ 8] = 0.8 * M_PI; vlon[ 9] = 1.2 * M_PI; vlon[10] = 1.6 * M_PI;
   vlon[11] =  0.0 * M_PI;
};

/*!
\author Vladimir Florinski
\date 01/05/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline constexpr void Polyhedron<PolyType::POLY_ICOSAHEDRON>::SetConnectivity(void)
{
   vert_vert[ 0][0] =  1; vert_vert[ 0][1] =  2; vert_vert[ 0][2] =  3; vert_vert[ 0][3] =  4; vert_vert[ 0][4] =  5;
   vert_vert[ 1][0] =  0; vert_vert[ 1][1] =  5; vert_vert[ 1][2] =  6; vert_vert[ 1][3] =  7; vert_vert[ 1][4] =  2;
   vert_vert[ 2][0] =  0; vert_vert[ 2][1] =  1; vert_vert[ 2][2] =  7; vert_vert[ 2][3] =  8; vert_vert[ 2][4] =  3;
   vert_vert[ 3][0] =  0; vert_vert[ 3][1] =  2; vert_vert[ 3][2] =  8; vert_vert[ 3][3] =  9; vert_vert[ 3][4] =  4;
   vert_vert[ 4][0] =  0; vert_vert[ 4][1] =  3; vert_vert[ 4][2] =  9; vert_vert[ 4][3] = 10; vert_vert[ 4][4] =  5;
   vert_vert[ 5][0] =  0; vert_vert[ 5][1] =  4; vert_vert[ 5][2] = 10; vert_vert[ 5][3] =  6; vert_vert[ 5][4] =  1;
   vert_vert[ 6][0] = 11; vert_vert[ 6][1] =  7; vert_vert[ 6][2] =  1; vert_vert[ 6][3] =  5; vert_vert[ 6][4] = 10;
   vert_vert[ 7][0] = 11; vert_vert[ 7][1] =  8; vert_vert[ 7][2] =  2; vert_vert[ 7][3] =  1; vert_vert[ 7][4] =  6;
   vert_vert[ 8][0] = 11; vert_vert[ 8][1] =  9; vert_vert[ 8][2] =  3; vert_vert[ 8][3] =  2; vert_vert[ 8][4] =  7;
   vert_vert[ 9][0] = 11; vert_vert[ 9][1] = 10; vert_vert[ 9][2] =  4; vert_vert[ 9][3] =  3; vert_vert[ 9][4] =  8;
   vert_vert[10][0] = 11; vert_vert[10][1] =  6; vert_vert[10][2] =  5; vert_vert[10][3] =  4; vert_vert[10][4] =  9;
   vert_vert[11][0] = 10; vert_vert[11][1] =  9; vert_vert[11][2] =  8; vert_vert[11][3] =  7; vert_vert[11][4] =  6;

   face_vert[ 0][0] =  0; face_vert[ 0][1] =  1; face_vert[ 0][2] =  2;
   face_vert[ 1][0] =  0; face_vert[ 1][1] =  2; face_vert[ 1][2] =  3;
   face_vert[ 2][0] =  0; face_vert[ 2][1] =  3; face_vert[ 2][2] =  4;
   face_vert[ 3][0] =  0; face_vert[ 3][1] =  4; face_vert[ 3][2] =  5;
   face_vert[ 4][0] =  0; face_vert[ 4][1] =  5; face_vert[ 4][2] =  1;
   face_vert[ 5][0] =  6; face_vert[ 5][1] =  1; face_vert[ 5][2] =  5;
   face_vert[ 6][0] =  1; face_vert[ 6][1] =  6; face_vert[ 6][2] =  7;
   face_vert[ 7][0] =  7; face_vert[ 7][1] =  2; face_vert[ 7][2] =  1;
   face_vert[ 8][0] =  2; face_vert[ 8][1] =  7; face_vert[ 8][2] =  8;
   face_vert[ 9][0] =  8; face_vert[ 9][1] =  3; face_vert[ 9][2] =  2;
   face_vert[10][0] =  3; face_vert[10][1] =  8; face_vert[10][2] =  9;
   face_vert[11][0] =  9; face_vert[11][1] =  4; face_vert[11][2] =  3;
   face_vert[12][0] =  4; face_vert[12][1] =  9; face_vert[12][2] = 10;
   face_vert[13][0] = 10; face_vert[13][1] =  5; face_vert[13][2] =  4;
   face_vert[14][0] =  5; face_vert[14][1] = 10; face_vert[14][2] =  6;
   face_vert[15][0] = 11; face_vert[15][1] =  6; face_vert[15][2] = 10;
   face_vert[16][0] = 11; face_vert[16][1] =  7; face_vert[16][2] =  6;
   face_vert[17][0] = 11; face_vert[17][1] =  8; face_vert[17][2] =  7;
   face_vert[18][0] = 11; face_vert[18][1] =  9; face_vert[18][2] =  8;
   face_vert[19][0] = 11; face_vert[19][1] = 10; face_vert[19][2] =  9;
};

};

#endif
