/*!
\file quadrature.hh
\brief Quadrature point locations for line and area elements
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_QUADRATURE
#define SPECTRUM_QUADRATURE

#ifdef GEO_DEBUG
#include <string>
#include <fstream>
#include <iomanip>
#endif
#include "common/definitions.hh"

namespace Spectrum {

//! Largest number of quadrature points of Gauss-Legendre type
#define MAX_QPOINTS_GLG 3

//! The largest common polynomial order for which exact quadrature rules are available
#define MAX_QORDER 4

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Gauss-Legendre quadrature points on [-1/2,1/2]
// Ref: http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Gauss-Legendre quadrature displacement, 1 point (exact to \f$x^1\f$)
const double qdispl_glg_1[1] = {0.0};

//! Gauss-Legendre quadrature weight, 1 point
const double weight_glg_1[1] = {1.0};

//! Gauss-Legendre quadrature displacements, 2 points (exact to \f$x^3\f$)
const double qdispl_glg_2[2] = {-sqrtthr / 6.0, sqrtthr / 6.0};

//! Gauss-Legendre quadrature weights, 2 points
const double weight_glg_2[2] = {1.0 / 2.0, 1.0 / 2.0};

//! Gauss-Legendre quadrature displacements, 3 points (exact to \f$x^5\f$)
const double qdispl_glg_3[3] = {-sqrtthr * sqrtfiv / 10.0, 0.0, sqrtthr * sqrtfiv / 10.0};

//! Gauss-Legendre quadrature weights, 3 points
const double weight_glg_3[3] = {5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};

//! Gauss-Legendre quadrature displacements
const double* const qdispl_gauss_legendre[MAX_QPOINTS_GLG] = {qdispl_glg_1, qdispl_glg_2, qdispl_glg_3};
                                                             
//! Gauss-Legendre quadrature weights
const double* const weight_gauss_legendre[MAX_QPOINTS_GLG] = {weight_glg_1, weight_glg_2, weight_glg_3};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Quadrature points on a unit interval
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 1st order polynomial exactly
//                 ------------0------------

//! Reference coordinates for 1st order rule on a line
const double qpoints_lin_1order_coordx[1] = {0.5 + qdispl_glg_1[0]};

//! Weights for 1st order rule on a line
const double qpoints_lin_1order_weight[1] = {weight_glg_1[0]};

// Integrates 3rd order polynomial exactly
//                 ------0-----------1------

//! Reference coordinates for 2nd/3rd order rule on a line
const double qpoints_lin_3order_coordx[2] = {0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1]};

//! Weights for 2nd/3rd order rule on a line
const double qpoints_lin_3order_weight[2] = {weight_glg_2[0], weight_glg_2[1]};

// Integrates 5th order polynomial exactly
//                 ---0--------1--------2---

//! Reference coordinates for 4th/5th order rule on a line
const double qpoints_lin_5order_coordx[3] = {0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[2]};

//! Weights for 4th/5th order rule on a line
const double qpoints_lin_5order_weight[3] = {weight_glg_3[0], weight_glg_3[1], weight_glg_3[2]};

//----------------------------------------------------------------------------------------------------------------------------------------------------

const double* qpoints_lin_1order[2] = {qpoints_lin_1order_coordx, qpoints_lin_1order_weight};
const double* qpoints_lin_2order[2] = {qpoints_lin_3order_coordx, qpoints_lin_3order_weight};
const double* qpoints_lin_3order[2] = {qpoints_lin_3order_coordx, qpoints_lin_3order_weight};
const double* qpoints_lin_4order[2] = {qpoints_lin_5order_coordx, qpoints_lin_5order_weight};

//! Number of line quadrature points
const int n_qpoints_lin[MAX_QORDER] = {1, 2, 2, 3};

//! All properties of line quadrature points
const double* const* qpoints_lin[MAX_QORDER] = {qpoints_lin_1order, qpoints_lin_2order, qpoints_lin_3order, qpoints_lin_4order};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Quadrature points on an equilateral unit triangle
// Ref: Dunavant, D. A., High degree efficient symmetrical Gaussian quadrature rules for the triangle, International Journal for Numerical Methods in Engineering, v. 21, p. 1129 (1985).
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 1st order bi-variate polynomial exactly.
//                             .
//                            / \
//                           /   \
//                          /     \
//                         /       \
//                        /         \
//                       /           \
//                      /             \
//                     /       0       \
//                    /                 \
//                   /                   \
//                  /                     \
//                 -------------------------

//! First reference coordinate for 1st order rule on a triangle
const double qpoints_tri_1order_coordx[] = {1.0 / 2.0};

//! Second reference coordinate for 1st order rule on a triangle
const double qpoints_tri_1order_coordy[] = {1.0 / (2.0 * sqrtthr)};

//! Weights for 1st order rule on a triangle
const double qpoints_tri_1order_weight[] = {1.0};

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 2nd order bi-variate polynomial exactly.
//                             .
//                            / \
//                           /   \
//                          /     \
//                         /   2   \
//                        /         \
//                       /           \
//                      /             \
//                     /               \
//                    /                 \
//                   /    0         1    \
//                  /                     \
//                 -------------------------

//! First reference coordinate for 2nd order rule on a triangle
const double qpoints_tri_2order_coordx[] = {1.0 / 4.0, 3.0 / 4.0, 2.0 / 4.0};

//! Second reference coordinate for 2nd order rule on a triangle
const double qpoints_tri_2order_coordy[] = {1.0 / (4.0 * sqrtthr), 1.0 / (4.0 * sqrtthr), 1.0 / (1.0 * sqrtthr)};

//! Weights for 2nd order rule on a triangle
const double qpoints_tri_2order_weight[] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 3rd order bi-variate polynomial exactly. This rule has a negative weight and should be avoided.
//                             .
//                            / \
//                           /   \
//                          /     \
//                         /       \
//                        /    3    \
//                       /           \
//                      /             \
//                     /       0       \
//                    /                 \
//                   /     1       2     \
//                  /                     \
//                 -------------------------

//! First reference coordinate for 3rd order rule on a triangle
const double qpoints_tri_3order_coordx[] = { 5.0 / 10.0, 3.0 / 10.0, 7.0 / 10.0, 5.0 / 10.0};

//! Second reference coordinate for 2nd order rule on a triangle
const double qpoints_tri_3order_coordy[] = { 5.0 / (10.0 * sqrtthr), 3.0 / (10.0 * sqrtthr), 3.0 / (10.0 * sqrtthr), 9.0 / (10.0 * sqrtthr)};

//! Weights for 2nd order rule on a triangle
const double qpoints_tri_3order_weight[] = {-27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0};

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 4th order bi-variate polynomial exactly.
//                             .
//                            / \
//                           /   \
//                          /  2  \
//                         /       \
//                        /         \
//                       /           \
//                      /   5     4   \
//                     /               \
//                    /                 \
//                   /         3         \
//                  /  0               1  \
//                 -------------------------

//! First reference coordinate for 4th order rule on a triangle
const double qpoints_tri_4order_coordx[] = {(8.0 - sqrtten) / 12.0 - sqrtten * sqrt(95.0 - 22.0 * sqrtten) / 60.0,
                                            (4.0 + sqrtten) / 12.0 + sqrtten * sqrt(95.0 - 22.0 * sqrtten) / 60.0,
                                            1.0 / 2.0,
                                            1.0 / 2.0,
                                            (8.0 - sqrtten) / 12.0 + sqrtten * sqrt(95.0 - 22.0 * sqrtten) / 60.0,
                                            (4.0 + sqrtten) / 12.0 - sqrtten * sqrt(95.0 - 22.0 * sqrtten) / 60.0};

//! Second reference coordinate for 3rd/4th order rule on a triangle
const double qpoints_tri_4order_coordy[] = {(8.0 - sqrtten) / (12.0 * sqrtthr) - sqrtten * sqrt(95.0 - 22.0 * sqrtten) / (60.0 * sqrtthr),
                                            (8.0 - sqrtten) / (12.0 * sqrtthr) - sqrtten * sqrt(95.0 - 22.0 * sqrtten) / (60.0 * sqrtthr),
                                            (1.0 + sqrtten) / ( 6.0 * sqrtthr) + sqrtten * sqrt(95.0 - 22.0 * sqrtten) / (30.0 * sqrtthr),
                                            (1.0 + sqrtten) / ( 6.0 * sqrtthr) - sqrtten * sqrt(95.0 - 22.0 * sqrtten) / (30.0 * sqrtthr),
                                            (8.0 - sqrtten) / (12.0 * sqrtthr) + sqrtten * sqrt(95.0 - 22.0 * sqrtten) / (60.0 * sqrtthr),
                                            (8.0 - sqrtten) / (12.0 * sqrtthr) + sqrtten * sqrt(95.0 - 22.0 * sqrtten) / (60.0 * sqrtthr)};

//! Weights for 3rd/4th order rule on a triangle
const double qpoints_tri_4order_weight[] = {1.0 / 6.0 - (45.0 - sqrtten) * sqrt(95.0 - 22.0 * sqrtten) / 3720.0,
                                            1.0 / 6.0 - (45.0 - sqrtten) * sqrt(95.0 - 22.0 * sqrtten) / 3720.0,
                                            1.0 / 6.0 - (45.0 - sqrtten) * sqrt(95.0 - 22.0 * sqrtten) / 3720.0,
                                            1.0 / 6.0 + (45.0 - sqrtten) * sqrt(95.0 - 22.0 * sqrtten) / 3720.0,
                                            1.0 / 6.0 + (45.0 - sqrtten) * sqrt(95.0 - 22.0 * sqrtten) / 3720.0,
                                            1.0 / 6.0 + (45.0 - sqrtten) * sqrt(95.0 - 22.0 * sqrtten) / 3720.0};

//----------------------------------------------------------------------------------------------------------------------------------------------------

const double* qpoints_tri_1order[3] = {qpoints_tri_1order_coordx, qpoints_tri_1order_coordy, qpoints_tri_1order_weight};
const double* qpoints_tri_2order[3] = {qpoints_tri_2order_coordx, qpoints_tri_2order_coordy, qpoints_tri_2order_weight};
const double* qpoints_tri_3order[3] = {qpoints_tri_4order_coordx, qpoints_tri_4order_coordy, qpoints_tri_4order_weight};
const double* qpoints_tri_4order[3] = {qpoints_tri_4order_coordx, qpoints_tri_4order_coordy, qpoints_tri_4order_weight};

//! Number of triangle quadrature points
const int n_qpoints_tri[MAX_QORDER] = {1, 3, 6, 6};

//! All properties of triangle quadrature points
const double* const* qpoints_tri[MAX_QORDER] = {qpoints_tri_1order, qpoints_tri_2order, qpoints_tri_3order, qpoints_tri_4order};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Quadrature points on a unit square
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 1st order bi-variate polynomial exactly.
//                 -------------------------
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |           0           |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 -------------------------

//! First reference coordinate for 1st order rule on a rectangle
const double qpoints_rec_1order_coordx[] = {0.5 + qdispl_glg_1[0]};

//! Second reference coordinate for 1st order rule on a rectangle
const double qpoints_rec_1order_coordy[] = {0.5 + qdispl_glg_1[0]};

//! Weights for 1st order rule on a rectangle
const double qpoints_rec_1order_weight[] = {weight_glg_1[0]};

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 2nd order bi-variate polynomial exactly. 
//                 -------------------------
//                 |                       |
//                 |                       |
//                 |           2           |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |     0           1     |
//                 |                       |
//                 |                       |
//                 -------------------------

//! First reference coordinate for 2nd order rule on a rectangle
const double qpoints_rec_2order_coordx[] = {0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1], 0.5 + qdispl_glg_1[0]};

//! Second reference coordinate for 2nd order rule on a rectangle
const double qpoints_rec_2order_coordy[] = {0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1]};

//! Weights for 2nd order rule on a rectangle
const double qpoints_rec_2order_weight[] = {weight_glg_2[0] * weight_glg_2[0], weight_glg_2[1] * weight_glg_2[0], weight_glg_1[0] * weight_glg_2[1]};
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 3rd order bi-variate polynomial exactly.
//                 -------------------------
//                 |                       |
//                 |                       |
//                 |     2           3     |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |     0           1     |
//                 |                       |
//                 |                       |
//                 -------------------------

//! First reference coordinate for 3rd order rule on a rectangle
const double qpoints_rec_3order_coordx[] = {0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1],
                                            0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1]};

//! Second reference coordinate for 3rd order rule on a rectangle
const double qpoints_rec_3order_coordy[] = {0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[0],
                                            0.5 + qdispl_glg_2[1], 0.5 + qdispl_glg_2[1]};

//! Weights for 3rd order rule on a rectangle
const double qpoints_rec_3order_weight[] = {weight_glg_2[0] * weight_glg_2[0], weight_glg_2[1] * weight_glg_2[0],
                                            weight_glg_2[0] * weight_glg_2[1], weight_glg_2[1] * weight_glg_2[1]};

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 4th order bi-variate polynomial exactly. 
//                 -------------------------
//                 |                       |
//                 |     4           5     |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |     2           3     |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |     0           1     |
//                 |                       |
//                 -------------------------

//! First reference coordinate for 4th order rule on a rectangle
const double qpoints_rec_4order_coordx[] = {0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1],
                                            0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1],
                                            0.5 + qdispl_glg_2[0], 0.5 + qdispl_glg_2[1]};

//! Second reference coordinate for 4th order rule on a rectangle
const double qpoints_rec_4order_coordy[] = {0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[0],
                                            0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[1],
                                            0.5 + qdispl_glg_3[2], 0.5 + qdispl_glg_3[2]};

//! Weights for 4th order rule on a rectangle
const double qpoints_rec_4order_weight[] = {weight_glg_2[0] * weight_glg_3[0], weight_glg_2[1] * weight_glg_3[0],
                                            weight_glg_2[0] * weight_glg_3[1], weight_glg_2[1] * weight_glg_3[1],
                                            weight_glg_2[0] * weight_glg_3[2], weight_glg_2[1] * weight_glg_3[2]};

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Integrates 5th order bi-variate polynomial exactly.
//                 -------------------------
//                 |                       |
//                 |  6        7        8  |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |  3        4        5  |
//                 |                       |
//                 |                       |
//                 |                       |
//                 |  0        1        2  |
//                 |                       |
//                 -------------------------

//! First reference coordinate for 5th order rule on a rectangle
const double qpoints_rec_5order_coordx[] = {0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[2],
                                            0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[2],
                                            0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[2]};

//! Second reference coordinate for 5th order rule on a rectangle
const double qpoints_rec_5order_coordy[] = {0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[0], 0.5 + qdispl_glg_3[0],
                                            0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[1], 0.5 + qdispl_glg_3[1],
                                            0.5 + qdispl_glg_3[2], 0.5 + qdispl_glg_3[2], 0.5 + qdispl_glg_3[2]};

//! Weights for 5th order rule on a rectangle
const double qpoints_rec_5order_weight[] = {weight_glg_3[0] * weight_glg_3[0], weight_glg_3[1] * weight_glg_3[0], weight_glg_3[2] * weight_glg_3[0],
                                            weight_glg_3[0] * weight_glg_3[1], weight_glg_3[1] * weight_glg_3[1], weight_glg_3[2] * weight_glg_3[1],
                                            weight_glg_3[0] * weight_glg_3[2], weight_glg_3[1] * weight_glg_3[2], weight_glg_3[2] * weight_glg_3[2]};

//----------------------------------------------------------------------------------------------------------------------------------------------------

const double* qpoints_rec_1order[3] = {qpoints_rec_1order_coordx, qpoints_rec_1order_coordy, qpoints_rec_1order_weight};
const double* qpoints_rec_2order[3] = {qpoints_rec_3order_coordx, qpoints_rec_3order_coordy, qpoints_rec_3order_weight};
const double* qpoints_rec_3order[3] = {qpoints_rec_3order_coordx, qpoints_rec_3order_coordy, qpoints_rec_3order_weight};
const double* qpoints_rec_4order[3] = {qpoints_rec_5order_coordx, qpoints_rec_5order_coordy, qpoints_rec_5order_weight};

//! Number of rectangle quadrature points
const int n_qpoints_rec[MAX_QORDER] = {1, 4, 4, 9};

//! All properties of rectangle quadrature points
const double* const* qpoints_rec[MAX_QORDER] = {qpoints_rec_1order, qpoints_rec_2order, qpoints_rec_3order, qpoints_rec_4order};

//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef GEO_DEBUG

//! Reference coordinates of the vertices
const double vert_coord[2][4][2] = {{{0.0, 0.0}, {1.0, 0.0}, {1.0 / 2.0, sqrtthr / 2.0}, {0.0, 0.0}},
                                    {{0.0, 0.0}, {1.0, 0.0}, {1.0      , 1.0          }, {0.0, 1.0}}};

//! File name for quadrature point visualization
const std::string qp_fname = "quadpoints.dat";

/*!
\brief Visualize area quadrature points on a reference element
\author Vladimir Florinski
\date 04/03/2024
\param[in] n_verts Number of vertices
\param[in] order   Order of accuracy
*/
inline void QuadVisualize(int n_verts, int order)
{
   const int* n_qpoints = (n_verts == 3 ? n_qpoints_tri : n_qpoints_rec);
   const double* const** qpoints = (n_verts == 3 ? qpoints_tri : qpoints_rec);

   std::ofstream qpfile;
   qpfile.open(qp_fname.c_str(), std::ofstream::out);
   qpfile << std::setprecision(8);

// The closed shape of the reference polygon
   for(auto iv = 0; iv <= n_verts; iv++) {
      qpfile << std::setw(16) << vert_coord[n_verts - 3][iv % n_verts][0]
             << std::setw(16) << vert_coord[n_verts - 3][iv % n_verts][1]
             << std::endl;
   };
   qpfile << "&\n";

// Quadrature points
   for(auto qpt = 0; qpt < n_qpoints[order - 1]; qpt++) {
      qpfile << std::setw(12) << qpoints[order - 1][0][qpt]
             << std::setw(12) << qpoints[order - 1][1][qpt]
             << std::endl;
   };
   qpfile.close();
};

#endif

};

#endif
