/*!
\file polynomial.hh
\brief Declares data and common operation on polynomials
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_POLYNOMIAL_HH
#define SPECTRUM_POLYNOMIAL_HH

#include <cstdint>
#include "common/vectors.hh"

namespace Spectrum {

//! The highest monomial order for which bitwise encoding is available
#define MONO_ORDER_HIGH 5

#define DOF_1D_0  1
#define DOF_1D_1  2
#define DOF_1D_2  3
#define DOF_1D_3  4
#define DOF_1D_4  5
#define DOF_1D_5  6

#define DOF_2D_0  1
#define DOF_2D_1  3
#define DOF_2D_2  6
#define DOF_2D_3 10
#define DOF_2D_4 15
#define DOF_2D_5 21

#define DOF_3D_0  1
#define DOF_3D_1  4
#define DOF_3D_2 10
#define DOF_3D_3 20
#define DOF_3D_4 35
#define DOF_3D_5 56

//! Number of degrees of freedom in 1D, 2D, and 3D
const int deg_freedom[3][MONO_ORDER_HIGH + 1]
= {{DOF_1D_0, DOF_1D_1, DOF_1D_2, DOF_1D_3, DOF_1D_4, DOF_1D_5},
   {DOF_2D_0, DOF_2D_1, DOF_2D_2, DOF_2D_3, DOF_2D_4, DOF_2D_5},
   {DOF_3D_0, DOF_3D_1, DOF_3D_2, DOF_3D_3, DOF_3D_4, DOF_3D_5}};

/*!
\brief Returns the number of degrees of freeedom in one dimension
\author Vladimir Florinski
\date 07/17/2019
\param[in] order The order of reconstruction polynomial
\return The number of DOF
*/
SPECTRUM_DEVICE_FUNC inline int DOF_1D(int order)
{
   return order + 1;
};

/*!
\brief Returns the number of degrees of freeedom in two dimensions
\author Vladimir Florinski
\date 07/17/2019
\param[in] order The order of reconstruction polynomial
\return The number of DOF
*/
SPECTRUM_DEVICE_FUNC inline int DOF_2D(int order)
{
   return (order + 1) * (order + 2) / 2;
};

/*!
\brief Returns the number of degrees of freeedom in three dimensions
\author Vladimir Florinski
\date 07/17/2019
\param[in] order The order of reconstruction polynomial
\return The number of DOF
*/
SPECTRUM_DEVICE_FUNC inline int DOF_3D(int order)
{
   return (order + 1) * (order + 2) * (order + 3) / 6;
};

//! Binomial coefficents, power 0
SPECTRUM_DEVICE_FUNC const int binomial0[1] = {1};

//! Binomial coefficents, power 1
SPECTRUM_DEVICE_FUNC const int binomial1[2] = {1, 1};

//! Binomial coefficents, power 2
SPECTRUM_DEVICE_FUNC const int binomial2[3] = {1, 2, 1};

//! Binomial coefficents, power 3
SPECTRUM_DEVICE_FUNC const int binomial3[4] = {1, 3, 3, 1};

//! Binomial coefficents, power 4
SPECTRUM_DEVICE_FUNC const int binomial4[5] = {1, 4, 6, 4, 1};

//! Binomial coefficents, power 5
const int binomial5[6] = {1, 5, 10, 10, 5, 1};

//! Assembled Pascal triange
const int* const binomial[MONO_ORDER_HIGH + 1] = {binomial0, binomial1, binomial2, binomial3, binomial4, binomial5};

//! Bitwise mask for x
const uint32_t xbitfield = 255 << 16;

//! Bitwise mask for y
const uint32_t ybitfield = 255 << 8;

//! Bitwise mask for z
const uint32_t zbitfield = 255;

//! Bitwise encoding of the monomial indices using reverse colexicographical ordering
const uint32_t moment_pl[]
= {0x000000,  //           1
   0x010000,  // x
   0x000100,  // y
   0x000001,  // z         4
   0x020000,  // xx
   0x010100,  // xy
   0x010001,  // xz
   0x000200,  // yy
   0x000101,  // yz
   0x000002,  // zz       10
   0x030000,  // xxx
   0x020100,  // xxy
   0x020001,  // xxz
   0x010200,  // xyy
   0x010101,  // xyz
   0x010002,  // xzz
   0x000300,  // yyy
   0x000201,  // yyz
   0x000102,  // yzz
   0x000003,  // zzz      20
   0x040000,  // xxxx
   0x030100,  // xxxy
   0x030001,  // xxxz
   0x020200,  // xxyy
   0x020101,  // xxyz
   0x020002,  // xxzz
   0x010300,  // xyyy
   0x010201,  // xyyz
   0x010102,  // xyzz
   0x010003,  // xzzz
   0x000400,  // yyyy
   0x000301,  // yyyz
   0x000202,  // yyzz
   0x000103,  // yzzz
   0x000004,  // zzzz     35
   0x050000,  // xxxxx
   0x040100,  // xxxxy
   0x040001,  // xxxxz
   0x030200,  // xxxyy
   0x030101,  // xxxyz
   0x030002,  // xxxzz
   0x020300,  // xxyyy
   0x020201,  // xxyyz
   0x020102,  // xxyzz
   0x020003,  // xxzzz
   0x010400,  // xyyyy
   0x010301,  // xyyyz
   0x010202,  // xyyzz
   0x010103,  // xyzzz
   0x010004,  // xzzzz
   0x000500,  // yyyyy
   0x000401,  // yyyyz
   0x000302,  // yyyzz
   0x000203,  // yyzzz
   0x000104,  // yzzzz
   0x000005}; // zzzzz    56

//! Inverse moment lookup table
const int moment_lookup[MONO_ORDER_HIGH + 1][MONO_ORDER_HIGH + 1][MONO_ORDER_HIGH + 1]
= {{{ 0,  1,  4, 10, 20, 35},   // y^0 z^0
    { 2,  5, 11, 21, 36, -1},   // y^1 z^0
    { 7, 13, 23, 38, -1, -1},   // y^2 z^0
    {16, 26, 41, -1, -1, -1},   // y^3 z^0
    {30, 45, -1, -1, -1, -1},   // y^4 z^0
    {50, -1, -1, -1, -1, -1}},  // y^5 z^0

   {{ 3,  6, 12, 22, 37, -1},   // y^0 z^1
    { 8, 14, 24, 39, -1, -1},   // y^1 z^1
    {17, 27, 42, -1, -1, -1},   // y^2 z^1
    {31, 46, -1, -1, -1, -1},   // y^3 z^1
    {51, -1, -1, -1, -1, -1},   // y^4 z^1
    {-1, -1, -1, -1, -1, -1}},  // y^5 z^1

   {{ 9, 15, 25, 40, -1, -1},   // y^0 z^2
    {18, 28, 43, -1, -1, -1},   // y^1 z^2
    {32, 47, -1, -1, -1, -1},   // y^2 z^2
    {52, -1, -1, -1, -1, -1},   // y^3 z^2
    {-1, -1, -1, -1, -1, -1},   // y^4 z^2
    {-1, -1, -1, -1, -1, -1}},  // y^5 z^2

   {{19, 29, 44, -1, -1, -1},   // y^0 z^3
    {33, 48, -1, -1, -1, -1},   // y^1 z^3
    {53, -1, -1, -1, -1, -1},   // y^2 z^3
    {-1, -1, -1, -1, -1, -1},   // y^3 z^3
    {-1, -1, -1, -1, -1, -1},   // y^4 z^3
    {-1, -1, -1, -1, -1, -1}},  // y^5 z^3

   {{34, 49, -1, -1, -1, -1},   // y^0 z^4
    {54, -1, -1, -1, -1, -1},   // y^1 z^4
    {-1, -1, -1, -1, -1, -1},   // y^2 z^4
    {-1, -1, -1, -1, -1, -1},   // y^3 z^4
    {-1, -1, -1, -1, -1, -1},   // y^4 z^4
    {-1, -1, -1, -1, -1, -1}},  // y^5 z^4

   {{55, -1, -1, -1, -1, -1},   // y^0 z^5
    {-1, -1, -1, -1, -1, -1},   // y^1 z^5
    {-1, -1, -1, -1, -1, -1},   // y^2 z^5
    {-1, -1, -1, -1, -1, -1},   // y^3 z^5
    {-1, -1, -1, -1, -1, -1},   // y^4 z^5
    {-1, -1, -1, -1, -1, -1}}}; // y^5 z^5

//! Returns a product of powers of x, y, and z using bitfield input
SPECTRUM_DEVICE_FUNC double CoordPower3D(const GeoVector& v, uint32_t pl);

//! Returns a power of r
SPECTRUM_DEVICE_FUNC double CoordPower3D(double r, uint32_t pl);

//! Returns a product of powers of x, y, and z
SPECTRUM_DEVICE_FUNC double CoordPower3D(const GeoVector& v, int l, int m, int n);

//! Evaluates a polynomial of one argument
SPECTRUM_DEVICE_FUNC double EvalPoly1D(int p, const double* c, double x);

//! Evaluates a polynomial of two arguments
SPECTRUM_DEVICE_FUNC double EvalPoly2D(int p, const double* c, double x, double y);

//! Evaluates a polynomial of three arguments
SPECTRUM_DEVICE_FUNC double EvalPoly3D(int p, const double* c, double x, double y, double z);

};

#endif
