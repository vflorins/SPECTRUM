/*!
\file definitions.hh
\brief Defines some useful macros, memory management, math, and search routines
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DEFINITIONS_HH
#define SPECTRUM_DEFINITIONS_HH

#include "common/gpu_config.hh"
#include <cmath>

namespace Spectrum {

#define LEFT  0
#define RIGHT 1

#define SZINT sizeof(int)
#define SZDBL sizeof(double)

#define GEO_INCR(var, T) var = static_cast<T>(var + 1) 
#define GEO_DECR(var, T) var = static_cast<T>(var - 1) 

#define is_odd(idx)   ((idx) & 1)
#define is_even(idx) !((idx) & 1)

//! +1 if the number is odd and -1 if the number is even
#define ODDPLUS(idx) ((idx) & 1 ? 1.0 : -1.0)

//! Raise given bits
#define RAISE_BITS(value, bits) ((value) |= (bits))

//! Lower given bits
#define LOWER_BITS(value, bits) ((value) &= ~(bits))

//! Keep the bits present in "bits" and lower the rest
#define KEEP_BITS(value, bits) ((value) &= (bits))

//! Bits are raised
#define BITS_RAISED(value, bits) ((value) & (bits))

//! Bits are lowered
#define BITS_LOWERED(value, bits) (!((value) & (bits)))

//! A large number, according to the code
#define large 1.0E20

//! A small number, according to the code
#define little 1.0E-5

//! A smaller number, according to the code
#define small 1.0E-8

//! An even smaller number, according to the code
#define miniscule 1.0E-12

//! A very small number, according to the code
#define tiny 1.0E-15

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Global floating point constants that are also available on the device
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! \f$2\pi\f$
const double twopi = 2.0 * M_PI;

//! \f$4\pi\f$
const double fourpi = 4.0 * M_PI;

//! \f$8\pi\f$
const double eightpi = 8.0 * M_PI;

//! \f$\sqrt{\pi}\f$
const double sqrtpi = sqrt(M_PI);

//! \f$\sqrt{2}\f$
const double sqrttwo = sqrt(2.0);

//! \f$\sqrt{3}\f$
const double sqrtthr = sqrt(3.0);

//! \f$1/\sqrt{3}\f$
const double oo_sqrtthr = 1.0 / sqrtthr;

//! \f$\sqrt{5}\f$
const double sqrtfiv = sqrt(5.0);

//! \f$\sqrt{6}\f$
const double sqrtsix = sqrt(6.0);

//! \f$\sqrt{7}\f$
const double sqrtsev = sqrt(7.0);

//! \f$\sqrt{8}\f$
const double sqrteig = sqrt(8.0);

//! \f$\sqrt{10}\f$
const double sqrtten = sqrt(10.0);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Some inlined functions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Converts from radians to degrees
\author Vladimir Florinski
\date 02/18/2018
\param[in] rad Angle in radians
\return Angle in degrees
*/
SPECTRUM_DEVICE_FUNC inline double RadToDeg(double rad)
{
   return rad * 180.0 / M_PI;
};

/*!
\brief Convert from degrees to radians
\author Vladimir Florinski
\date 02/18/2018
\param[in] rad Angle in degrees
\return Angle in radians
*/
SPECTRUM_DEVICE_FUNC inline double DegToRad(double deg)
{
   return deg / 180.0 * M_PI;
};

/*!
\brief Power of two
\author Vladimir Florinski
\date 07/28/2016
\param[in] n The exponent
\return \f$2^n\f$
*/
SPECTRUM_DEVICE_FUNC inline int Pow2(int n) {return 1 << n;};

/*!
\brief Checks whether an integer is a power of two
\author Sean Eron Anderson
\author Vladimir Florinski
\date 05/14/2021
\param[in] n The integer to check
\return True if the number is a power of 2

Reference: http://www.graphics.stanford.edu/~seander/bithacks.htm
*/
SPECTRUM_DEVICE_FUNC inline bool IsPow2(int n)
{
   return !(n & (n - 1));
};

/*!
\brief Integer logarithm of base 2
\author Sean Eron Anderson
\author Vladimir Florinski
\date 05/14/2021
\param[in] n The integer to find logarithm of
\return \f$\log_2 n\f$

Reference: http://www.graphics.stanford.edu/~seander/bithacks.htm
*/
SPECTRUM_DEVICE_FUNC inline int Log2(int n)
{
   int log2 = 0;
   while(n >>= 1) log2++;
   return log2;
};

/*!
\brief Computes the square
\author Vladimir Florinski
\date 07/28/2016
\param[in] x The argument
\return \f$x^2\f$
*/
template <typename T> SPECTRUM_DEVICE_FUNC inline T Sqr(T x)
{
   return x * x;
};

/*!
\brief Computes the cube
\author Vladimir Florinski
\date 07/28/2016
\param[in] x The argument
\return \f$x^3\f$
*/
template <typename T> SPECTRUM_DEVICE_FUNC inline T Cube(T x)
{
   return x * x * x;
};

/*!
\brief Computes the fourth power
\author Vladimir Florinski
\date 07/22/2019
\param[in] x The argument
\return \f$x^4\f$
*/
template <typename T> SPECTRUM_DEVICE_FUNC inline T Quad(T x)
{
   return x * x * x * x;
};

/*!
\brief Computes the fifth power
\author Vladimir Florinski
\date 09/08/2023
\param[in] x The argument
\return \f$x^5\f$
*/
template <typename T> SPECTRUM_DEVICE_FUNC inline T Penta(T x)
{
   return x * x * x * x * x;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Some non-inlined functions
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Creates a two-dimensional array of arbitrary type
template <typename T> SPECTRUM_DEVICE_FUNC T** Create2D(int n, int m);

//! Creates a two-dimensional array of arbitrary type with index offset
template <typename T> SPECTRUM_DEVICE_FUNC T** Create2D_IO(int n, int m, int fidx);

//! Releases memory allocated by Create2D()
template <typename T> SPECTRUM_DEVICE_FUNC void Delete2D(T** array);

//! Creates pointers into a two-dimensional array previously allocated
template <typename T> SPECTRUM_DEVICE_FUNC T** Map2D(int n, int m, T* start, int fidx);

//! Deletes the pointers used by Map2D()
template <typename T> SPECTRUM_DEVICE_FUNC void Unmap2D(T** array);

//! Integer power
template <typename T> SPECTRUM_DEVICE_FUNC T IntPow(T x, int n);

//! Minmod operation
template <typename T> SPECTRUM_DEVICE_FUNC T MinMod(T x, T y);

//! Find a matching element in an array
template <typename T> SPECTRUM_DEVICE_FUNC int InList(int size, const T* array, T val);

//! Find the interval in an ascending array containing the given number
template <typename T> SPECTRUM_DEVICE_FUNC int LocateInArray(int l1, int l2, const T* array, T val, bool limit = false);

//! Solve a quadratic equation (real roots only)
template <typename T> SPECTRUM_DEVICE_FUNC bool QuadraticSolve(T a, T b, T c, T& x1, T& x2);

//! Solve a cubic equation (real roots only)
template <typename T> SPECTRUM_DEVICE_FUNC bool CubicSolve(T a, T b, T c, T d, T& x1, T& x2, T& x3);

//! Reduce a periodic argument to lie between 0 and "period"
template <typename T> SPECTRUM_DEVICE_FUNC T MakePeriodic(T x, T period);

};

#endif
