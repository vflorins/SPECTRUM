/*!
\file definitions.hh
\brief Defines some useful macros, memory management, math, and search routines
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DEFINITIONS_HH
#define SPECTRUM_DEFINITIONS_HH

#include <cmath>
#include <complex>

#include "config.h"
#include "common/gpu_config.hh"

namespace Spectrum {

//! An empty structure
struct EmptyStruct{};

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
#define sp_large 1.0E+20

//! A small number, according to the code
#define sp_little 1.0E-5

//! A smaller number, according to the code
#define sp_small 1.0E-8

//! An even smaller number, according to the code
#define sp_miniscule 1.0E-12

//! A very small number (machine precision), according to the code
#define sp_tiny 1.0E-15

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Global floating point constants that are also available on the device
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef M_PI
#define M_PI      3.141592653589793238462643383279502884E+0
#define M_PI_2    1.570796326794896619231321691639751442E+0
#define M_PI_4    7.853981633974483096156608458198757205E-1
#define M_SQRT2   1.414213562373095048801688724209698079E+0
#endif

#define M_PI_8    3.926990816987241548078304229099378605E-1
#define M_2PI     6.283185307179586476925286766559005768E+0
#define M_4PI     1.256637061435917295385057353311801154E+1
#define M_8PI     2.513274122871834590770114706623602307E+1
#define M_SQRTPI  1.772453850905516027298167483341145183E+0
#define M_SQRT3   1.732050807568877293527446341505872367E+0
#define M_SQRT1_3 5.773502691896257645091487805019574556E-1
#define M_SQRT5   2.236067977499789696409173668731276235E+0
#define M_SQRT6   2.449489742783178098197284074705891392E+0
#define M_SQRT7   2.645751311064590590501615753639260426E+0
#define M_SQRT8   2.828427124746190097603377448419396157E+0
#define M_SQRT10  3.162277660168379331998893544432718534E+0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Some commonly used functions
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
SPECTRUM_DEVICE_FUNC inline int Pow2(int n)
{
   return 1 << n;
};

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
   while (n >>= 1) log2++;
   return log2;
};

/*!
\brief Sign of a variable
\author Vladimir Florinski
\date 04/29/2024
\param[in] x The argument
\return Sign of x (-1 or 1, never 0)
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T sign(T x)
{
   return (x >= 0 ? 1 : -1);
};

/*!
\brief Computes the square
\author Vladimir Florinski
\date 07/28/2016
\param[in] x The argument
\return \f$x^2\f$
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T Sqr(T x)
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
template <typename T>
SPECTRUM_DEVICE_FUNC inline T Cube(T x)
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
template <typename T>
SPECTRUM_DEVICE_FUNC inline T Quad(T x)
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
template <typename T>
SPECTRUM_DEVICE_FUNC inline T Quint(T x)
{
   return x * x * x * x * x;
};

/*!
\brief Creates a two-dimensional array of arbitrary type
\author Vladimir Florinski
\date 07/28/2016
\param[in] n First dimension
\param[in] m Second dimension
\return Pointer to the storage
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T** Create2D(int n, int m)
{
   T** array = new T*[n];
   array[0] = new T[n * m];
   for (auto i = 1; i < n; i++) array[i] = array[i - 1] + m;
   return array;
};

/*!
\brief Releases memory allocated by Create2D()
\author Vladimir Florinski
\date 07/28/2016
\param[in,out] array Pointer to the storage
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline void Delete2D(T** array)
{
   if (array) {
      delete[] array[0];
      delete[] array;
   };
};

/*!
\brief Creates a three-dimensional array of arbitrary type
\author Vladimir Florinski
\date 03/04/2024
\param[in] n First dimension
\param[in] m Second dimension
\param[in] l Third dimension
\return Pointer to the storage
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T*** Create3D(int n, int m, int l)
{
   T*** array = new T**[n];
   array[0] = new T*[n * m];
   array[0][0] = new T[n * m * l];
   for (auto i = 1; i < n; i++) array[i] = array[i - 1] + m;
   for (auto j = 1; j < n * m; j++) array[0][j] = array[0][j - 1] + l;
   return array;
};

/*!
\brief Releases memory allocated by Create3D()
\author Vladimir Florinski
\date 03/04/2024
\param[in,out] array Pointer to the storage
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline void Delete3D(T*** array)
{
   if (array) {
      delete[] array[0][0];
      delete[] array[0];
      delete[] array;
   };
};

/*!
\brief Integer power
\author Vladimir Florinski
\date 02/18/2018
\param[in] x The argument
\param[in] n The exponent
\return \f$x^n\f$, or -1 if a negative power of zero is requested
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T IntPow(T x, int n)
{
   int i;
   T res = 1.0;

   if (n < 0) {
      if (x == 0) return -1.0;
      for (i = 0; i > n; i--) res /= x;
   }
   else {
      for (i = 0; i < n; i++) res *= x;
   };

   return res;
};

/*!
\brief Minmod operation
\author Vladimir Florinski
\date 07/22/2019
\param[in] x First argument
\param[in] y Second argument
\return The result of a minmod operation
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T MinMod(T x, T y)
{
   if (x * y <= 0.0) return 0.0;
   else if (x * x > y * y) return y;
   else return x;
};

/*!
\brief Minmod operation with three arguments
\author Vladimir Florinski
\date 11/16/2023
\param[in] x First argument
\param[in] y Second argument
\param[in] z Third argument
\return The result of a 3-argument minmod operation
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T MinMod(T x, T y, T z)
{
   T absx, absy, absz;
   absx = std::abs(x);
   absy = std::abs(y);
   absz = std::abs(z);
   if (x * y <= 0.0) return 0.0;
   if (absx > absy) {
      if (absy > absz) return z;
      else return y;
   }
   else {
      if (absx > absz) return z;
      else return x;
   };
};

/*!
\brief Find a matching element in an array
\author Vladimir Florinski
\date 07/22/2019
\param[in] size  Size of array
\param[in] array Array of numbers
\param[in] val   Number to match
\return Index in the array, or -1 if not found
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline int InList(int size, const T* array, T val)
{
   int idx = 0;
   while (idx < size && array[idx] != val) idx++;
   return (idx == size ? -1 : idx);
};

/*!
\brief Find the interval in an ascending array containing the given number
\author Vladimir Florinski
\date 02/18/2018
\param[in] l1    Lower starting index
\param[in] l2    Upper starting index
\param[in] array Array of numbers in ascending order
\param[in] val   Value to locate in the array
\param[in] limit If true, return "l2" for val > array[l2]
\return Interval containing the value "val", -1 if outside the limits
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline int LocateInArray(int l1, int l2, const T* array, T val, bool limit)
{
   if (l2 <= l1) return -1;
   if (val < array[l1]) return -1;
   if (val > array[l2]) {
      if (limit) return l2;
      else return -1;
   };

   int i1, i2, i3;
   i1 = i2 = l1;
   i3 = l2;

// Bisection algorithm
   while (i3 - i1 > 1) {
      i2 = (i1 + i3) >> 1;
      (val > array[i2] ? i1 : i3) = i2;
   };

// The interval has the same index as the _left_ interface
   return i1;
};

/*!
\brief Number of edges at a vertex in an infinite tesselation
\author Vladimir Florinski
\date 07/02/2024
\return Number of edges meeting at a vertex (3->6, 4->4, 6->3)
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC static constexpr int EdgesAtVert(void)
{
   return 2 * verts_per_face / (verts_per_face - 2);
};

/*!
\brief Side length of an inscribed polygon
\author Vladimir Florinski
\date 10/07/2025
\return Side length of a regular polygon inscribed in a unit circle 
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC static constexpr double InscribedPolygonSide(void)
{
   return 2.0 * sin(M_PI / verts_per_face);
};

/*!
\brief Solve a quadratic equation (real roots only)
\author Vladimir Florinski
\date 08/11/2022
\param[in] a  Coefficient of x^2
\param[in] b  Coefficient of x^1
\param[in] c  Coefficient of x^0
\param[in] x1 First root
\param[in] x2 Second root
\return True if the roots are real
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline bool QuadraticSolve(T a, T b, T c, T& x1, T& x2)
{
   T q = Sqr(b) - 4.0 * a * c;
   if (q >= 0.0) {
      q = sqrt(q);
      x1 = (-b + q) / (2.0 * a);
      x2 = (-b - q) / (2.0 * a);
      return true;
   }
   else {
     x1 = 0.0;
     x2 = 0.0;
     return false;
  };
};

/*!
\brief Solve a quadratic equation
\author Vladimir Florinski
\date 08/21/2024
\param[in] a  Coefficient of x^2
\param[in] b  Coefficient of x^1
\param[in] c  Coefficient of x^0
\param[in] x1 First root
\param[in] x2 Second root
\return True if the roots are real
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline bool QuadraticSolve(std::complex<T> a, std::complex<T> b, std::complex<T> c, std::complex<T>& x1, std::complex<T>& x2)
{
   std::complex<T> q = sqrt(Sqr(b) - 4.0 * a * c);
   x1 = (-b + q) / (2.0 * a);
   x2 = (-b - q) / (2.0 * a);
   return (q.imag() == 0.0);
};

/*!
\brief Solve a cubic equation
\author Vladimir Florinski
\date 08/21/2024
\param[in] a  Coefficient of x^3
\param[in] b  Coefficient of x^2
\param[in] c  Coefficient of x^1
\param[in] d  Coefficient of x^0
\param[in] x1 First root
\param[in] x2 Second root
\param[in] x2 Third root
\return True if all roots are real
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline bool CubicSolve(T a, T b, T c, T d, T& x1, T& x2, T& x3)
{
   T a1, b1, c1, p, q, Q, R, D, theta;

   a1 = b / a;
   b1 = c / a;
   c1 = d / a;
   p = b1 - a1 * a1 / 3.0;
   q = a1 * b1 / 3.0 - c1 - 2.0 * a1 * a1 * a1 / 27.0;
   Q = p / 3.0;
   R = q / 2.0;
   D = Q * Q * Q + R * R;
   if (D <= 0.0) {
      theta = acos(R / sqrt(-Q * Q * Q));
      Q = sqrt(-Q);
      x1 = 2.0 * Q * cos(theta / 3.0) - a1 / 3.0;
      x2 = 2.0 * Q * cos((theta + 2.0 * M_PI) / 3.0) - a1 / 3.0;
      x3 = 2.0 * Q * cos((theta + 4.0 * M_PI) / 3.0) - a1 / 3.0;
      return true;
   }
   else {
      x1 = 0.0;
      x2 = 0.0;
      x3 = 0.0;
      return false;
   };
};

/*!
\brief Solve a cubic equation (real roots only)
\author Vladimir Florinski
\date 08/21/2024
\param[in] a  Coefficient of x^3
\param[in] b  Coefficient of x^2
\param[in] c  Coefficient of x^1
\param[in] d  Coefficient of x^0
\param[in] x1 First root
\param[in] x2 Second root
\param[in] x2 Third root
\return True if all roots are real
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline bool CubicSolve(std::complex<T> a, std::complex<T> b, std::complex<T> c, std::complex<T> d,
                                            std::complex<T>& x1, std::complex<T>& x2, std::complex<T>& x3)
{
   std::complex<T> a1, b1, c1, p, q, Q, R, D, theta;

   a1 = b / a;
   b1 = c / a;
   c1 = d / a;
   p = b1 - a1 * a1 / 3.0;
   q = a1 * b1 / 3.0 - c1 - 2.0 * a1 * a1 * a1 / 27.0;
   Q = p / 3.0;
   R = q / 2.0;
   D = -Q * Q * Q + R * R;

   theta = acos(R / sqrt(-Q * Q * Q));
   Q = sqrt(-Q);
   x1 = 2.0 * Q * cos(theta / 3.0) - a1 / 3.0;
   x2 = 2.0 * Q * cos((theta + 2.0 * M_PI) / 3.0) - a1 / 3.0;
   x3 = 2.0 * Q * cos((theta + 4.0 * M_PI) / 3.0) - a1 / 3.0;

   return (sqrt(D).imag() == 0);
};

/*!
\brief Reduce a periodic argument to lie between 0 and "period"
\author Vladimir Florinski
\date 08/22/2023
\param[in] x      A number to be reduced
\param[in] period The reduced number is between 0 and this number (must be positive)
\return Reduced number
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T MakePeriodic(T x, T period)
{
   int n = x / period;
   if (x < 0.0) return x - (n - 1) * period;
   else return x - n * period;
};

/*!
\brief Average values in array using arithmetic or geometric mean
\author Juan G Alonso Guzman
\date 01/29/2025
\param[in] size  Size of array
\param[in] x     Array of values to average
\param[in] arith Whether to calculate the arithmetic or geometric average
\return Average of values in array
\note The geometric average is done by exponentiating the average of the logged values for numerical stability
*/
template <typename T>
SPECTRUM_DEVICE_FUNC inline T Average(int size, const T* x, bool arith)
{
   T avg = 0.0;
   for (auto i = 0; i < size; i++) avg += (arith ? x[i] : log(x[i]));
   avg /= size;
   return (arith ? avg : exp(avg));
};

};

#endif
