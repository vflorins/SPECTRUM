/*!
\file definitions.cc
\brief Defines some useful macros, memory management, math, and search routines
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/definitions.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 07/28/2016
\param[in] n First dimension
\param[in] m Second dimension
\return Pointer to the storage
*/
template <typename T> SPECTRUM_DEVICE_FUNC T** Create2D(int n, int m)
{
   T** array = new T*[n];
   array[0] = new T[n * m];
   for(auto i = 1; i < n; i++) array[i] = array[i - 1] + m;
   return array;
};

/*!
\author Vladimir Florinski
\date 07/29/2019
\param[in] n    First dimension
\param[in] m    Second dimension
\param[in] fidx Index to start numbering elements (0 or 1)
\return Pointer to the storage
*/
template <typename T> SPECTRUM_DEVICE_FUNC T** Create2D_IO(int n, int m, int fidx)
{
   int offset = (fidx ? 1 : 0);
   T** array = new T*[n + offset];

// The array always has one extra element at the end, so it can be used with second index starting from 0 or 1.
   array[0] = new T[n * m + 1];
   if(offset) array[1] = array[0];
   for(auto i = 1 + offset; i < n + offset; i++) array[i] = array[i - 1] + m;
   return array;
};

/*!
\author Vladimir Florinski
\date 07/28/2016
\param[in,out] array Pointer to the storage
*/
template <typename T> SPECTRUM_DEVICE_FUNC void Delete2D(T** array)
{
   if(array) {
      delete[] array[0];
      delete[] array;
   };
};

/*!
\author Vladimir Florinski
\date 07/28/2016
\param[in] n     The first dimension
\param[in] m     The second dimension
\param[in] start Pointer to the start of the array
\param[in] fidx  Index to start numbering elements (0 or 1)
\return Pointer to the storage
*/
template <typename T> SPECTRUM_DEVICE_FUNC T** Map2D(int n, int m, T* start, int fidx)
{
   int offset = (fidx ? 1 : 0);
   T** array = new T*[n + offset];
   array[0] = start;
   if(offset) array[1] = array[0];
   for(auto i = 1 + offset; i < n + offset; i++) array[i] = array[i - 1] + m;
   return array;
};

/*!
\author Vladimir Florinski
\date 03/04/2024
\param[in] n First dimension
\param[in] m Second dimension
\param[in] l Third dimension
\return Pointer to the storage
*/
template <typename T> SPECTRUM_DEVICE_FUNC T*** Create3D(int n, int m, int l)
{
   T*** array = new T**[n];
   array[0] = new T*[n * m];
   array[0][0] = new T[n * m * l];
   for(auto i = 1; i < n; i++) array[i] = array[i - 1] + m;
   for(auto j = 1; j < n * m; j++) array[0][j] = array[0][j - 1] + l;
   return array;
};

/*!
\author Vladimir Florinski
\date 03/04/2024
\param[in,out] array Pointer to the storage
*/
template <typename T> SPECTRUM_DEVICE_FUNC void Delete3D(T*** array)
{
   if(array) {
      delete[] array[0][0];
      delete[] array[0];
      delete[] array;
   };
};

/*!
\author Vladimir Florinski
\date 02/18/2018
\param[in] x The argument
\param[in] n The exponent
\return \f$x^n\f$, or -1 if a negative power of zero is requested
*/
template <typename T> SPECTRUM_DEVICE_FUNC T IntPow(T x, int n)
{
   int i;
   T res = 1.0;

   if(n < 0) {
      if(x == 0) return -1.0;
      for(i = 0; i > n; i--) res /= x;
   }
   else {
      for(i = 0; i < n; i++) res *= x;
   };

   return res;
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] x First argument
\param[in] y Second argument
\return The result of a minmod operation
*/
template <typename T> SPECTRUM_DEVICE_FUNC T MinMod(T x, T y)
{
   if(x * y <= 0.0) return 0.0;
   else if(x * x > y * y) return y;
   else return x;
};

/*!
\author Vladimir Florinski
\date 11/16/2023
\param[in] x First argument
\param[in] y Second argument
\param[in] z Third argument
\return The result of a 3-argument minmod operation
*/
template <typename T> SPECTRUM_DEVICE_FUNC T MinMod(T x, T y, T z)
{
   T absx, absy, absz;
   absx = std::abs(x);
   absy = std::abs(y);
   absz = std::abs(z);
   if(x * y <= 0.0) return 0.0;
   if(absx > absy) {
      if(absy > absz) return z;
      else return y;
   }
   else {
      if(absx > absz) return z;
      else return x;
   };
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] size  Size of array
\param[in] array Array of numbers
\param[in] val   Number to match
\return Index in the array, or -1 if not found
*/
template <typename T> SPECTRUM_DEVICE_FUNC int InList(int size, const T* array, T val)
{
   int idx = 0;
   while(idx < size && array[idx] != val) idx++;
   return (idx == size ? -1 : idx);
};

/*!
\author Vladimir Florinski
\date 02/18/2018
\param[in] l1    Lower starting index
\param[in] l2    Upper starting index
\param[in] array Array of numbers in ascending order
\param[in] val   Value to locate in the array
\param[in] limit If true, return "l2" for val > array[l2]
\return Interval containing the value "val", -1 if outside the limits
*/
template <typename T> SPECTRUM_DEVICE_FUNC int LocateInArray(int l1, int l2, const T* array, T val, bool limit)
{
   if(l2 <= l1) return -1;
   if(val < array[l1]) return -1;
   if(val > array[l2]) {
      if(limit) return l2;
      else return -1;
   };

   int i1, i2, i3;
   i1 = i2 = l1;
   i3 = l2;

// Bisection algorithm
   while(i3 - i1 > 1) {
      i2 = (i1 + i3) >> 1;
      (val > array[i2] ? i1 : i3) = i2;
   };

// The interval has the same index as the _left_ interface
   return i1;
};

/*!
\author Vladimir Florinski
\date 08/11/2022
\param[in] a  Coefficient of x^2
\param[in] b  Coefficient of x^1
\param[in] c  Coefficient of x^0
\param[in] x1 First root
\param[in] x2 Second root
\return True if the roots are real
*/
template <typename T> SPECTRUM_DEVICE_FUNC bool QuadraticSolve(T a, T b, T c, T& x1, T& x2)
{
   T q;

   q = Sqr(b) - 4.0 * a * c;
   if(q >= 0.0) {
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
\author Vladimir Florinski
\date 08/11/2022
\param[in] a  Coefficient of x^3
\param[in] b  Coefficient of x^2
\param[in] c  Coefficient of x^1
\param[in] d  Coefficient of x^0
\param[in] x1 First root
\param[in] x2 Second root
\param[in] x2 Third root
\return True if all roots are real
*/
template <typename T> SPECTRUM_DEVICE_FUNC bool CubicSolve(T a, T b, T c, T d, T& x1, T& x2, T& x3)
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
   if(D <= 0.0) {
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
\author Vladimir Florinski
\date 08/22/2023
\param[in] x      A number to be reduced
\param[in] period The reduced number is between 0 and this number (must be positive)
\return Reduced number
*/
template <typename T> SPECTRUM_DEVICE_FUNC T MakePeriodic(T x, T period)
{
   int n = x / period;
   if(x < 0.0) return x - (n - 1) * period;
   else return x - n * period;
};

/*!
\author Vladimir Florinski
\date 02/15/2024
\param[in,out] tasks   Number of tasks to perform
\param[in,out] workers Number of workers available/required
\param[out]    tpw     Number of tasks per worker (except the last)
\param[out]    tplw    Number of tasks for the last worker
*/
template <typename T> SPECTRUM_DEVICE_FUNC void DistributeTasks(T& tasks, T& workers, T& tpw, T& tplw)
{
   T excess_tasks;

// Invalid input, but the return must make sense
   if((tasks <= 0) || (workers <= 0)) {
      tasks = 0;
      workers = 0;
      tpw = tplw = 0;
   }

// Workers exceed tasks
   else if(workers >= tasks) {
      workers = tasks;
      tpw = tplw = 1;
   }

// Normal situation
   else {

// Initial estimate
      tpw = tasks / workers;
      if(tasks % workers) tpw++;

// Redundant workers and the last worker load
      excess_tasks = tpw * workers - tasks;
      workers -= excess_tasks / tpw;
      tplw = tpw - excess_tasks % tpw;
   };
};

template SPECTRUM_DEVICE_FUNC double IntPow <double>(double x, int n);
template SPECTRUM_DEVICE_FUNC int InList <int>(int size, const int* array, int val);
template SPECTRUM_DEVICE_FUNC int LocateInArray <double>(int l1, int l2, const double* array, double val, bool limit);
template SPECTRUM_DEVICE_FUNC bool QuadraticSolve <double>(double a, double b, double c, double& x1, double& x2);
template SPECTRUM_DEVICE_FUNC bool CubicSolve <double>(double a, double b, double c, double d, double& x1, double& x2, double& x3);
template SPECTRUM_DEVICE_FUNC double MakePeriodic <double>(double x, double period);
template SPECTRUM_DEVICE_FUNC double** Create2D <double>(int n, int m);
template SPECTRUM_DEVICE_FUNC int** Create2D <int>(int n, int m);
template SPECTRUM_DEVICE_FUNC void Delete2D <double>(double** array);
template SPECTRUM_DEVICE_FUNC void Delete2D <int>(int** array);
template SPECTRUM_DEVICE_FUNC double*** Create3D <double>(int n, int m, int l);
template SPECTRUM_DEVICE_FUNC void Delete3D <double>(double*** array);
template SPECTRUM_DEVICE_FUNC double MinMod <double>(double x, double y);
template SPECTRUM_DEVICE_FUNC double MinMod <double>(double x, double y, double z);
template SPECTRUM_DEVICE_FUNC void DistributeTasks(int& tasks, int& workers, int& tpw, int& tplw);

};
