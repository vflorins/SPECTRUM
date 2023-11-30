/*!
\file math_extra.hh
\brief Defines some advanced mathematical routines
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MATH_EXTRA_HH
#define SPECTRUM_MATH_EXTRA_HH

#include <cmath>

namespace Spectrum {

//! Creates a two-dimensional array of arbitrary type
template <typename T1, typename T2> T2 HypergeometricSeries1F2(T1 a, T1 b1, T1 b2, double z)

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
