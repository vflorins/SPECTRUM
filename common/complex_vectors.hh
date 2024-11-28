/*!
\file complex_vectors.hh
\brief Declares and defines three-component complex-valued vectors for wave dispersion solvers
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COMPLEX_VECTORS_HH
#define SPECTRUM_COMPLEX_VECTORS_HH

#include <complex>

#include "common/simple_array.hh"
#include "common/print_warn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ComplexVector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A skeleton complex vector class that obeys the algebraic rules for _real_ vectors, for use in dispersion solvers
\author Vladimir Florinski
*/
struct ComplexVector : public SimpleArray<std::complex<double>, 3>
{
//! Default constructor
   constexpr ComplexVector(void) {};

//! Constructor from a single value
   explicit constexpr ComplexVector(std::complex<double> a);

//! Constructor from components
   constexpr ComplexVector(std::complex<double> x_in, std::complex<double> y_in, std::complex<double> z_in);

//! Computes the norm of this vector
   std::complex<double> Norm(void) const;

//! Vector-multiply this by another vector
   ComplexVector& operator ^=(const ComplexVector& other);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ComplexVector inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/23/2024
\param[in] a Number to be asigned to each index
*/
inline constexpr ComplexVector::ComplexVector(std::complex<double> a)
                              : SimpleArray(a)
{
};

/*!
\author Vladimir Florinski
\date 11/22/2024
\param[in] x First component
\param[in] y Second component
\param[in] z Third component
*/
inline constexpr ComplexVector::ComplexVector(std::complex<double> x_in, std::complex<double> y_in, std::complex<double> z_in)
{
   x = x_in;
   y = y_in;
   z = z_in;
};

/*!
\author Vladimir Florinski
\date 11/22/2024
\return \f$|v|\f$
*/
inline std::complex<double> ComplexVector::Norm(void) const
{
   return sqrt(Sqr(x) + Sqr(y) + Sqr(z));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on vectors
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Return a scalar product of two vectors
\author Vladimir Florinski
\date 11/22/2024
\param[in] vect_l Left operand \f$\mathbf{v}_1\f$
\param[in] vect_r Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1\cdot\mathbf{v}_2\f$
*/
inline std::complex<double> operator *(const ComplexVector& vect_l, const ComplexVector& vect_r)
{
   return vect_l.ScalarProd(vect_r);
};

/*!
\brief Compute a vector product of two vectors
\author Vladimir Florinski
\date 11/22/2024
\param[in] vect_l Left operand \f$\mathbf{v}_1\f$
\param[in] vect_r Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1\times\mathbf{v}_2\f$
*/
inline ComplexVector operator ^(const ComplexVector& vect_l, const ComplexVector& vect_r)
{
   ComplexVector vect_tmp;
   vect_tmp.x = vect_l.y * vect_r.z - vect_l.z * vect_r.y;
   vect_tmp.y = vect_l.z * vect_r.x - vect_l.x * vect_r.z;
   vect_tmp.z = vect_l.x * vect_r.y - vect_l.y * vect_r.x;
   return vect_tmp;
};

/*!
\brief Compute a vector product of two vectors
\author Vladimir Florinski
\date 11/23/2024
\param[in] M Three row vectors representing the matrix
\return Null vector of the matrix
*/
inline ComplexVector NullVector(const ComplexVector M[3])
{
   std::complex<double> m1, m2;
   ComplexVector null_vec(0.0, 0.0, 0.0);
   
// Generate the reduced second row vector (the first element is not used)
   m1 = M[1][1] - M[1][0] * M[0][1] / M[0][0];
   m2 = M[1][2] - M[1][0] * M[0][2] / M[0][0];

// Test if the matrix is indeed rank deficient
   if (std::abs(M[2][2] - M[2][1] * m2 / m1) > sp_miniscule) PrintError(__FILE__, __LINE__, "Matrix not rank deficient", true);

// Compute the null vector
   else {
      null_vec[2] = 1.0;
      null_vec[1] = -m2 / m1;
      null_vec[0] = -(M[0][2] + M[0][1] * null_vec[1]) / M[0][0];
      null_vec /= null_vec.Norm();
   };

   return null_vec;
};

};

#endif
