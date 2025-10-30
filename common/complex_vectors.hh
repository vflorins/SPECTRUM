/*!
\file complex_vectors.hh
\brief Declares and defines three-component complex-valued vectors suitable for wave dispersion solvers
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COMPLEX_VECTORS_HH
#define SPECTRUM_COMPLEX_VECTORS_HH

#include "common/vectors.hh"

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
   SPECTRUM_DEVICE_FUNC constexpr ComplexVector(void) {};

//! Constructor from a single value
   SPECTRUM_DEVICE_FUNC explicit constexpr ComplexVector(std::complex<double> val);

//! Constructor from components
   SPECTRUM_DEVICE_FUNC constexpr ComplexVector(std::complex<double> x_in, std::complex<double> y_in, std::complex<double> z_in);

//! Constructor from a GeoVector
   SPECTRUM_DEVICE_FUNC constexpr ComplexVector(const GeoVector& other);

//! Computes the norm of this vector
   SPECTRUM_DEVICE_FUNC double Norm(void) const;

//! Get the real part
   SPECTRUM_DEVICE_FUNC GeoVector real(void) const;

//! Get the imaginary part
   SPECTRUM_DEVICE_FUNC GeoVector imag(void) const;

//! Convert components from the standard basis to a different basis
   SPECTRUM_DEVICE_FUNC void ChangeToBasis(const GeoVector* basis);

//! Convert components to the standard basis from a different basis
   SPECTRUM_DEVICE_FUNC void ChangeFromBasis(const GeoVector* basis);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ComplexVector inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/09/2024
\param[in] val Value to be asigned to each component
*/
SPECTRUM_DEVICE_FUNC inline constexpr ComplexVector::ComplexVector(std::complex<double> val)
                                                   : SimpleArray(val)
{
};

/*!
\author Vladimir Florinski
\date 12/09/2024
\param[in] x_in First component
\param[in] y_in Second component
\param[in] z_in Third component
*/
SPECTRUM_DEVICE_FUNC inline constexpr ComplexVector::ComplexVector(std::complex<double> x_in, std::complex<double> y_in, std::complex<double> z_in)
{
   x = x_in;
   y = y_in;
   z = z_in;
};

/*!
\author Vladimir Florinski
\date 12/06/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline constexpr ComplexVector::ComplexVector(const GeoVector& other)
{
   x = other.x;
   y = other.y;
   z = other.z;
};

/*!
\author Vladimir Florinski
\date 11/22/2024
\return \f$|v|\f$
*/
SPECTRUM_DEVICE_FUNC inline double ComplexVector::Norm(void) const
{
   return sqrt(std::norm(x) + std::norm(y) + std::norm(z));
};

/*!
\author Vladimir Florinski
\date 02/20/2025
\return Real part
*/
SPECTRUM_DEVICE_FUNC inline GeoVector ComplexVector::real(void) const
{
   return GeoVector(x.real(), y.real(), z.real());
};

/*!
\author Vladimir Florinski
\date 02/20/2025
\return Imaginary part
*/
SPECTRUM_DEVICE_FUNC inline GeoVector ComplexVector::imag(void) const
{
   return GeoVector(x.imag(), y.imag(), z.imag());
};

/*!
\author Vladimir Florinski
\date 02/21/2025
\param[in] basis Three new basis vectors
*/
SPECTRUM_DEVICE_FUNC inline void ComplexVector::ChangeToBasis(const GeoVector* basis)
{
   ComplexVector vect_dif;
   vect_dif.x = x * basis[0].x + y * basis[0].y + z * basis[0].z;
   vect_dif.y = x * basis[1].x + y * basis[1].y + z * basis[1].z;
   vect_dif.z = x * basis[2].x + y * basis[2].y + z * basis[2].z;
   std::memcpy(data, vect_dif.data, sizeof(ComplexVector));
};

/*!
\author Vladimir Florinski
\date 02/20/2025
\param[in] basis Three new basis vectors
*/
SPECTRUM_DEVICE_FUNC inline void ComplexVector::ChangeFromBasis(const GeoVector* basis)
{
   ComplexVector vect_std;
   vect_std.x = x * basis[0].x + y * basis[1].x + z * basis[2].x;
   vect_std.y = x * basis[0].y + y * basis[1].y + z * basis[2].y;
   vect_std.z = x * basis[0].z + y * basis[1].z + z * basis[2].z;
   std::memcpy(data, vect_std.data, sizeof(ComplexVector));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on vectors
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Multiply a real vector by a complex number (as left operand)
\author Vladimir Florinski
\date 12/09/2024
\param[in] a      Left operand \f$a\f$
\param[in] vect_r Right operand \f$\mathbf{v}\f$
\return \f$a\mathbf{v}\f$
*/
SPECTRUM_DEVICE_FUNC inline ComplexVector operator *(std::complex<double> a, const GeoVector& vect_r)
{
   ComplexVector retval(vect_r);
   retval.x *= a;
   retval.y *= a;
   retval.z *= a;
   return retval;
};

/*!
\brief Return a scalar product of a complex and a real vector
\author Vladimir Florinski
\date 02/20/2025
\param[in] cvec_l Left operand \f$\mathbf{v}_1\f$
\param[in] vec_r  Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1\cdot\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline std::complex<double> operator *(const ComplexVector& cvec_l, const GeoVector& vec_r)
{
   return cvec_l.x * vec_r.x + cvec_l.y * vec_r.y + cvec_l.z * vec_r.z;
};

/*!
\brief Compute a vector product of a complex and a real vector
\author Vladimir Florinski
\date 12/09/2024
\param[in] cvec_l Left operand \f$\mathbf{v}_1\f$
\param[in] vec_r  Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1\times\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline ComplexVector operator ^(const ComplexVector& cvec_l, const GeoVector& vec_r)
{
   ComplexVector retval;
   retval.x = cvec_l.y * vec_r.z - cvec_l.z * vec_r.y;
   retval.y = cvec_l.z * vec_r.x - cvec_l.x * vec_r.z;
   retval.z = cvec_l.x * vec_r.y - cvec_l.y * vec_r.x;
   return retval;
};

/*!
\brief Compute a null vector of a rank-deficient matrix
\author Vladimir Florinski
\date 12/09/2024
\param[in] M Three row vectors representing the matrix
\return Null vector of the matrix with unit complex amplitude
*/
SPECTRUM_DEVICE_FUNC inline ComplexVector NullVectorZ(const ComplexVector M[3])
{
   std::complex<double> m1, m2;
   ComplexVector null_vec(0.0, 0.0, 0.0);
   
// Generate the reduced second row vector (the first element is not used)
   m1 = M[1][1] - M[1][0] * M[0][1] / M[0][0];
   m2 = M[1][2] - M[1][0] * M[0][2] / M[0][0];

// Compute the null vector
   null_vec[2] = 1.0;
   null_vec[1] = -m2 / m1;
   null_vec[0] = -(M[0][2] + M[0][1] * null_vec[1]) / M[0][0];
   null_vec /= null_vec.Norm();

   return null_vec;
};

};

#endif
