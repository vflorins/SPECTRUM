/*!
\file matrix.hh
\brief Declares and defines a 3x3 matrix class
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MATRIX_HH
#define SPECTRUM_MATRIX_HH

#include "vectors.hh"

namespace Spectrum {

//! Size of this class (should be 72)
#define SZGM sizeof(GeoMatrix)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoMatrix class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A 3x3 component matrix class
\author Vladimir Florinski
\author Juan G Alonso Guzman

A class representing a square 3x3 matrix used in vector transformations
*/
struct GeoMatrix {

//! Storage
   union {
      GeoVector row[3];
      double linear[9];
   };

//! Default constructor
   SPECTRUM_DEVICE_FUNC constexpr GeoMatrix(void) {};

//! Constructor from a single value
   SPECTRUM_DEVICE_FUNC explicit constexpr GeoMatrix(double val);

//! Constructor from row vectors
   SPECTRUM_DEVICE_FUNC constexpr GeoMatrix(const GeoVector& v1, const GeoVector& v2, const GeoVector& v3);

//! Access to the data for reading
   SPECTRUM_DEVICE_FUNC const double* Data(void) const;

//! Access to the data for writing
   SPECTRUM_DEVICE_FUNC double* Data(void);

//! Trace of the matrix
   SPECTRUM_DEVICE_FUNC double Trace(void) const;

//! Access to the rows as a vector array
   SPECTRUM_DEVICE_FUNC const GeoVector* VectorArray(void) const;

//! Access to rows for reading
   SPECTRUM_DEVICE_FUNC const GeoVector& operator [](int i) const;

//! Access to rows for writing
   SPECTRUM_DEVICE_FUNC GeoVector& operator [](int i);

//! Set all nine components to the given value
   SPECTRUM_DEVICE_FUNC constexpr GeoMatrix& operator =(double val);

//! Transpose this matrix
   SPECTRUM_DEVICE_FUNC GeoMatrix& Transpose(void);

//! Set this matrix to the transpose of a given matrix
   SPECTRUM_DEVICE_FUNC void Transpose(const GeoMatrix& matr_in);

//! Add another matrix to this
   SPECTRUM_DEVICE_FUNC GeoMatrix& operator +=(const GeoMatrix& matr_r);

//! Subtract another matrix from this
   SPECTRUM_DEVICE_FUNC GeoMatrix& operator -=(const GeoMatrix& matr_r);

//! Multiply this matrix by a scalar from the right
   SPECTRUM_DEVICE_FUNC GeoMatrix& operator *=(double sclr_r);

//! Divide this matrix by a scalar
   SPECTRUM_DEVICE_FUNC GeoMatrix& operator /=(double sclr_r);

//! Compute the determinant
   SPECTRUM_DEVICE_FUNC double Det(void) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Make a dyadic product out of the same vector
   SPECTRUM_DEVICE_FUNC void Dyadic(const GeoVector& vect);

//! Make a dyadic product out of two vectors
   SPECTRUM_DEVICE_FUNC void Dyadic(const GeoVector& vect_l, const GeoVector& vect_r);

//! Add to a correlation matrix
   SPECTRUM_DEVICE_FUNC void IncrementCov(const GeoVector& vect_l, const GeoVector& vect_r);

//! Generate a basis with the z-axis along the given direction
   SPECTRUM_DEVICE_FUNC void AxisymmetricBasis(const GeoVector& ez);

//! Convert rows from the standard basis to a different basis
   SPECTRUM_DEVICE_FUNC void ChangeToBasis(const GeoMatrix& basis);

//! Convert components from the standard basis to a different basis
   SPECTRUM_DEVICE_FUNC void ChangeToBasis(const GeoVector* basis);

//! Convert rows to the standard basis from a different basis
   SPECTRUM_DEVICE_FUNC void ChangeFromBasis(const GeoMatrix& basis);

//! Convert components to the standard basis from a different basis
   SPECTRUM_DEVICE_FUNC void ChangeFromBasis(const GeoVector* basis);

//! Compute a minor
   SPECTRUM_DEVICE_FUNC double Minor(int i, int j) const;

//! Compute the inverse
   SPECTRUM_DEVICE_FUNC GeoMatrix Inverse(void) const;

//! Compute the eigenvalues
   SPECTRUM_DEVICE_FUNC GeoVector Eigenvalues(void) const;

//! Compute the eigenvalues and left eigenvectors
   SPECTRUM_DEVICE_FUNC GeoVector Eigensystem(GeoMatrix& evec) const;

//! Complete the matrix to make it symmetric
   SPECTRUM_DEVICE_FUNC void CompleteSymmetric(void);

//! Complete the matrix to make it anti-symmetric
   SPECTRUM_DEVICE_FUNC void CompleteAntiSymmetric(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoMatrix inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/29/2024
\param[in] a Number to be asigned to each index
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoMatrix::GeoMatrix(double val)
{
   operator =(val);
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] v1 First row
\param[in] v2 Second row
\param[in] v3 Third row
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoMatrix::GeoMatrix(const GeoVector& v1, const GeoVector& v2, const GeoVector& v3)
                                               : row{v1, v2, v3}
{
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return "data" array
*/
SPECTRUM_DEVICE_FUNC inline const double* GeoMatrix::Data(void) const
{
   return linear;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return "data" array
*/
SPECTRUM_DEVICE_FUNC inline double* GeoMatrix::Data(void)
{
   return linear;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\return Sum of diagonal components
*/
SPECTRUM_DEVICE_FUNC inline double GeoMatrix::Trace(void) const
{
   return linear[0] + linear[4] + linear[8];
};

/*!
\author Vladimir Florinski
\date 08/15/2022
\return Access to the storage as an array of type GeoVector
*/
SPECTRUM_DEVICE_FUNC inline const GeoVector* GeoMatrix::VectorArray(void) const
{
   return row;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] i The desired row
\return \f$M_i\f$
*/
SPECTRUM_DEVICE_FUNC inline const GeoVector& GeoMatrix::operator [](int i) const
{
   return row[i];
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] i The desired row
\return \f$M_i\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoMatrix::operator [](int i)
{
   return row[i];
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] val A number to be assigned to all nine components
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoMatrix& GeoMatrix::operator =(double val)
{
   row[0] = row[1] = row[2] = val;
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix& GeoMatrix::Transpose(void)
{
   std::swap(linear[1], linear[3]);
   std::swap(linear[2], linear[6]);
   std::swap(linear[5], linear[7]);
   return *this;
};

/*!
\author Juan G Alonso Guzman
\date 06/12/2024
*/
SPECTRUM_DEVICE_FUNC inline void GeoMatrix::Transpose(const GeoMatrix& matr_in)
{
   if (this != &matr_in) memcpy(linear, matr_in.linear, 9 * SZDBL);
   this->Transpose();
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_r right operand \f$\mathbf{M}_1\f$
\return \f$\mathbf{M}+\mathbf{M}_1\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix& GeoMatrix::operator +=(const GeoMatrix& matr_r)
{
   row[0] += matr_r.row[0];
   row[1] += matr_r.row[1];
   row[2] += matr_r.row[2];
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_r right operand \f$\mathbf{M}_1\f$
\return \f$\mathbf{M}-\mathbf{M}_1\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix& GeoMatrix::operator -=(const GeoMatrix& matr_r)
{
   row[0] -= matr_r.row[0];
   row[1] -= matr_r.row[1];
   row[2] -= matr_r.row[2];
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] sclr_r right operand \f$a\f$
\return \f$a\mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix& GeoMatrix::operator *=(double sclr_r)
{
   row[0] *= sclr_r;
   row[1] *= sclr_r;
   row[2] *= sclr_r;
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] sclr_r right operand \f$a\f$
\return \f$a^{-1}\mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix& GeoMatrix::operator /=(double sclr_r)
{
   row[0] /= sclr_r;
   row[1] /= sclr_r;
   row[2] /= sclr_r;
   return *this;
};

/*!
\author Vladimir Florinski
\date 08/11/2022
\return The determinant of this matrix
*/
SPECTRUM_DEVICE_FUNC inline double GeoMatrix::Det(void) const
{
   return row[0][0] * row[1][1] * row[2][2] + row[0][1] * row[1][2] * row[2][0] + row[0][2] * row[1][0] * row[2][1]
        - row[0][0] * row[1][2] * row[2][1] - row[0][1] * row[1][0] * row[2][2] - row[0][2] * row[1][1] * row[2][0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on GeoMatrix
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/12/2023
\param[in] matr Matrix to multiply by -1 \f$\mathbf{M}\f$
\return \f$-\mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator -(const GeoMatrix& matr)
{
   GeoMatrix matr_tmp;
   matr_tmp.row[0] = -matr.row[0];
   matr_tmp.row[1] = -matr.row[1];
   matr_tmp.row[2] = -matr.row[2];
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_l left operand \f$\mathbf{M}_1\f$
\param[in] matr_r right operand \f$\mathbf{M}_2\f$
\return \f$\mathbf{M}_1+\mathbf{M}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator +(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] += matr_r.row[0];
   matr_tmp.row[1] += matr_r.row[1];
   matr_tmp.row[2] += matr_r.row[2];
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_l left operand \f$\mathbf{M}_1\f$
\param[in] matr_r right operand \f$\mathbf{M}_2\f$
\return \f$\mathbf{M}_1-\mathbf{M}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator -(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] -= matr_r.row[0];
   matr_tmp.row[1] -= matr_r.row[1];
   matr_tmp.row[2] -= matr_r.row[2];
   return matr_tmp;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] matr_l left operand \f$\mathbf{M}\f$
\param[in] sclr_r right operand \f$a\f$
\return \f$\mathbf{M} + a\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator +(const GeoMatrix& matr_l, double sclr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] += sclr_r;
   matr_tmp.row[1] += sclr_r;
   matr_tmp.row[2] += sclr_r;
   return matr_tmp;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] sclr_l left operand \f$a\f$
\param[in] matr_r right operand \f$\mathbf{v}\f$
\return \f$a + \mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator +(double sclr_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp(matr_r);
   matr_tmp.row[0] += sclr_l;
   matr_tmp.row[1] += sclr_l;
   matr_tmp.row[2] += sclr_l;
   return matr_tmp;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2025
\param[in] matr_l left operand \f$\mathbf{M}\f$
\param[in] sclr_r right operand \f$a\f$
\return \f$\mathbf{M} - a\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator -(const GeoMatrix& matr_l, double sclr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] -= sclr_r;
   matr_tmp.row[1] -= sclr_r;
   matr_tmp.row[2] -= sclr_r;
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_l left operand \f$\mathbf{M}\f$
\param[in] sclr_r right operand \f$a\f$
\return \f$a\mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator *(const GeoMatrix& matr_l, double sclr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] *= sclr_r;
   matr_tmp.row[1] *= sclr_r;
   matr_tmp.row[2] *= sclr_r;
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] sclr_l left operand \f$a\f$
\param[in] matr_r right operand \f$\mathbf{v}\f$
\return \f$a\mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator *(double sclr_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp(matr_r);
   matr_tmp.row[0] *= sclr_l;
   matr_tmp.row[1] *= sclr_l;
   matr_tmp.row[2] *= sclr_l;
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_l left operand \f$\mathbf{M}\f$
\param[in] sclr_r right operand \f$a\f$
\return \f$a^{-1}\mathbf{M}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator /(const GeoMatrix& matr_l, double sclr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] /= sclr_r;
   matr_tmp.row[1] /= sclr_r;
   matr_tmp.row[2] /= sclr_r;
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/20/2023
\param[in] matr_l left operand \f$\mathbf{M}_1\f$
\param[in] vect_r right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{M}_1\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator *(const GeoMatrix& matr_l, const GeoVector& vect_r)
{
   return GeoVector(matr_l.row[0] * vect_r, matr_l.row[1] * vect_r, matr_l.row[2] * vect_r);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/20/2023
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] matr_r right operand \f$\mathbf{M}_2\f$
\return \f$\mathbf{v}^T_1\mathbf{M}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator *(const GeoVector& vect_l, const GeoMatrix& matr_r)
{
   GeoVector vect_tmp;
   vect_tmp[0] = vect_l[0] * matr_r.row[0][0] + vect_l[1] * matr_r.row[1][0] + vect_l[2] * matr_r.row[2][0];
   vect_tmp[1] = vect_l[0] * matr_r.row[0][1] + vect_l[1] * matr_r.row[1][1] + vect_l[2] * matr_r.row[2][1];
   vect_tmp[2] = vect_l[0] * matr_r.row[0][2] + vect_l[1] * matr_r.row[1][2] + vect_l[2] * matr_r.row[2][2];
   return vect_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_l left operand \f$\mathbf{M}_1\f$
\param[in] matr_r right operand \f$\mathbf{M}_2\f$
\return \f$\mathbf{M}_1\mathbf{M}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator *(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp, matr_rt(matr_r);
   matr_rt.Transpose();
   for (auto i = 0; i < 3; i++) {
      for (auto j = 0; j < 3; j++) {
         matr_tmp.row[i][j] = matr_l[i] * matr_rt[j];
      };
   };
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 02/27/2023
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] matr_r right operand \f$\mathbf{M}_2\f$
\return \f$\mathbf{v}_1\times\mathbf{M}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator ^(const GeoVector& vect_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp;
   for (auto i = 0; i < 3; i++) matr_tmp[i] = vect_l ^ matr_r[i];
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 02/27/2023
\param[in] matr_l left operand \f$\mathbf{M}_1\f$
\param[in] vect_r right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{M}_1\times\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoMatrix operator ^(const GeoMatrix& matr_l, const GeoVector& vect_r)
{
   GeoMatrix matr_tmp;
   for (auto i = 0; i < 3; i++) matr_tmp[i] = matr_l[i] ^ vect_r;
   return matr_tmp;
};

/*!
\author Juan G Alonso Guzman
\date 10/06/2023
\param[in] matr_l left operand \f$\mathbf{M}_1\f$
\param[in] matr_r right operand \f$\mathbf{M}_2\f$
\return \f$\mathbf{M}_1 : \mathbf{M}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline double operator %(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
{
   double ip = 0.0;
   for (auto i = 0; i < 3; i++) ip += matr_l[i] * matr_r[i];
   return ip;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Compute a covariance matrix
SPECTRUM_DEVICE_FUNC GeoMatrix Dyadic(const GeoVector& vect_l, const GeoVector& vect_r);

//! Make a matrix object out of components of a vector
SPECTRUM_DEVICE_FUNC GeoMatrix Dyadic(const GeoVector& vect);

//! Stream insertion operator
std::ostream& operator <<(std::ostream& os, const GeoMatrix& matr_r);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Types and constants
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! The unit matrix
SPECTRUM_CONSTEXPR GeoMatrix gm_unit = {gv_nx, gv_ny, gv_nz};

//! A matrix with zero components
SPECTRUM_CONSTEXPR GeoMatrix gm_zeros = {gv_zeros, gv_zeros, gv_zeros};

//! A matrix with all components equal to one
SPECTRUM_CONSTEXPR GeoMatrix gm_ones = {gv_ones, gv_ones, gv_ones};

};

#endif
