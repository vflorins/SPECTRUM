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

#define SZGM sizeof(GeoMatrix)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoMatrix class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A 3x3 component matrix class
\author Vladimir Florinski

A class representing a square 3x3 matrix used in vector transformations
*/
struct GeoMatrix {

//! Storage
   union {
      GeoVector row[3];
      double linear[9];
   };

//! Default constructor
   GeoMatrix(void);

//! Constructor from row vectors
   GeoMatrix(const GeoVector& v1, const GeoVector& v2, const GeoVector& v3);

//! Copy constructor
   GeoMatrix(const GeoMatrix& matr_r);

//! Access to the data for reading
   const double* Data(void) const;

//! Access to the data for writing
   double* Data(void);

//! Trace of the matrix
   double Trace(void) const;

//! Access to the rows as a vector array
   const GeoVector* VectorArray(void) const;

//! Access to rows for reading
   const GeoVector& operator [](int i) const;

//! Access to rows for writing
   GeoVector& operator [](int i);

//! Assignment operator from another matrix
   GeoMatrix& operator =(const GeoMatrix& matr_r);

//! Set all three components to the given value
   GeoMatrix& operator =(double val);

//! Transpose this matrix
   GeoMatrix& Transpose(void);

//! Add another matrix to this
   GeoMatrix& operator +=(const GeoMatrix& matr_r);

//! Subtract another matrix from this
   GeoMatrix& operator -=(const GeoMatrix& matr_r);

//! Multiply this matrix by a scalar from the left
   GeoMatrix& operator *=(double sclr_r);

//! Divide this matrix by a scalar
   GeoMatrix& operator /=(double sclr_r);

//! Make a dyadic product out of the same vector
   void Dyadic(const GeoVector& vect);

//! Make a dyadic product out of two vectors
   void Dyadic(const GeoVector& vect_l, const GeoVector& vect_r);

//! Add to a correlation matrix
   void IncrementCov(const GeoVector& vect_l, const GeoVector& vect_r);

//! Generate a basis with the z-axis along the given direction
   void AxisymmetricBasis(const GeoVector& ez);

//! Convert components from the standard basis to a different basis
   void ChangeToBasis(const GeoMatrix& basis);

//! Convert components to the standard basis from a different basis
   void ChangeFromBasis(const GeoMatrix& basis);

//! Compute the determinant
   double Det(void) const;

//! Compute a minor
   double Minor(int i, int j) const;

//! Compute the inverse
   GeoMatrix Inverse(void) const;

//! Compute the eigenvalues
   GeoVector Eigenvalues(void) const;

//! Compute the eigenvalues and left eigenvectors
   GeoVector Eigensystem(GeoMatrix& evec) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Negative of a matrix
   friend GeoMatrix operator -(const GeoMatrix& matr);

//! Add two matrices together
   friend GeoMatrix operator +(const GeoMatrix& matr_l, const GeoMatrix& matr_r);

//! Subtract one matrix from another
   friend GeoMatrix operator -(const GeoMatrix& matr_l, const GeoMatrix& matr_r);

//! Multiply a matrix by a scalar from the right
   friend GeoMatrix operator *(const GeoMatrix& matr_l, double sclr_r);

//! Multiply a matrix by a scalar from the left
   friend GeoMatrix operator *(double sclr_l, const GeoMatrix& matr_r);

//! Divide a matrix by a scalar
   friend GeoMatrix operator /(const GeoMatrix& matr_l, double sclr_r);

//! Return a product of a matrix and a column vectors
   friend GeoVector operator *(const GeoMatrix& matr_l, const GeoVector& vect_r);

//! Return a product of a row vector and a matrix
   friend GeoVector operator *(const GeoVector& vect_l, const GeoMatrix& matr_r);

//! Return a product of two matrices
   friend GeoMatrix operator *(const GeoMatrix& matr_l, const GeoMatrix& matr_r);

//! Return a vector product of a vector and a matrix
   friend GeoMatrix operator ^(const GeoVector& vect_l, const GeoMatrix& matr_r);

//! Return a vector product of a matrix and a vector
   friend GeoMatrix operator ^(const GeoMatrix& matr_l, const GeoVector& vect_r);

//! Inner product between two matrices together
   friend double operator %(const GeoMatrix& matr_l, const GeoMatrix& matr_r);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoMatrix inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/02/2021
*/
inline GeoMatrix::GeoMatrix(void)
{
}

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] v1 First row
\param[in] v2 Second row
\param[in] v3 Third row
*/
inline GeoMatrix::GeoMatrix(const GeoVector& v1, const GeoVector& v2, const GeoVector& v3) : row{v1, v2, v3}
{
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] vect_r Matrix to create a copy of
*/
inline GeoMatrix::GeoMatrix(const GeoMatrix& matr_r)
{
   operator =(matr_r);
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return "data" array
*/
inline const double* GeoMatrix::Data(void) const
{
   return linear;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return "data" array
*/
inline double* GeoMatrix::Data(void)
{
   return linear;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\return Sum of diagonal components
*/
inline double GeoMatrix::Trace(void) const
{
   return linear[0] + linear[4] + linear[8];
};

/*!
\author Vladimir Florinski
\date 08/15/2022
\return Access to the storage as an array of type GeoVector
*/
inline const GeoVector* GeoMatrix::VectorArray(void) const
{
   return row;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] i The desired row
\return \f$M_i\f$
*/
inline const GeoVector& GeoMatrix::operator [](int i) const
{
   return row[i];
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] i The desired row
\return \f$M_i\f$
*/
inline GeoVector& GeoMatrix::operator [](int i)
{
   return row[i];
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_r A matrix that will be copied into this matrix
\return Reference to this object
*/
inline GeoMatrix& GeoMatrix::operator =(const GeoMatrix& matr_r)
{
   if(this != &matr_r) memcpy(linear, matr_r.linear, 9 * SZDBL);
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] val A number to be assigned to all nine components
\return Reference to this object
*/
inline GeoMatrix& GeoMatrix::operator =(double val)
{
   row[0] = row[1] = row[2] = val;
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\return Reference to this object
*/
inline GeoMatrix& GeoMatrix::Transpose(void)
{
   std::swap(linear[1], linear[3]);
   std::swap(linear[2], linear[6]);
   std::swap(linear[5], linear[7]);
   return *this;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_r right operand \f$\mathbf{M}_1\f$
\return \f$\mathbf{M}+\mathbf{M}_1\f$
*/
inline GeoMatrix& GeoMatrix::operator +=(const GeoMatrix& matr_r)
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
inline GeoMatrix& GeoMatrix::operator -=(const GeoMatrix& matr_r)
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
inline GeoMatrix& GeoMatrix::operator *=(double sclr_r)
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
inline GeoMatrix& GeoMatrix::operator /=(double sclr_r)
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
inline double GeoMatrix::Det(void) const
{
   return row[0][0] * row[1][1] * row[2][2] + row[0][1] * row[1][2] * row[2][0] + row[0][2] * row[1][0] * row[2][1]
        - row[0][0] * row[1][2] * row[2][1] - row[0][1] * row[1][0] * row[2][2] - row[0][2] * row[1][1] * row[2][0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Friend methods of GeoMatrix
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/12/2023
\param[in] matr Matrix to multiply by -1 \f$\mathbf{M}\f$
\return \f$-\mathbf{M}\f$
*/
inline GeoMatrix operator -(const GeoMatrix& matr)
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
inline GeoMatrix operator +(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
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
inline GeoMatrix operator -(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
{
   GeoMatrix matr_tmp(matr_l);
   matr_tmp.row[0] -= matr_r.row[0];
   matr_tmp.row[1] -= matr_r.row[1];
   matr_tmp.row[2] -= matr_r.row[2];
   return matr_tmp;
};

/*!
\author Vladimir Florinski
\date 11/02/2021
\param[in] matr_l left operand \f$\mathbf{M}\f$
\param[in] sclr_r right operand \f$a\f$
\return \f$a\mathbf{M}\f$
*/
inline GeoMatrix operator *(const GeoMatrix& matr_l, double sclr_r)
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
inline GeoMatrix operator *(double sclr_l, const GeoMatrix& matr_r)
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
\return \f$a^{-1}\mathbf{v}\f$
*/
inline GeoMatrix operator /(const GeoMatrix& matr_l, double sclr_r)
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
inline GeoVector operator *(const GeoMatrix& matr_l, const GeoVector& vect_r)
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
inline GeoVector operator *(const GeoVector& vect_l, const GeoMatrix& matr_r)
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
inline GeoMatrix operator *(const GeoMatrix& matr_l, const GeoMatrix& matr_r)
{
   int i, j;
   GeoMatrix matr_tmp, matr_rt(matr_r);
   matr_rt.Transpose();
   for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
         matr_tmp.row[i][j] = matr_l[i] * matr_rt[j];
      };
   };
   return matr_tmp;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on matrices
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Compute a covariance matrix
GeoMatrix CovMatrix(const GeoVector& vect_l, const GeoVector& vect_r);

//! Make a matrix object out of components of a vector
GeoMatrix Dyadic(const GeoVector& vect);

//! Stream insertion operator
std::ostream& operator <<(std::ostream& os, const GeoMatrix& matr_r);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Types and constants
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! The unit matrix
const GeoMatrix gm_unit(gv_nx, gv_ny, gv_nz);

//! A matrix with zero components
const GeoMatrix gm_zeros(gv_zeros, gv_zeros, gv_zeros);

//! A matrix with all components equal to one
const GeoMatrix gm_ones(gv_ones, gv_ones, gv_ones);

};

#endif
