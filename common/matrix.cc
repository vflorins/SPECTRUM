/*!
\file matrix.cc
\brief Defines some non-trivial functions operating on matrices
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "matrix.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Methods of GeoMatrix
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/15/2022
\param[in] vect Vector
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::Dyadic(const GeoVector& vect)
{
   for (auto i = 0; i < 3; i++) {
      for (auto j = i; j < 3; j++) {
         row[i][j] = vect[i] * vect[j];
      };
   };
   row[1][0] = row[0][1];
   row[2][0] = row[0][2];
   row[2][1] = row[1][2];
};

/*!
\author Vladimir Florinski
\date 09/15/2022
\param[in] vect_l First vector \f$\mathbf{v}_1\f$
\param[in] vect_r Second vector \f$\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::Dyadic(const GeoVector& vect_l, const GeoVector& vect_r)
{
   for (auto i = 0; i < 3; i++) {
      for (auto j = 0; j < 3; j++) {
         row[i][j] = vect_l[i] * vect_r[j];
      };
   };
};

/*!
\author Vladimir Florinski
\date 09/09/2022
\param[in] vect_l First vector \f$\mathbf{v}_1\f$
\param[in] vect_r Second vector \f$\mathbf{v}_2\f$

\note mean can be subtracted by calling IncrementCov(Bmean, -Bmean);
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::IncrementCov(const GeoVector& vect_l, const GeoVector& vect_r)
{
   for (auto i = 0; i < 3; i++) {
      for (auto j = 0; j < 3; j++) {
         row[i][j] += vect_l[i] * vect_r[j];
      };
   };
};

/*!
\author Vladimir Florinski
\date 08/30/2022
\param[in] ez A vector in the new z-direction
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::AxisymmetricBasis(const GeoVector& ez)
{
   row[2] = UnitVec(ez);
   row[0] = GetSecondUnitVec(row[2]);
   row[1] = row[2] ^ row[0];
};

/*!
\author Vladimir Florinski
\date 08/31/2022
\param[in] basis A matrix whose rows are new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::ChangeToBasis(const GeoMatrix& basis)
{
   for (auto uvw = 0; uvw < 3; uvw++) row[uvw].ChangeToBasis(basis.VectorArray());
};

/*!
\author Juan G Alonso Guzman
\date 07/24/2025
\param[in] basis A matrix whose components are new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::ChangeToBasis(const GeoVector* basis)
{
// TODO write this function
};

/*!
\author Vladimir Florinski
\date 08/31/2022
\param[in] basis A matrix whose rows are new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::ChangeFromBasis(const GeoMatrix& basis)
{
   for (auto uvw = 0; uvw < 3; uvw++) row[uvw].ChangeFromBasis(basis.VectorArray());
};

/*!
\author Juan G Alonso Guzman
\date 07/24/2025
\param[in] basis A matrix whose components are new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::ChangeFromBasis(const GeoVector* basis)
{
   GeoMatrix mat_tmp;
   memcpy(mat_tmp.linear, linear, 9 * SZDBL);
   for (auto xyz = 0; xyz < 3; xyz++) {
      for (auto uvw = 0; uvw < 3; uvw++) row[xyz][uvw] = basis[xyz] * mat_tmp * basis[uvw];
   }
};

/*!
\author Vladimir Florinski
\date 08/11/2022
\param[in] i First index
\param[in] j Second index
\return A minor corresponding to (i,j)
*/
SPECTRUM_DEVICE_FUNC double GeoMatrix::Minor(int i, int j) const
{
   int i1, i2, j1, j2;
   
   i1 = (i + 1) % 3;
   i2 = (i + 2) % 3;
   j1 = (j + 1) % 3;
   j2 = (j + 2) % 3;

   return row[i1][j1] * row[i2][j2] - row[i1][j2] * row[i2][j1];
};

/*!
\author Vladimir Florinski
\date 08/11/2022
\return The inverse of this matrix
*/
SPECTRUM_DEVICE_FUNC GeoMatrix GeoMatrix::Inverse(void) const
{
   int i, j;
   GeoMatrix inverse;
   double det = Det();
   
   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
         inverse[i][j] = Minor(j, i) / det;
      };
   };
   return inverse;
};

/*!
\author Vladimir Florinski
\date 08/11/2022
\return Three eigenvalues as a vector
*/
SPECTRUM_DEVICE_FUNC GeoVector GeoMatrix::Eigenvalues(void) const
{
   double a, b, c, d;
   GeoVector eval;
   
   a = 1.0;
   b = -Trace();
   c = row[0][0] * row[1][1] + row[1][1] * row[2][2] + row[0][0] * row[2][2]
     - row[0][1] * row[1][0] - row[0][2] * row[2][0] - row[1][2] * row[2][1];
   d = -Det();

   CubicSolve(a, b, c, d, eval[0], eval[1], eval[2]);
   return eval;
};

/*!
\author Vladimir Florinski
\date 08/11/2022
\param[out] evec Normalized left eigenvectors, but arranged in rows
\return Three eigenvalues as a vector
*/
SPECTRUM_DEVICE_FUNC GeoVector GeoMatrix::Eigensystem(GeoMatrix& evec) const
{
   int i, k, ei, ei1, ei2;
   GeoVector eval;
   GeoMatrix mchar[3];

// Find eigenvalues
   eval = Eigenvalues();
   if (eval.Norm() < sp_tiny) {
      evec = gm_zeros;
      return gv_zeros;
   };

// Build characteristic matrices
   for (ei = 0; ei < 3; ei++) mchar[ei] = *this - eval[ei] * gm_unit;

// Apply the operator C(ei1)C(ei2) to the initial vector (1,1,1) to eliminate all components except that parallel to "ei". This could fail if one of the eigenvectors corresponding to ei1 or ei2 is parallel to the initial vector.
   for (ei = 0; ei < 3; ei++) {
      ei1 = (ei + 1) % 3;
      ei2 = (ei + 2) % 3;

      for (i = 0; i < 3; i++) {
         evec[ei][i] = 0.0;
         for (k = 0; k < 3; k++) {
            evec[ei][i] += mchar[ei1][i][k] * (mchar[ei2][k][0] + mchar[ei2][k][1] + mchar[ei2][k][2]);
         };
      };

// Normalize the eigenvector
      evec[ei].Normalize();
   };
   return eval;
};

/*!
\author Vladimir Florinski
\date 10/13/2025
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::CompleteSymmetric(void)
{
   row[1][0] = row[0][1];
   row[2][0] = row[0][2];
   row[2][1] = row[1][2];
};

/*!
\author Vladimir Florinski
\date 10/13/2025
*/
SPECTRUM_DEVICE_FUNC void GeoMatrix::CompleteAntiSymmetric(void)
{
   row[0][0] = 0.0;
   row[1][0] = -row[0][1];
   row[1][1] = 0.0;
   row[2][0] = -row[0][2];
   row[2][1] = -row[1][2];
   row[2][2] = 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on matrices
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/09/2022
\param[in] vect_l First vector \f$\mathbf{v}_1\f$
\param[in] vect_r Second vector \f$\mathbf{v}_2\f$
\return Dyadic matrix
*/
SPECTRUM_DEVICE_FUNC GeoMatrix Dyadic(const GeoVector& vect_l, const GeoVector& vect_r)
{
   GeoMatrix cov_matr;

   for (auto i = 0; i < 3; i++) {
      for (auto j = 0; j < 3; j++) {
         cov_matr[i][j] = vect_l[i] * vect_r[j];
      };
   };
   
   return cov_matr;
};

/*!
\author Vladimir Florinski
\date 08/02/2023
\param[in] vect Vector to build dyadic from
\return Dyadic matrix
*/
SPECTRUM_DEVICE_FUNC GeoMatrix Dyadic(const GeoVector& vect)
{
   GeoMatrix dyadic;
   dyadic.Dyadic(vect);
   return dyadic;
};

/*!
\author Vladimir Florinski
\date 11/03/2021
\param[in] matr_r right operand
\return Modified "ostream" object
*/
std::ostream& operator <<(std::ostream& os, const GeoMatrix& matr_r)
{
   os << matr_r[0] << "\n" << matr_r[1] << "\n" << matr_r[2];
   return os;
};

};
