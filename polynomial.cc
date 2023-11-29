/*!
\file polynomial.cc
\brief Defines common operation on polynomials
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "polynomial.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] v  Position vector
\param[in] pl Three power law indices encoded as a single number
\return \f$x^l \cdot y^m \cdot z^n\f$
*/
double CoordPower3D(const GeoVector& v, uint32_t pl)
{
   uint32_t i;
   double res = 1.0;

// To compute the exponent, all bits except those relevant to the desired vector component are set to zero by applying the mask.
   for(i = 0; i < (pl & xbitfield) >> 16; i++) res *= v[0];
   for(i = 0; i < (pl & ybitfield) >> 8 ; i++) res *= v[1];
   for(i = 0; i < (pl & zbitfield)      ; i++) res *= v[2];

   return res;
};

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] r  Distance from the center of the sphere
\param[in] pl Three power law indices encoded as a single number
\return \f$r^{l+m+n}\f$
*/
double CoordPower3D(double r, uint32_t pl)
{
   int i, pl_tot;
   double res = 1.0;

   pl_tot = ((pl & xbitfield) >> 16) + ((pl & ybitfield) >> 8) + (pl & zbitfield);
   for(i = 0; i < pl_tot; i++) res *= r;
   return res;
};

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] v Position vector
\param[in] l Power law for x
\param[in] m Power law for y
\param[in] n Power law for z
\return \f$x^l \cdot y^m \cdot z^n\f$
*/
double CoordPower3D(const GeoVector& v, int l, int m, int n)
{
   int i;
   double res = 1.0;

   for(i = 0; i < l; i++) res *= v[0];
   for(i = 0; i < m; i++) res *= v[1];
   for(i = 0; i < n; i++) res *= v[2];
   
   return res;
};

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] p Degree of the polynomial
\param[in] c Coefficients of the polynomial
\param[in] x Argument
\return \f$P(x)=c_0+c_1x+c_2x^2+...\f$
*/
double EvalPoly1D(int p, const double* c, double x)
{
   int i;
   double sum = c[p];

   for(i = p - 1; i >= 0; i--) sum = sum * x + c[i];
   return sum;
};

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] p Degree of the polynomial
\param[in] c Coefficients of the polynomial
\param[in] x First argument
\param[in] y Second argument
\return \f$P(x,y)=c_0+c_1x+c_2y+c_3x^2+c_4xy+c_5y^2+...\f$
*/
double EvalPoly2D(int p, const double* c, double x, double y)
{
   int i, j, idx = 0;
   double psum, sum = c[0];

// Precompute powers of y
   double powy[p + 1];
   powy[0] = 1.0;
   for(i = 1; i <= p; i++) powy[i] = powy[i - 1] * y;

// Loop on powers of xy
   for(j = 1; j <= p; j++) {
      psum = c[++idx];

// Loop on powers of y
      for(i = 1; i <= j; i++) psum = psum * x + c[++idx] * powy[i];
      sum += psum;
   };

   return sum;
};

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] p Degree of the polynomial
\param[in] c Coefficients of the polynomial
\param[in] x First argument
\param[in] y Second argument
\param[in] z Third argument
\return \f$P(x,y,z)=c_0+c_1x+c_2y+c_3z+c_4x^2+c_5xy+c_6xz+c_7y^2+c_8yz+c_9z^2+...\f$
*/
double EvalPoly3D(int p, const double* c, double x, double y, double z)
{
   int i, j, k, idx = 0;
   double psum, qsum, sum = c[0];

// Precompute powers of z
   double powz[p + 1];
   powz[0] = 1.0;
   for(i = 1; i <= p; i++) powz[i] = powz[i - 1] * z;

// Loop on powers of xyz
   for(j = 1; j <= p; j++) {
      psum = c[++idx];

// Loop on powers of yz
      for(i = 1; i <= j; i++) {
         qsum = c[++idx];

// Loop on powers of z
         for(k = 1; k <= i; k++) qsum = qsum * y + c[++idx] * powz[k];
         psum = psum * x + qsum;
      };
      sum += psum;
   };

   return sum;
};

};
