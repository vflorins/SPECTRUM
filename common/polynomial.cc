/*!
\file polynomial.cc
\brief Defines common operation on polynomials
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifdef GEO_DEBUG
#include <iostream>
#include <iomanip>
#endif
#include "common/polynomial.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] v    Position vector
\param[in] mlex Three power law indices encoded as a single number
\return \f$x^l \cdot y^m \cdot z^n\f$
*/
SPECTRUM_DEVICE_FUNC double CoordPower3D(const GeoVector& v, uint32_t mlex)
{
   uint32_t i;
   double res = 1.0;

// To compute the exponent, all bits except those relevant to the desired vector component are set to zero by applying the mask.
   for (i = 0; i < (mlex & xbitfield) >> 16; i++) res *= v[0];
   for (i = 0; i < (mlex & ybitfield) >> 8 ; i++) res *= v[1];
   for (i = 0; i < (mlex & zbitfield)      ; i++) res *= v[2];

   return res;
};

/*!
\author Vladimir Florinski
\date 07/17/2019
\param[in] x    A number to be raised
\param[in] mlex Three power law indices encoded as a single number
\return \f$x^{l+m+n}\f$
*/
SPECTRUM_DEVICE_FUNC double CoordPower3D(double x, uint32_t mlex)
{
   int i, pl_tot;
   double res = 1.0;

   pl_tot = ((mlex & xbitfield) >> 16) + ((mlex & ybitfield) >> 8) + (mlex & zbitfield);
   for (i = 0; i < pl_tot; i++) res *= x;
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
SPECTRUM_DEVICE_FUNC double CoordPower3D(const GeoVector& v, int l, int m, int n)
{
   int i;
   double res = 1.0;

   for (i = 0; i < l; i++) res *= v[0];
   for (i = 0; i < m; i++) res *= v[1];
   for (i = 0; i < n; i++) res *= v[2];
   
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
SPECTRUM_DEVICE_FUNC double EvalPoly1D(int p, const double* c, double x)
{
   int i;
   double sum = c[p];

   for (i = p - 1; i >= 0; i--) sum = sum * x + c[i];
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
SPECTRUM_DEVICE_FUNC double EvalPoly2D(int p, const double* c, double x, double y)
{
   int i, j, idx = 0;
   double psum, sum = c[0];

// Precompute powers of y
   double powy[p + 1];
   powy[0] = 1.0;
   for (i = 1; i <= p; i++) powy[i] = powy[i - 1] * y;

// Loop on powers of xy
   for (j = 1; j <= p; j++) {
      psum = c[++idx];

// Loop on powers of y
      for (i = 1; i <= j; i++) psum = psum * x + c[++idx] * powy[i];
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
SPECTRUM_DEVICE_FUNC double EvalPoly3D(int p, const double* c, double x, double y, double z)
{
   int i, j, k, idx = 0;
   double psum, qsum, sum = c[0];

// Precompute powers of z
   double powz[p + 1];
   powz[0] = 1.0;
   for (i = 1; i <= p; i++) powz[i] = powz[i - 1] * z;

// Loop on powers of xyz
   for (j = 1; j <= p; j++) {
      psum = c[++idx];

// Loop on powers of yz
      for (i = 1; i <= j; i++) {
         qsum = c[++idx];

// Loop on powers of z
         for (k = 1; k <= i; k++) qsum = qsum * y + c[++idx] * powz[k];
         psum = psum * x + qsum;
      };
      sum += psum;
   };

   return sum;
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 05/10/2024
*/
void PrintMomentTables(void)
{
   Polynomial poly;
   std::cerr << std::setiosflags(std::ios::showbase);

   std::cerr << "Printing the binomial coefficients\n";   
   for (auto n = 0; n <= MONO_DEGREE_HIGH; n++) {
      for (auto k = 0; k <= n; k++) {
         std::cerr << std::dec << std::setw(3) << poly.binomial[n][k];
      };
      std::cerr << std::endl;
   };

   std::cerr << "\nPrinting the moment numbering table\n";   
   for (auto idx = 0; idx < poly_table_length; idx++) {
      std::cerr << std::dec << std::setw(3) << idx;
      std::cerr << std::hex << std::setw(12) << poly.moment_pl[idx] << std::endl;
   };

   std::cerr << "\nPrinting the moment lookup table\n";   
   for (auto k = 0; k <= MONO_DEGREE_HIGH; k++) {
      for (auto j = 0; j <= MONO_DEGREE_HIGH; j++) {
         for (auto i = 0; i <= MONO_DEGREE_HIGH; i++) {
            std::cerr << std::dec << std::setw(3) << poly.moment_lu[k][j][i];
         };
         std::cerr << std::endl;
      };
      std::cerr << std::endl;
   };
};

#endif

};
