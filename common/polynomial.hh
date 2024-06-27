/*!
\file polynomial.hh
\brief Declares data and common operation on polynomials
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_POLYNOMIAL_HH
#define SPECTRUM_POLYNOMIAL_HH

#include <cstdint>
#include "common/vectors.hh"

namespace Spectrum {

//! The highest monomial degree for which bitwise encoding is available
#define MONO_DEGREE_HIGH 5

//! Length of tables
constexpr int poly_table_length = (MONO_DEGREE_HIGH + 1) * (MONO_DEGREE_HIGH + 2) * (MONO_DEGREE_HIGH + 3) / 6;

//! Bitwise mask for x
constexpr uint32_t xbitfield = 0xFF << 16;

//! Bitwise mask for y
constexpr uint32_t ybitfield = 0xFF << 8;

//! Bitwise mask for z
constexpr uint32_t zbitfield = 0xFF;

/*!
\brief A helper class to simplify working with polynomials
\author Vladimir Florinski
*/
struct Polynomial
{
//! Assembled Pascal triange
   static int binomial[MONO_DEGREE_HIGH + 1][MONO_DEGREE_HIGH + 1];

//! Bitwise encoding of the monomial indices using reverse colexicographical ordering
   static uint32_t moment_pl[poly_table_length];

//! Inverse moment lookup table
   static int moment_lu[MONO_DEGREE_HIGH + 1][MONO_DEGREE_HIGH + 1][MONO_DEGREE_HIGH + 1];

//! Default constructor
   SPECTRUM_DEVICE_FUNC Polynomial(void);

//! Copy constructor is not needed because the class is all constants
   SPECTRUM_DEVICE_FUNC Polynomial(const Polynomial& other) = delete;

//! Calculate the polynomial access tables
   SPECTRUM_DEVICE_FUNC static constexpr void Setup(void);

//! Returns the number of degrees of freeedom (coefficients) in a polynomial
   SPECTRUM_DEVICE_FUNC int DOF(int dim, int poly_deg) const;
};

/*!
\author Vladimir Florinski
\date 05/09/2024
*/
SPECTRUM_DEVICE_FUNC inline Polynomial::Polynomial(void)
{
   Setup();
};

/*!
\author Vladimir Florinski
\date 05/09/2024
\param[in] dim      The number of dimensions
\param[in] poly_deg The polynomial degree
\return The number of coefficients in the polynomial
*/
SPECTRUM_DEVICE_FUNC inline int Polynomial::DOF(int dim, int poly_deg) const
{
   if(dim == 1) return poly_deg + 1;
   else if(dim == 2) return (poly_deg + 1) * (poly_deg + 2) / 2;
   else if(dim == 3) return (poly_deg + 1) * (poly_deg + 2) * (poly_deg + 3) / 6;
};

/*!
\author Vladimir Florinski
\date 05/09/2024
*/
SPECTRUM_DEVICE_FUNC inline constexpr void Polynomial::Setup(void)
{
// Compute the binomial coefficients
   for(auto n = 0; n <= MONO_DEGREE_HIGH; n++) {
      binomial[n][0] = 1;
      binomial[n][n] = 1;
      for(auto k = 1; k < n; k++) {
         binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
      };
   };

// Build the moment numbering table
   int midx = 0;
   for(auto pow_hi = 0; pow_hi <= MONO_DEGREE_HIGH; pow_hi++) {
      for(auto i = pow_hi; i >= 0; i--) {
         for(auto j = pow_hi - i ; j >= 0; j--) {
            moment_pl[midx++] = i * (1 << 16) + j * (1 << 8) + pow_hi - i - j;
         };
      };
   };

// Build the moment lookup table
   uint32_t mlex = 0x0;
   for(auto k = 0; k <= MONO_DEGREE_HIGH; k++) {
      for(auto j = 0; j <= MONO_DEGREE_HIGH; j++) {
         for(auto i = 0; i <= MONO_DEGREE_HIGH; i++) {
            moment_lu[k][j][i] = InList(poly_table_length, moment_pl, mlex);
            mlex += 1 << 16;
         };
         LOWER_BITS(mlex, xbitfield);
         mlex += 1 << 8;
      };
      LOWER_BITS(mlex, ybitfield);
      mlex += 1;
   };
};

//! Returns a product of powers of x, y, and z using bitfield input
SPECTRUM_DEVICE_FUNC double CoordPower3D(const GeoVector& v, uint32_t mlex);

//! Returns a power of r
SPECTRUM_DEVICE_FUNC double CoordPower3D(double x, uint32_t mlex);

//! Returns a product of powers of x, y, and z
SPECTRUM_DEVICE_FUNC double CoordPower3D(const GeoVector& v, int l, int m, int n);

//! Evaluates a polynomial of one argument
SPECTRUM_DEVICE_FUNC double EvalPoly1D(int p, const double* c, double x);

//! Evaluates a polynomial of two arguments
SPECTRUM_DEVICE_FUNC double EvalPoly2D(int p, const double* c, double x, double y);

//! Evaluates a polynomial of three arguments
SPECTRUM_DEVICE_FUNC double EvalPoly3D(int p, const double* c, double x, double y, double z);

#ifdef GEO_DEBUG
//! Prints the moments and the lookup table
void PrintMomentTables(void);
#endif

};

#endif
