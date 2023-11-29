/*!
\file vectors.hh
\brief Declares and defines a three-component vector class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_VECTORS_HH
#define SPECTRUM_VECTORS_HH

#include "definitions.hh"
#include "multi_index.hh"

namespace Spectrum {

//! Size of this class (should be 24)
#define SZGV sizeof(GeoVector)

//! The unit vector in the x direction
#define gv_nx GeoVector(1.0, 0.0, 0.0)

//! The unit vector in the y direction
#define gv_ny GeoVector(0.0, 1.0, 0.0)

//! The unit vector in the z direction
#define gv_nz GeoVector(0.0, 0.0, 1.0)

//! The unit vector in the r direction
#define gv_nr gv_nx

//! The unit vector in the theta direction
#define gv_nt gv_ny

//! The unit vector in the phi direction
#define gv_np gv_nz

//! A vector with zero components
#define gv_zeros GeoVector(0.0, 0.0, 0.0)

//! A vector with unit components
#define gv_ones GeoVector(1.0, 1.0, 1.0)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoVector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three component vector class
\author Vladimir Florinski

A class operating on a three-component vector. The following operations are provided: copying, reflection, normalization, rescaling, scalar product, vector product. Almost all functions are inlined.
*/
struct GeoVector {

//! Storage
   union {
      struct {
         double x, y, z;
      };
      double data[3];
   };

//! Default constructor
   SPECTRUM_DEVICE_FUNC GeoVector(void);

//! Constructor from components
   SPECTRUM_DEVICE_FUNC constexpr GeoVector(double x_in, double y_in, double z_in);

//! Constructor from an array
   SPECTRUM_DEVICE_FUNC GeoVector(const double* data_r);

//! Copy constructor
   SPECTRUM_DEVICE_FUNC GeoVector(const GeoVector& other);

//! Access to the data for reading
   SPECTRUM_DEVICE_FUNC const double* Data(void) const;

//! Access to the data for writing
   SPECTRUM_DEVICE_FUNC double* Data(void);

//! Access to components for reading
   SPECTRUM_DEVICE_FUNC const double& operator [](int i) const;

//! Access to components for writing
   SPECTRUM_DEVICE_FUNC double& operator [](int i);

//! Store the content of the vector into an array
   SPECTRUM_DEVICE_FUNC void Store(double* data_l) const;

//! Store the content of the vector into three separate components
   SPECTRUM_DEVICE_FUNC void Store(double& x, double& y, double& z) const;

//! Assignment operator from an array
   SPECTRUM_DEVICE_FUNC GeoVector& operator =(const double* data_r);

//! Assignment operator from another vector
   SPECTRUM_DEVICE_FUNC GeoVector& operator =(const GeoVector& vect_r);

//! Set all three components to the given value
   SPECTRUM_DEVICE_FUNC GeoVector& operator =(double val);

//! Convertion operator to MultiIndex
   SPECTRUM_DEVICE_FUNC operator MultiIndex(void) const;

//! Add another vector to this
   SPECTRUM_DEVICE_FUNC GeoVector& operator +=(const GeoVector& vect_r);

//! Add a multi-index to this
   SPECTRUM_DEVICE_FUNC GeoVector& operator +=(const MultiIndex& midx_r);

//! Subtract another vector from this
   SPECTRUM_DEVICE_FUNC GeoVector& operator -=(const GeoVector& vect_r);

//! Subtract a multi-index from this
   SPECTRUM_DEVICE_FUNC GeoVector& operator -=(const MultiIndex& midx_r);

//! Multiply this vector by a scalar
   SPECTRUM_DEVICE_FUNC GeoVector& operator *=(double sclr_r);

//! Divide this vector by a scalar
   SPECTRUM_DEVICE_FUNC GeoVector& operator /=(double sclr_r);

//! Computes the norm of this vector
   SPECTRUM_DEVICE_FUNC double Norm(void) const;

//! Computes the square of the norm of this vector
   SPECTRUM_DEVICE_FUNC double Norm2(void) const;

//! Makes this a unit vector
   SPECTRUM_DEVICE_FUNC GeoVector& Normalize(void);

//! Makes this a unit vector but saves the norm
   SPECTRUM_DEVICE_FUNC GeoVector& Normalize(double& norm);

//! Sum of the components
   SPECTRUM_DEVICE_FUNC double Sum(void) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Smallest component
   SPECTRUM_DEVICE_FUNC double Smallest(void) const;

//! Polar angle of a position vector
   SPECTRUM_DEVICE_FUNC double Theta(void);

//! Azimuthal angle of a position vector
   SPECTRUM_DEVICE_FUNC double Phi(void);

//! Converts a position vector from r,theta,phi to x,y,z
   SPECTRUM_DEVICE_FUNC void RTP_XYZ(void);

//! Converts a position vector from x,y,z to r,theta,phi
   SPECTRUM_DEVICE_FUNC void XYZ_RTP(void);

//! Converts a vector to spherical coordinates
   SPECTRUM_DEVICE_FUNC void ToSpherical(double sintheta, double costheta, double sinphi, double cosphi);

//! Converts a vector to Cartesian coordinates
   SPECTRUM_DEVICE_FUNC void ToCartesian(double sintheta, double costheta, double sinphi, double cosphi);

//! Rotate about a given unit vector
   SPECTRUM_DEVICE_FUNC void Rotate(const GeoVector& n, double sina, double cosa);

//! Rotate about an axis
   SPECTRUM_DEVICE_FUNC void Rotate(const GeoVector& axis, double angle);

//! Convert components from the standard basis to a different basis
   SPECTRUM_DEVICE_FUNC void ChangeToBasis(const GeoVector* basis);

//! Convert components of a variance vector from the standard basis to a different basis
   SPECTRUM_DEVICE_FUNC void ChangeToBasis2(const GeoVector* basis);

//! Convert components to the standard basis from a different basis
   SPECTRUM_DEVICE_FUNC void ChangeFromBasis(const GeoVector* basis);

//! Project a point on a sphere onto a tangent plane
   SPECTRUM_DEVICE_FUNC void ProjectToPlane(int projection, const GeoVector& normal, const GeoVector& north, double& xi, double& eta) const;

//! Check if a vector is inside a wedge
   SPECTRUM_DEVICE_FUNC bool InsideWedge(int n_verts, const GeoVector* verts, double tol) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Computes the reflected vector
   SPECTRUM_DEVICE_FUNC friend GeoVector operator -(const GeoVector& vect);

//! Add an array and a vector
   SPECTRUM_DEVICE_FUNC friend GeoVector operator +(const double* data_l, const GeoVector& vect_r);

//! Add a vector and an array
   SPECTRUM_DEVICE_FUNC friend GeoVector operator +(const GeoVector& vect_l, const double* data_r);

//! Add two vectors together
   SPECTRUM_DEVICE_FUNC friend GeoVector operator +(const GeoVector& vect_l, const GeoVector& vect_r);

//! Add a vector to scaled ones
   SPECTRUM_DEVICE_FUNC friend GeoVector operator +(double sclr_l, const GeoVector& vect_r);

//! Increment all components by an equal amount
   SPECTRUM_DEVICE_FUNC friend GeoVector operator +(const GeoVector& vect_l, double sclr_r);

//! Subtract a vector from an array
   SPECTRUM_DEVICE_FUNC friend GeoVector operator -(const double* data_l, const GeoVector& vect_r);

//! Subtract an array from a vector
   SPECTRUM_DEVICE_FUNC friend GeoVector operator -(const GeoVector& vect_l, const double* data_r);

//! Subtract one vector from another
   SPECTRUM_DEVICE_FUNC friend GeoVector operator -(const GeoVector& vect_l, const GeoVector& vect_r);

//! Subtract a vector from scaled ones
   SPECTRUM_DEVICE_FUNC friend GeoVector operator -(double sclr_l, const GeoVector& vect_r);

//! Decrement all components by an equal amount
   SPECTRUM_DEVICE_FUNC friend GeoVector operator -(const GeoVector& vect_l, double sclr_r);

//! Multiply a vector by a scalar from the right
   SPECTRUM_DEVICE_FUNC friend GeoVector operator *(const GeoVector& vect_l, double sclr_r);

//! Multiply a vector by a scalar from the left
   SPECTRUM_DEVICE_FUNC friend GeoVector operator *(double sclr_l, const GeoVector& vect_r);

//! Divide a vector by a scalar
   SPECTRUM_DEVICE_FUNC friend GeoVector operator /(const GeoVector& vect_l, double sclr_r);

//! Divide a vector by another vector component-wise
   SPECTRUM_DEVICE_FUNC friend GeoVector operator /(const GeoVector& vect_l, const GeoVector& vect_r);

//! Multiply a vector by a multi-index component-wise
   SPECTRUM_DEVICE_FUNC friend GeoVector operator *(const GeoVector& vect_l, const MultiIndex& midx_r);

//! Multiply a multi-index by a vector component-wise
   SPECTRUM_DEVICE_FUNC friend GeoVector operator *(const MultiIndex& midx_l, const GeoVector& vect_r);

//! Divide a vector by a multi-index component-wise
   SPECTRUM_DEVICE_FUNC friend GeoVector operator /(const GeoVector& vect_l, const MultiIndex& midx_r);

//! Return a scalar product of a vector and an array
   SPECTRUM_DEVICE_FUNC friend double operator *(const GeoVector& vect_l, const double* data_r);

//! Return a scalar product of an array and a vector
   SPECTRUM_DEVICE_FUNC friend double operator *(const double* data_l, const GeoVector& vect_r);

//! Return a scalar product of two vectors
   SPECTRUM_DEVICE_FUNC friend double operator *(const GeoVector& vect_l, const GeoVector& vect_r);

//! Compute a vector product of a vector and an array
   SPECTRUM_DEVICE_FUNC friend GeoVector operator ^(const GeoVector& vect_l, const double* data_r);

//! Compute a vector product of an array and a vector
   SPECTRUM_DEVICE_FUNC friend GeoVector operator ^(const double* data_l, const GeoVector& vect_r);

//! Compute a vector product of two vectors
   SPECTRUM_DEVICE_FUNC friend GeoVector operator ^(const GeoVector& vect_l, const GeoVector& vect_r);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoVector inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/19/2023
*/
SPECTRUM_DEVICE_FUNC inline GeoVector::GeoVector(void)
{
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] x First component
\param[in] y Second component
\param[in] z Third component
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector::GeoVector(double x_in, double y_in, double z_in)
                                               : x(x_in),
                                                 y(y_in),
                                                 z(z_in)
{
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] data_r Array to initialize the vector from
*/
SPECTRUM_DEVICE_FUNC inline GeoVector::GeoVector(const double* data_r)
{
   std::memcpy(data, data_r, 3 * SZDBL);
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] other Vector to create a copy of
*/
SPECTRUM_DEVICE_FUNC inline GeoVector::GeoVector(const GeoVector& other)
{
   std::memcpy(data, &other, 3 * SZDBL);;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return "data" array
*/
SPECTRUM_DEVICE_FUNC inline const double* GeoVector::Data(void) const
{
   return data;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return "data" array
*/
SPECTRUM_DEVICE_FUNC inline double* GeoVector::Data(void)
{
   return data;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] i The desired component
\return \f$v_i\f$
*/
SPECTRUM_DEVICE_FUNC inline const double& GeoVector::operator [](int i) const
{
   return data[i];
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] i The desired component
\return \f$v_i\f$
*/
SPECTRUM_DEVICE_FUNC inline double& GeoVector::operator [](int i)
{
   return data[i];
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[out] data_l Array to store the vector in
*/
SPECTRUM_DEVICE_FUNC inline void GeoVector::Store(double* data_l) const
{
   std::memcpy(data_l, data, 3 * SZDBL);
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[out] x First component
\param[out] y Second component
\param[out] z Third component
*/
SPECTRUM_DEVICE_FUNC inline void GeoVector::Store(double& x, double& y, double& z) const
{
   x = data[0];
   y = data[1];
   z = data[2];
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] data_r Array to copy into this vector
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator =(const double* data_r)
{
   std::memcpy(data, data_r, 3 * SZDBL);
   return *this;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] vect_r A vector that will be copied into this vector
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator =(const GeoVector& vect_r)
{
   if(this != &vect_r) std::memcpy(data, vect_r.data, 3 * SZDBL);
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] val A number to be assigned to all three components
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator =(double val)
{
   data[0] = data[1] = data[2] = val;
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/15/2021
\return A MultiIndex whose components are truncated copies of the vector
*/
SPECTRUM_DEVICE_FUNC inline GeoVector::operator MultiIndex(void) const
{
   return MultiIndex(data[0], data[1], data[2]);
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] vect_r right operand \f$\mathbf{v}_1\f$
\return \f$\mathbf{v}+\mathbf{v}_1\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator +=(const GeoVector& vect_r)
{
   data[0] += vect_r.data[0];
   data[1] += vect_r.data[1];
   data[2] += vect_r.data[2];
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] midx_r right operand
\return Sum of this vector and a multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator +=(const MultiIndex& midx_r)
{
   data[0] += midx_r.i;
   data[1] += midx_r.j;
   data[2] += midx_r.k;
   return *this;
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] data_r right operand \f$\mathbf{v}_1\f$
\return \f$\mathbf{v}-\mathbf{v}_1\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator -=(const GeoVector& vect_r)
{
   data[0] -= vect_r.data[0];
   data[1] -= vect_r.data[1];
   data[2] -= vect_r.data[2];
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] midx_r right operand
\return Difference between this vector and a multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator -=(const MultiIndex& midx_r)
{
   data[0] -= midx_r.i;
   data[1] -= midx_r.j;
   data[2] -= midx_r.k;
   return *this;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] sclr_r right operand \f$a\f$
\return \f$a\mathbf{v}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator *=(double sclr_r)
{
   data[0] *= sclr_r;
   data[1] *= sclr_r;
   data[2] *= sclr_r;
   return *this;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] sclr_r right operand \f$a\f$
\return \f$a^{-1}\mathbf{v}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator /=(double sclr_r)
{
   data[0] /= sclr_r;
   data[1] /= sclr_r;
   data[2] /= sclr_r;
   return *this;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\return \f$|v|\f$
*/
SPECTRUM_DEVICE_FUNC inline double GeoVector::Norm(void) const
{
   return sqrt(fmax(Sqr(data[0]) + Sqr(data[1]) + Sqr(data[2]), 0.0));
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\return \f$|v|^2\f$
*/
SPECTRUM_DEVICE_FUNC inline double GeoVector::Norm2(void) const
{
   return Sqr(data[0]) + Sqr(data[1]) + Sqr(data[2]);
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::Normalize(void)
{
   double norm = Norm();
   data[0] /= norm;
   data[1] /= norm;
   data[2] /= norm;
   return *this;
};

/*!
\author Vladimir Florinski
\date 12/02/2020
\param[out] norm \f$|v|\f$
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::Normalize(double& norm)
{
   norm = Norm();
   data[0] /= norm;
   data[1] /= norm;
   data[2] /= norm;
   return *this;
};

/*!
\author Vladimir Florinski
\date 05/14/2021
\return Sum of components
*/
SPECTRUM_DEVICE_FUNC inline double GeoVector::Sum(void) const
{
   return data[0] + data[1] + data[2];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on vectors
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Convert input to a unit vector
SPECTRUM_DEVICE_FUNC GeoVector UnitVec(const GeoVector& vect);

//! Compute a unit vector bisecting the angle formed by two unit vectors
SPECTRUM_DEVICE_FUNC GeoVector Bisect(const GeoVector& vect_l, const GeoVector& vect_r);

//! Area of a flat triangle defined with three vectors
SPECTRUM_DEVICE_FUNC double TriangleArea(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r);

//! Unit normal vector to a plane defined with three vectors
SPECTRUM_DEVICE_FUNC GeoVector PlaneNormal(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r);

//! Geometric center of a flat triangle defined with three vectors
SPECTRUM_DEVICE_FUNC GeoVector TriangleCenter(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r);

//! Height of a cone defined with three vectors (= circumcenter of the triangle)
SPECTRUM_DEVICE_FUNC GeoVector CircumCenter(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r);

//! Divide the angle between two vectors in a given proportion
SPECTRUM_DEVICE_FUNC GeoVector DivideAngle(const GeoVector& vect_l, const GeoVector& vect_r, double frac);

//! Barycentric coordinates of a point inside a flat triangle
SPECTRUM_DEVICE_FUNC GeoVector BarycentricCoords(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3, const GeoVector& vect);

//! Moments of a flat triangle
SPECTRUM_DEVICE_FUNC void TriangleMoments(const GeoVector v[3], double& C0, double C1[3], double C2[3][3], double C3[3][3][3], double C4[3][3][3][3]);

//! Cosine of the vertex angle at the middle operand
SPECTRUM_DEVICE_FUNC double VertexAngle(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r);

//! Area of a spherical triangle
SPECTRUM_DEVICE_FUNC double SphTriArea(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3);

//! Area of a spherical quadrilateral
SPECTRUM_DEVICE_FUNC double SphQuadArea(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3, const GeoVector& vect4);

//! Geometric center of a spherical triangle
SPECTRUM_DEVICE_FUNC GeoVector SphTriCenter(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3);

//! Calculates a point of intersection of two great circles
SPECTRUM_DEVICE_FUNC GeoVector GreatCircleInt(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3, const GeoVector& vect4);

//! Find a unit vector normal to a given vector
SPECTRUM_DEVICE_FUNC GeoVector GetSecondUnitVec(const GeoVector& first);

//! Find a unit vector normal to a given vector in the plane given by the second vector
SPECTRUM_DEVICE_FUNC GeoVector GetSecondUnitVec(const GeoVector& vect_l, const GeoVector& vect_r);

//! Stream insertion operator (host only)
std::ostream& operator <<(std::ostream& os, const GeoVector& vect_r);

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Three unit vectors in an array
#ifdef __CUDA_ARCH__
__device__ constexpr GeoVector cart_unit_vec[] = {gv_nx, gv_ny, gv_nz};
#else
constexpr GeoVector cart_unit_vec[] = {gv_nx, gv_ny, gv_nz};
#endif

};

#endif
