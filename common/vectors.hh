/*!
\file vectors.hh
\brief Declares and defines a three-component vector class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_VECTORS_HH
#define SPECTRUM_VECTORS_HH

#include "common/multi_index.hh"

namespace Spectrum {

//! Size of this class (should be 24)
#define SZGV sizeof(GeoVector)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoVector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three component vector class
\author Vladimir Florinski
\author Juan G Alonso Guzman
*/
struct GeoVector : public SimpleArray<double, 3>
{
   using SimpleArray::operator=;
   using SimpleArray::operator+=;
   using SimpleArray::operator-=;
   using SimpleArray::operator*=;
   using SimpleArray::operator/=;

//! Default constructor
   SPECTRUM_DEVICE_FUNC constexpr GeoVector(void) {};

//! Constructor from a single value
   SPECTRUM_DEVICE_FUNC explicit constexpr GeoVector(double val);

//! Constructor from an array
   SPECTRUM_DEVICE_FUNC explicit constexpr GeoVector(const double* other);

//! Constructor from components
   SPECTRUM_DEVICE_FUNC constexpr GeoVector(double x_in, double y_in, double z_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC constexpr GeoVector(const SimpleArray<double, 3>& other);

//! Constructor from a multi-index
   SPECTRUM_DEVICE_FUNC constexpr GeoVector(const MultiIndex& other);

//! Store the content of the vector into three separate components
   SPECTRUM_DEVICE_FUNC void Store(double& x_out, double& y_out, double& z_out) const;

//! Conversion operator to MultiIndex
   SPECTRUM_DEVICE_FUNC operator MultiIndex(void) const;

//! Computes the norm of this vector
   SPECTRUM_DEVICE_FUNC double Norm(void) const;

//! Makes this a unit vector
   SPECTRUM_DEVICE_FUNC GeoVector& Normalize(void);

//! Makes this a unit vector but saves the norm
   SPECTRUM_DEVICE_FUNC GeoVector& Normalize(double& norm);

//! Vector-multiply this by another vector
   SPECTRUM_DEVICE_FUNC GeoVector& operator ^=(const GeoVector& other);

//! Add a multi-index to this
   SPECTRUM_DEVICE_FUNC GeoVector& operator +=(const MultiIndex& other);

//! Subtract a multi-index from this
   SPECTRUM_DEVICE_FUNC GeoVector& operator -=(const MultiIndex& other);

//! Multiply component-wise by a multi-index
   SPECTRUM_DEVICE_FUNC GeoVector& operator *=(const MultiIndex& other);

//! Divide component-wise by a multi-index
   SPECTRUM_DEVICE_FUNC GeoVector& operator /=(const MultiIndex& other);

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Rotate about a given unit vector
   SPECTRUM_DEVICE_FUNC void Rotate(const GeoVector& n, double sina, double cosa);

//! Rotate about an axis
   SPECTRUM_DEVICE_FUNC void Rotate(const GeoVector& axis, double angle);

//! Project onto normal space of a vector
   SPECTRUM_DEVICE_FUNC void SubtractParallel(const GeoVector& axis);

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
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeoVector inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/10/2024
\param[in] val Value to be asigned to each component
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector::GeoVector(double val)
                                               : SimpleArray(val)
{
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] other Array to initialize the vector from
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector::GeoVector(const double* other)
{
   std::memcpy(data, other, 3 * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] x_in First component
\param[in] y_in Second component
\param[in] z_in Third component
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector::GeoVector(double x_in, double y_in, double z_in)
{
   x = x_in;
   y = y_in;
   z = z_in;
};

/*!
\author Vladimir Florinski
\date 03/10/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector::GeoVector(const SimpleArray<double, 3>& other)
{
   memcpy(data, other.data, 3 * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/10/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector::GeoVector(const MultiIndex& other)
{
   x = other.i;
   y = other.j;
   z = other.k;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[out] x First component
\param[out] y Second component
\param[out] z Third component
*/
SPECTRUM_DEVICE_FUNC inline void GeoVector::Store(double& x_out, double& y_out, double& z_out) const
{
   x_out = x;
   y_out = y;
   z_out = z;
};

/*!
\author Vladimir Florinski
\date 06/15/2021
\return A MultiIndex whose components are truncated copies of the vector
*/
SPECTRUM_DEVICE_FUNC inline GeoVector::operator MultiIndex(void) const
{
   return MultiIndex(x, y, z);
};

/*!
\author Vladimir Florinski
\date 04/24/2024
\return \f$|v|\f$
*/
SPECTRUM_DEVICE_FUNC inline double GeoVector::Norm(void) const
{
   return sqrt(Norm2());
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::Normalize(void)
{
   double norm = Norm();
   x /= norm;
   y /= norm;
   z /= norm;
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
   x /= norm;
   y /= norm;
   z /= norm;
   return *this;
};

/*!
\author Vladimir Florinski
\date 04/25/2024
\param[in] other Right operand \f$\mathbf{v}_1\f$
\return \f$\mathbf{v}\times\mathbf{v}_1\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator ^=(const GeoVector& other)
{
   GeoVector vect_tmp(*this);
   x = vect_tmp.y * other.z - vect_tmp.z * other.y;
   y = vect_tmp.z * other.x - vect_tmp.x * other.z;
   z = vect_tmp.x * other.y - vect_tmp.y * other.x;
   return *this;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/16/2024
\param[in] other Right operand (multi-index)
\return The result of a summation with a multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator +=(const MultiIndex& other)
{
   x += other.x;
   y += other.y;
   z += other.z;
   return *this;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/16/2024
\param[in] other Right operand (multi-index)
\return The result of a subtraction of a multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator -=(const MultiIndex& other)
{
   x -= other.x;
   y -= other.y;
   z -= other.z;
   return *this;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/16/2024
\param[in] other Right operand (multi-index)
\return The result of a component-wise multiplication by a multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator *=(const MultiIndex& other)
{
   x *= other.x;
   y *= other.y;
   z *= other.z;
   return *this;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/16/2024
\param[in] other Right operand (multi-index)
\return The result of a component-wise division by a multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector& GeoVector::operator /=(const MultiIndex& other)
{
   x /= other.x;
   y /= other.y;
   z /= other.z;
   return *this;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on vectors
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Computes the reflected vector
\author Vladimir Florinski
\date 04/25/2024
\param[in] vect Vector to reflect \f$\mathbf{v}\f$
\return \f$-\mathbf{v}\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator -(const GeoVector& vect)
{
   GeoVector vect_tmp(vect);
   vect_tmp.Negate();
   return vect_tmp;
};

/*!
\brief Return a scalar product of two vectors
\author Vladimir Florinski
\date 04/25/2024
\param[in] sarr_l Left operand \f$\mathbf{v}_1\f$
\param[in] sarr_r Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1\cdot\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline double operator *(const GeoVector& vect_l, const GeoVector& vect_r)
{
   return vect_l.ScalarProd(vect_r);
};

/*!
\brief Compute a vector product of two vectors
\author Vladimir Florinski
\date 04/25/2024
\param[in] vect_l Left operand \f$\mathbf{v}_1\f$
\param[in] vect_r Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1\times\mathbf{v}_2\f$
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator ^(const GeoVector& vect_l, const GeoVector& vect_r)
{
   GeoVector vect_tmp(vect_l);
   vect_tmp ^= vect_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise addition of a vector and a multi-index
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/07/2024
\param[in] vect_l Left operand (vector)
\param[in] midx_r Right operand (multi-index)
\return Addition of vector and multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator +(const GeoVector& vect_l, const MultiIndex& midx_r)
{
   GeoVector vect_tmp(vect_l);
   vect_tmp += midx_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise addition of a multi-index and a vector
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/07/2024
\param[in] midx_l Left operand (multi-index)
\param[in] vect_r Right operand (vector)
\return Addition of vector and multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator +(const MultiIndex& midx_l, const GeoVector& vect_r)
{
   GeoVector vect_tmp(midx_l);
   vect_tmp += vect_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise subtraction of a multi-index from a vector
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/07/2024
\param[in] vect_l Left operand (vector)
\param[in] midx_r Right operand (multi-index)
\return Subtraction of vector and multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator -(const GeoVector& vect_l, const MultiIndex& midx_r)
{
   GeoVector vect_tmp(vect_l);
   vect_tmp -= midx_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise subtraction of a vector from a multi-index
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/07/2024
\param[in] midx_l Left operand (multi-index)
\param[in] vect_r Right operand (vector)
\return Subtraction of vector and multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator -(const MultiIndex& midx_l, const GeoVector& vect_r)
{
   GeoVector vect_tmp(midx_l);
   vect_tmp -= vect_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise multiplication of a vector by a multi-index
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/07/2024
\param[in] vect_l Left operand (vector)
\param[in] midx_r Right operand (multi-index)
\return Vector multiplied by the multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator *(const GeoVector& vect_l, const MultiIndex& midx_r)
{
   GeoVector vect_tmp(vect_l);
   vect_tmp *= midx_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise multiplication of a multi-index by a vector
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/07/2024
\param[in] midx_l Left operand (multi-index)
\param[in] vect_r Right operand (vector)
\return Vector multiplied by the multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator *(const MultiIndex& midx_l, const GeoVector& vect_r)
{
   GeoVector vect_tmp(midx_l);
   vect_tmp *= vect_r;
   return vect_tmp;
};

/*!
\brief Compute a component-wise division of a vector by a multi-index
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 04/26/2024
\param[in] vect_l Left operand (vector)
\param[in] midx_r Right operand (multi-index)
\return Vector divided by the multi-index
*/
SPECTRUM_DEVICE_FUNC inline GeoVector operator /(const GeoVector& vect_l, const MultiIndex& midx_r)
{
   GeoVector vect_tmp(vect_l);
   vect_tmp /= midx_r;
   return vect_tmp;
};

/*!
\brief Convert input to a unit vector
\author Vladimir Florinski
\date 07/10/2019
\param[in] vect Vector to rescale \f$\mathbf{v}\f$
\return Normalized vector
*/
SPECTRUM_DEVICE_FUNC inline GeoVector UnitVec(const GeoVector& vect)
{
   GeoVector vect_tmp(vect);
   return vect_tmp.Normalize();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

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

//! Length of a circular arc
SPECTRUM_DEVICE_FUNC double CircArcLength(const GeoVector& vect1, const GeoVector& vect2);

//! Geometric center of a circular arc
SPECTRUM_DEVICE_FUNC GeoVector CircArcCenter(const GeoVector& vect1, const GeoVector& vect2);

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

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Unit vector in the x direction
SPECTRUM_CONSTEXPR GeoVector gv_nx = {1.0, 0.0, 0.0};

//! Unit vector in the y direction
SPECTRUM_CONSTEXPR GeoVector gv_ny = {0.0, 1.0, 0.0};

//! Unit vector in the z direction
SPECTRUM_CONSTEXPR GeoVector gv_nz = {0.0, 0.0, 1.0};

//! A vector with zero components
SPECTRUM_CONSTEXPR GeoVector gv_zeros = {0.0, 0.0, 0.0};

//! A vector with unit components
SPECTRUM_CONSTEXPR GeoVector gv_ones = {1.0, 1.0, 1.0};

//! Three unit vectors in an array
SPECTRUM_CONSTEXPR GeoVector cart_unit_vec[3] = {gv_nx, gv_ny, gv_nz};

};

#endif
