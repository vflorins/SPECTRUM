/*!
\file coordinates.hh
\brief Declares some routines facilitating conversion between orthogonal coordinate systems
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COORDINATES_HH
#define SPECTRUM_COORDINATES_HH

#include "common/vectors.hh"

namespace Spectrum {

//! Coordinate system - defines curvilinear coordinates
enum class CoordinateSystem {
   Cartesian,
   Cylindrical,
   Spherical,
   Custom
};

/*!
\brief A collection of routines to convert between curvilinear and Cartesian coordinate systems 
\author Vladimir Florinski

The class only includes static routines and so can be used without an object.
*/
template <CoordinateSystem coord_system>
class Metric
{
public:

//! Convert a position vector from Cartesian to curvilinear
   SPECTRUM_DEVICE_FUNC static GeoVector PosToCurv(const GeoVector& pos) {return pos;};

//! Convert a position vector from Cartesian to curvilinear in place
   SPECTRUM_DEVICE_FUNC static void PosToCurv(GeoVector& pos) {};

//! Convert a position vector from curvilinear to Cartesian 
   SPECTRUM_DEVICE_FUNC static GeoVector PosToCart(const GeoVector& pos) {return pos;};

//! Convert a position vector from curvilinear to Cartesian in place 
   SPECTRUM_DEVICE_FUNC static void PosToCart(GeoVector& pos) {};

//! Convert the components of a vector from Cartesian to curvilinear
   SPECTRUM_DEVICE_FUNC static GeoVector CompToCurv(const GeoVector& pos, const GeoVector& comp) {return comp;};

//! Convert the components of a vector from curvilinear to Cartesian 
   SPECTRUM_DEVICE_FUNC static GeoVector CompToCart(const GeoVector& pos, const GeoVector& comp) {return comp;};

//! Scale factors 
   SPECTRUM_DEVICE_FUNC static GeoVector ScaleFactors(const GeoVector& pos) {return gv_ones;};

//! Jacobian 
   SPECTRUM_DEVICE_FUNC static double Jacobian(const GeoVector& pos) {return 1.0;};
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Cartesian coordinates
\return Position in Cylindrical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cylindrical>::PosToCurv(const GeoVector& pos)
{
   double phi = atan2(pos.y, pos.x);
   if (phi < 0.0) phi += M_2PI;
   return GeoVector(sqrt(Sqr(pos.x) + Sqr(pos.y)), phi, pos.z);
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Cartesian coordinates
\return Position in Spherical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Spherical>::PosToCurv(const GeoVector& pos)
{
   double r = pos.Norm();
   double phi = atan2(pos.y, pos.x);
   if (phi < 0.0) phi += M_2PI;
   return GeoVector(r, acos(pos.z / r), phi);
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in,out] pos Position vector to be converted
*/
template <>
SPECTRUM_DEVICE_FUNC inline void Metric<CoordinateSystem::Cylindrical>::PosToCurv(GeoVector& pos)
{
   double phi = atan2(pos.y, pos.x);
   if (phi < 0.0) phi += M_2PI;
   pos[0] = sqrt(Sqr(pos.x) + Sqr(pos.y));
   pos[1] = phi;
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in,out] pos Position vector to be converted
*/
template <>
SPECTRUM_DEVICE_FUNC inline void Metric<CoordinateSystem::Spherical>::PosToCurv(GeoVector& pos)
{
   double r = pos.Norm();
   double phi = atan2(pos.y, pos.x);
   if (phi < 0.0) phi += M_2PI;
   pos[0] = r;
   pos[1] = acos(pos.z / r);
   pos[2] = phi;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Cylindrical coordinates
\return Position in Cartesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cylindrical>::PosToCart(const GeoVector& pos)
{
   return GeoVector(pos[0] * cos(pos[1]), pos[0] * sin(pos[1]), pos[2]);
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Spherical coordinates
\return Position in Cartesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Spherical>::PosToCart(const GeoVector& pos)
{
   double r = pos.Norm();
   double s = r * sin(pos[1]);
   return GeoVector(s * cos(pos[2]), s * sin(pos[2]), r * cos(pos[1]));
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in,out] pos Position vector to be converted
*/
template <>
SPECTRUM_DEVICE_FUNC inline void Metric<CoordinateSystem::Cylindrical>::PosToCart(GeoVector& pos)
{
   double s = pos[0];
   pos.x = s * cos(pos[1]);
   pos.y = s * sin(pos[1]);
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in,out] pos Position vector to be converted
*/
template <>
SPECTRUM_DEVICE_FUNC inline void Metric<CoordinateSystem::Spherical>::PosToCart(GeoVector& pos)
{
   double r = pos.Norm();
   double s = r * sin(pos[1]);
   double z = r * cos(pos[1]);
   pos.x = s * cos(pos[2]);
   pos.y = s * sin(pos[2]);
   pos.z = z;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Cylindrical coordinates
\param[in] comp Components of a vector in Cartesian coordinates
\return Components of the vector in Cylindrical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cylindrical>::CompToCurv(const GeoVector& pos, const GeoVector& comp)
{
   double sinphi = sin(pos[1]);
   double cosphi = cos(pos[1]);
   return GeoVector(comp.x * cosphi + comp.y * sinphi, -comp.x * sinphi + comp.y * cosphi, comp.z);
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Spherical coordinates
\param[in] comp Components of a vector in Cartesian coordinates
\return Components of the vector in Spherical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Spherical>::CompToCurv(const GeoVector& pos, const GeoVector& comp)
{
   double sintheta = sin(pos[1]);
   double costheta = cos(pos[1]);
   double sinphi = sin(pos[2]);
   double cosphi = cos(pos[2]);
   double vxy = comp.x * cosphi + comp.y * sinphi;
   return GeoVector(vxy * sintheta + comp.z * costheta, vxy * costheta - comp.z * sintheta, -comp.x * sinphi + comp.y * cosphi);
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Cylindrical coordinates
\param[in] comp Components of a vector in Cylindrical coordinates
\return Components of the vector in Cartesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cylindrical>::CompToCart(const GeoVector& pos, const GeoVector& comp)
{
   double sinphi = sin(pos[1]);
   double cosphi = cos(pos[1]);
   return GeoVector(comp[0] * cosphi - comp[1] * sinphi, comp[0] * sinphi + comp[1] * cosphi, comp.z);
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Spherical coordinates
\param[in] comp Components of a vector in Spherical coordinates
\return Components of the vector in Cartesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Spherical>::CompToCart(const GeoVector& pos, const GeoVector& comp)
{
   double sintheta = sin(pos[1]);
   double costheta = cos(pos[1]);
   double sinphi = sin(pos[2]);
   double cosphi = cos(pos[2]);
   double vrt = comp[0] * sintheta + comp[1] * costheta;
   return GeoVector(vrt * cosphi - comp[2] * sinphi, vrt * sinphi + comp[2] * cosphi, comp[0] * costheta - comp[1] * sintheta);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Cylindrical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cylindrical>::ScaleFactors(const GeoVector& pos)
{
   return GeoVector(1.0, pos[0], 1.0);
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Spherical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Spherical>::ScaleFactors(const GeoVector& pos)
{
   return GeoVector(1.0, pos[0], pos[0] * sin(pos[1]));
};

/*!
\author Vladimir Florinski
\date 10/20/2025
\param[in] pos Position vector in Cylindrical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline double Metric<CoordinateSystem::Cylindrical>::Jacobian(const GeoVector& pos)
{
   return pos[0];
};

/*!
\author Vladimir Florinski
\date 10/20/2025
\param[in] pos Position vector in v coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline double Metric<CoordinateSystem::Spherical>::Jacobian(const GeoVector& pos)
{
   return Sqr(pos[0]) * sin(pos[1]);
};

};

#endif
