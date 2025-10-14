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
protected:

public:

//! Convert a position vector from Cartesian to curvilinear
   SPECTRUM_DEVICE_FUNC static GeoVector PosToCurv(const GeoVector& pos);

//! Convert a position vector from curvilinear to Cartesian 
   SPECTRUM_DEVICE_FUNC static GeoVector PosToCart(const GeoVector& pos);

//! Convert the components of a vector from Cartesian to curvilinear
   SPECTRUM_DEVICE_FUNC static GeoVector CompToCurv(const GeoVector& pos, const GeoVector& comp);

//! Convert the components of a vector from curvilinear to Cartesian 
   SPECTRUM_DEVICE_FUNC static GeoVector CompToCart(const GeoVector& pos, const GeoVector& comp);

//! Scale factors 
   SPECTRUM_DEVICE_FUNC static GeoVector ScaleFactors(const GeoVector& pos);
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Caretesian coordinates
\return Position in Caretesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cartesian>::PosToCurv(const GeoVector& pos)
{
   return pos;
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Caretesian coordinates
\return Position in Cylindrical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cylindrical>::PosToCurv(const GeoVector& pos)
{
   double phi = atan2(pos.y, pos.x);
   if(phi < 0.0) phi += M_2PI;
   return GeoVector(sqrt(Sqr(pos.x) + Sqr(pos.y)), phi, pos.z);
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Caretesian coordinates
\return Position in Spherical coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Spherical>::PosToCurv(const GeoVector& pos)
{
   double r = pos.Norm();
   double theta = acos(pos.z / r);
   double phi = atan2(pos.y, pos.x);
   if(phi < 0.0) phi += M_2PI;
   return GeoVector(r, theta, phi);
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Caretesian coordinates
\return Position in Caretesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cartesian>::PosToCart(const GeoVector& pos)
{
   return pos;
};

/*!
\author Vladimir Florinski
\date 10/13/2025
\param[in] pos Position in Cylindrical coordinates
\return Position in Caretesian coordinates
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
\return Position in Caretesian coordinates
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
\param[in] pos Position vector (unused)
\param[in] comp Components of a vector in Caretesian coordinates
\return Components of the vector in Caretesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cartesian>::CompToCurv(const GeoVector& pos, const GeoVector& comp)
{
   return comp;
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Cylindrical coordinates
\param[in] comp Components of a vector in Caretesian coordinates
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
\param[in] comp Components of a vector in Caretesian coordinates
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
\param[in] pos Position vector (unused)
\param[in] comp Components of a vector in Caretesian coordinates
\return Components of the vector in Caretesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cartesian>::CompToCart(const GeoVector& pos, const GeoVector& comp)
{
   return comp;
};

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Cylindrical coordinates
\param[in] comp Components of a vector in Cylindrical coordinates
\return Components of the vector in Caretesian coordinates
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
\return Components of the vector in Caretesian coordinates
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

/*!
\author Vladimir Florinski
\date 10/14/2025
\param[in] pos Position vector in Caretesian coordinates
*/
template <>
SPECTRUM_DEVICE_FUNC inline GeoVector Metric<CoordinateSystem::Cartesian>::ScaleFactors(const GeoVector& pos)
{
   return gv_ones;
};

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

};

#endif