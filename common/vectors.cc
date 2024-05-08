/*!
\file vectors.cc
\brief Defines some non-trivial functions operating on vectors
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Methods of GeoVector
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/01/2018
\return Polar angle
*/
SPECTRUM_DEVICE_FUNC double GeoVector::Theta(void)
{
   double theta = atan2(sqrt(Sqr(data[0]) + Sqr(data[1])), data[2]);
   return (theta >= 0.0 ? theta : theta + twopi);
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\return Azimuthal angle
*/
SPECTRUM_DEVICE_FUNC double GeoVector::Phi(void)
{
   double phi = atan2(data[1], data[0]);
   return (phi >= 0.0 ? phi : phi + twopi);
};

/*!
\author Vladimir Florinski
\date 07/23/2019
*/
SPECTRUM_DEVICE_FUNC void GeoVector::RTP_XYZ(void)
{
   double xy = data[0] * sin(data[1]);
   double z =  data[0] * cos(data[1]);
   data[0] = xy * cos(data[2]);
   data[1] = xy * sin(data[2]);
   data[2] = z;
};

/*!
\author Vladimir Florinski
\date 07/23/2019
*/
SPECTRUM_DEVICE_FUNC void GeoVector::XYZ_RTP(void)
{
   double r = Norm();
   double theta = acos(data[2] / r);
   double phi = atan2(data[1], data[0]);
   if(phi < 0.0) phi += twopi;
   data[0] = r;
   data[1] = theta;
   data[2] = phi;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] sintheta \f$\sin\theta\f$
\param[in] costheta \f$\cos\theta\f$
\param[in] sinphi   \f$\sin\phi\f$
\param[in] cosphi   \f$\cos\phi\f$
*/
SPECTRUM_DEVICE_FUNC void GeoVector::ToSpherical(double sintheta, double costheta, double sinphi, double cosphi)
{
   double uxy =  data[0] * cosphi + data[1] * sinphi;
   double up  = -data[0] * sinphi + data[1] * cosphi;
   data[0] = uxy * sintheta + data[2] * costheta;
   data[1] = uxy * costheta - data[2] * sintheta;
   data[2] = up;
};

/*!
\author Vladimir Florinski
\date 05/01/2018
\param[in] sintheta \f$\sin\theta\f$
\param[in] costheta \f$\cos\theta\f$
\param[in] sinphi   \f$\sin\phi\f$
\param[in] cosphi   \f$\cos\phi\f$
*/
SPECTRUM_DEVICE_FUNC void GeoVector::ToCartesian(double sintheta, double costheta, double sinphi, double cosphi)
{
   double urt = data[0] * sintheta + data[1] * costheta;
   double uz  = data[0] * costheta - data[1] * sintheta;
   data[0] = urt * cosphi - data[2] * sinphi;
   data[1] = urt * sinphi + data[2] * cosphi;
   data[2] = uz;
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] n    Axis of rotation \f$\mathbf{n}\f$
\param[in] sina Sine of the rotation angle \f$\sin\alpha\f$
\param[in] cosa Cosine of the rotation angle \f$\cos\alpha\f$

The formula is \f$\mathbf{v}'=\cos\alpha\mathbf{v}+(1-\cos\alpha)(\mathbf{v}\cdot\mathbf{n})\mathbf{n}+\sin\alpha(\mathbf{n}\times\mathbf{v})\f$
*/
SPECTRUM_DEVICE_FUNC void GeoVector::Rotate(const GeoVector& n, double sina, double cosa)
{
   *this = *this * cosa + n * (*this * n) * (1.0 - cosa) + (n ^ *this) * sina;
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] axis  Axis of rotation (may not be a unit vector) \f$\sim\mathbf{n}\f$
\param[in] angle Angle of rotation \f$\alpha\f$
*/
SPECTRUM_DEVICE_FUNC void GeoVector::Rotate(const GeoVector& axis, double angle)
{
   GeoVector n = UnitVec(axis);
   Rotate(n, sin(angle), cos(angle));
};

/*!
\author Vladimir Florinski
\date 02/13/2024
\param[in] axis  Unit vector whose normal space we project to \f$\sim\mathbf{n}\f$
*/
SPECTRUM_DEVICE_FUNC void GeoVector::SubtractParallel(const GeoVector& axis)
{
   *this -= (*this * axis) * axis;
};

/*!
\author Vladimir Florinski
\date 09/12/2019
\param[in] basis Three new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoVector::ChangeToBasis(const GeoVector* basis)
{
   GeoVector vect_dif;
   vect_dif[0] = data[0] * basis[0].data[0] + data[1] * basis[0].data[1] + data[2] * basis[0].data[2];
   vect_dif[1] = data[0] * basis[1].data[0] + data[1] * basis[1].data[1] + data[2] * basis[1].data[2];
   vect_dif[2] = data[0] * basis[2].data[0] + data[1] * basis[2].data[1] + data[2] * basis[2].data[2];
   std::memcpy(data, vect_dif.data, 3 * SZDBL);
};

/*!
\author Vladimir Florinski
\date 09/09/2022
\param[in] basis Three new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoVector::ChangeToBasis2(const GeoVector* basis)
{
   GeoVector vect_dif;
   vect_dif[0] = data[0] * Sqr(basis[0].data[0]) + data[1] * Sqr(basis[0].data[1]) + data[2] * Sqr(basis[0].data[2]);
   vect_dif[1] = data[0] * Sqr(basis[1].data[0]) + data[1] * Sqr(basis[1].data[1]) + data[2] * Sqr(basis[1].data[2]);
   vect_dif[2] = data[0] * Sqr(basis[2].data[0]) + data[1] * Sqr(basis[2].data[1]) + data[2] * Sqr(basis[2].data[2]);
   std::memcpy(data, vect_dif.data, 3 * SZDBL);
};

/*!
\author Vladimir Florinski
\date 09/12/2019
\param[in] basis Three new basis vectors
*/
SPECTRUM_DEVICE_FUNC void GeoVector::ChangeFromBasis(const GeoVector* basis)
{
   GeoVector vect_std;
   vect_std[0] = data[0] * basis[0].data[0] + data[1] * basis[1].data[0] + data[2] * basis[2].data[0];
   vect_std[1] = data[0] * basis[0].data[1] + data[1] * basis[1].data[1] + data[2] * basis[2].data[1];
   vect_std[2] = data[0] * basis[0].data[2] + data[1] * basis[1].data[2] + data[2] * basis[2].data[2];
   std::memcpy(data, vect_std.data, 3 * SZDBL);
};

/*!
\author Vladimir Florinski
\date 03/06/2019
\param[in]  projection Type of projection (gnomonic, orthographic, equidistant)
\param[in]  normal     Normal to the plane
\param[in]  north      Up direction in the projection plane
\param[out] xi         The \f$\xi\f$ local coordinate of the point
\param[out] eta        The \f$\eta\f$ local coordinate of the point

The projections are given by:

gnomonic     \f$\xi=\tan\theta\cos\phi, \eta=\tan\theta\sin\phi\f$
orthographic \f$\xi=\sin\theta\cos\phi, \eta=\sin\theta\sin\phi\f$
equidistant  \f$\xi=\theta\cos\phi,     \eta=\theta\sin\phi\f$
*/
SPECTRUM_DEVICE_FUNC void GeoVector::ProjectToPlane(int projection, const GeoVector& normal, const GeoVector& north, double& xi, double& eta) const
{
   double rho, cos_dist, cp;

// Get the scalar product and check for bounds
   cos_dist = normal * *this;
   cos_dist = fmax(-1.0, fmin(1.0, cos_dist));

// Distance from the center
   if(projection == 1) rho = tan(acos(cos_dist));
   else if(projection == 2) rho = sin(acos(cos_dist));
   else if(projection == 3) rho = acos(cos_dist);
   else rho = 0.0;

// Azimuthal angle cosine is the same for all projections.
   cp = VertexAngle(north, normal, *this);
   xi = rho * cp;
   eta = rho * sqrt(fmax(0.0, 1.0 - Sqr(cp)));

// The vector product of center and North gives the West direction. The triple product is negative if v points East.
   if(*this * (normal ^ north) < 0.0) eta = -eta;
};

/*!
\author Vladimir Florinski
\date 05/01/2020
\param[in] n_verts Number of vertices forming the wedge
\param[in] verts   Vertices, must be CC ordered
\param[in] tol     Tolerance (include points if the scalar product is small)
\return True if the vector is in the interior, including the boundary, of the wedge (solid angle)
*/
SPECTRUM_DEVICE_FUNC bool GeoVector::InsideWedge(int n_verts, const GeoVector* verts, double tol) const
{
   double scp;
   GeoVector edge_norm;

// Check the sign of the projection of "v" onto planes through each vector pair. This relies on CC ordering of the vertices.
   for(auto iv = 0; iv < n_verts; iv++) {
      edge_norm = verts[iv] ^ verts[(iv + 1) % n_verts];
      scp = *this * edge_norm;
      if(tol > 0.0) {
         if((scp < 0.0) && (scp < -tol)) return false;
      }
      else {
         if(scp < 0.0) return false;
         else if(scp < -tol) return false;
      };
   };
   return true;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on vectors
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_r right operand \f$\mathbf{v}_2\f$
\return \f$\sim(\mathbf{v}_1+\mathbf{v}_2)/2\f$
*/
SPECTRUM_DEVICE_FUNC GeoVector Bisect(const GeoVector& vect_l, const GeoVector& vect_r)
{
   GeoVector sum = vect_l + vect_r;
   return sum.Normalize();
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_m middle operand \f$\mathbf{v}_2\f$
\param[in] vect_r right operand \f$\mathbf{v}_3\f$
\return Area of a triangle defined with three vertices
*/
SPECTRUM_DEVICE_FUNC double TriangleArea(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r)
{
   return 0.5 * ((vect_m - vect_l) ^ (vect_r - vect_l)).Norm();
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_m middle operand \f$\mathbf{v}_2\f$
\param[in] vect_r right operand \f$\mathbf{v}_3\f$
\return \f$\sim(\mathbf{v}_2-\mathbf{v}_1)\times(\mathbf{v}_3-\mathbf{v}_1)\f$
*/
SPECTRUM_DEVICE_FUNC GeoVector PlaneNormal(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r)
{
   GeoVector n = (vect_m - vect_l) ^ (vect_r - vect_l);
   return n.Normalize();
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_m middle operand \f$\mathbf{v}_2\f$
\param[in] vect_r right operand \f$\mathbf{v}_3\f$
\return \f$(\mathbf{v}_1+\mathbf{v}_2+\mathbf{v}_3)/3\f$
*/
SPECTRUM_DEVICE_FUNC GeoVector TriangleCenter(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r)
{
   return (vect_l + vect_m + vect_r) / 3.0;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_m middle operand \f$\mathbf{v}_2\f$
\param[in] vect_r right operand \f$\mathbf{v}_3\f$
\return \f$(\mathbf{n}\cdot\mathbf{v}_1)\mathbf{n}\f$
*/
SPECTRUM_DEVICE_FUNC GeoVector CircumCenter(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r)
{
   GeoVector n = (vect_m - vect_l) ^ (vect_r - vect_l);
   n.Normalize();
   return (n * vect_l) * n;
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_r right operand \f$\mathbf{v}_2\f$
\param[in] frac   The fraction of the total angle, measured from the left
\return \f$\sim\alpha\mathbf{v}_1+\beta\mathbf{v}_2\f$

The function returns the unit vector in the v1-v2 plane that makes an angle with v1 that is a fraction "frac" of the total.
*/
SPECTRUM_DEVICE_FUNC GeoVector DivideAngle(const GeoVector& vect_l, const GeoVector& vect_r, double frac)
{
   double theta, costheta;

   costheta = vect_l * vect_r;
   theta = acos(costheta);
   return (sin((1.0 - frac) * theta) * vect_l + sin(frac * theta) * vect_r) / sqrt(1.0 - Sqr(costheta));
};

/*!
\author Vladimir Florinski
\date 05/14/2021
\param[in] vect1 First triangle vertex
\param[in] vect2 Second triangle vertex
\param[in] vect3 Third triangle vertex
\param[in] vect  Vector in the plane of the triangle
\return Relative areas of the three sub-triangles (add up to unity)
*/
SPECTRUM_DEVICE_FUNC GeoVector BarycentricCoords(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3, const GeoVector& vect)
{
   GeoVector areas;

   areas[0] = TriangleArea(vect, vect2, vect3);
   areas[1] = TriangleArea(vect, vect3, vect1);
   areas[2] = TriangleArea(vect, vect1, vect2);

   return areas / areas.Sum();
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in]  v  Three vertices of the trtiangle
\param[out] C0 Zeroth moment
\param[out] C1 First moments
\param[out] C2 Second moments
\param[out] C3 Third moments
\param[out] C4 Fourth moments

Ref: flat_quad_moments_3d.pdf
*/
SPECTRUM_DEVICE_FUNC void TriangleMoments(const GeoVector v[3], double& C0, double C1[3], double C2[3][3], double C3[3][3][3], double C4[3][3][3][3])
{
   int i, j, k, l, p;

   C0 = TriangleArea(v[0], v[1], v[2]);

// First order moments
   for(i = 0; i < 3; i++) {
      C1[i] = 0.0;
      for(p = 0; p < 3; p++) {
         C1[i] += C0 / 3.0 * v[p][i];
      };
   };

// Second order moments
   for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
         C2[i][j] = 9.0 * C1[i] * C1[j] / Sqr(C0);
         for(p = 0; p < 3; p++) {
            C2[i][j] += v[p][i] * v[p][j];
         };
         C2[i][j] *= C0 / 12.0;
      };
   };

// Third order moments
   for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
         for(k = 0; k < 3; k++) {
            C3[i][j][k] = -27.0 * C1[i] * C1[j] * C1[k] / Cube(C0)
               + 18.0 / Sqr(C0) * (C1[i] * C2[j][k] + C1[j] * C2[i][k] + C1[k] * C2[i][j]);
            for(p = 0; p < 3; p++) {
               C3[i][j][k] += v[p][i] * v[p][j] * v[p][k];
            };
            C3[i][j][k] *= C0 / 30.0;
         };
      };
   };

// Fourth order moments
   for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
         for(k = 0; k < 3; k++) {
            for(l = 0; l < 3; l++) {
               C4[i][j][k][l] = -81.0 * C1[i] * C1[j] * C1[k] * C1[l] / Quad(C0)
                  + 18.0 / Cube(C0) * (C1[i] * C1[j] * C2[k][l] + C1[i] * C1[k] * C2[j][l]
                                     + C1[i] * C1[l] * C2[j][k] + C1[j] * C1[k] * C2[i][l]
                                     + C1[j] * C1[l] * C2[i][k] + C1[k] * C1[l] * C2[i][j]);
               for(p = 0; p < 3; p++) {
                  C4[i][j][k][l] += v[p][i] * v[p][j] * v[p][k] * v[p][l];
               };
               C4[i][j][k][l] *= C0 / 30.0;
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] vect_l left operand \f$\mathbf{v}_1\f$
\param[in] vect_m middle operand \f$\mathbf{v}_2\f$
\param[in] vect_r right operand \f$\mathbf{v}_3\f$
\return \f$[(\mathbf{v}_3\cdot\mathbf{v}_1)-(\mathbf{v}_1\cdot\mathbf{v}_2)(\mathbf{v}_2\cdot\mathbf{v}_3)]/[|\mathbf{v}_1\times\mathbf{v}_2||\mathbf{v}_2\times\mathbf{v}_3|]\f$
*/
SPECTRUM_DEVICE_FUNC double VertexAngle(const GeoVector& vect_l, const GeoVector& vect_m, const GeoVector& vect_r)
{
   return ((vect_r * vect_l) - (vect_l * vect_m) * (vect_m * vect_r)) / ((vect_l ^ vect_m).Norm() * (vect_m ^ vect_r).Norm());
};

/*!
\author Vladimir Florinski
\date 01/08/2020
\param[in] vect1 first vertex \f$\mathbf{v}_1\f$
\param[in] vect2 second vertex \f$\mathbf{v}_2\f$
\param[in] vect3 third vertex \f$\mathbf{v}_3\f$
\return Area of the spherical triangle on a unit sphere

Ref: Tuynman, G. M., "Areas of spherical and hyperbolic triangles in terms of their midpoints", https://arxiv.org/pdf/1307.2567.pdf
*/
SPECTRUM_DEVICE_FUNC double SphTriArea(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3)
{
   double v1v2 = vect1 * vect2;
   double v2v3 = vect2 * vect3;
   double v3v1 = vect3 * vect1;
   return 2.0 * acos((1.0 + v1v2 + v2v3 + v3v1) / sqrt(2.0 * (1.0 + v1v2) * (1.0 + v2v3) * (1.0 + v3v1)));
};

/*!
\author Vladimir Florinski
\date 01/08/2020
\param[in] vect1 first vertex \f$\mathbf{v}_1\f$
\param[in] vect2 second vertex \f$\mathbf{v}_2\f$
\param[in] vect3 third vertex \f$\mathbf{v}_3\f$
\param[in] vect4 fourth vertex \f$\mathbf{v}_4\f$
\return Area of the spherical quadrilateral on a unit sphere
*/
SPECTRUM_DEVICE_FUNC double SphQuadArea(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3, const GeoVector& vect4)
{
   return SphTriArea(vect1, vect2, vect3) + SphTriArea(vect3, vect4, vect1);
};

/*!
\author Vladimir Florinski
\date 01/08/2020
\param[in] vect1 first vertex \f$\mathbf{v}_1\f$
\param[in] vect2 second vertex \f$\mathbf{v}_2\f$
\param[in] vect3 third vertex \f$\mathbf{v}_3\f$
\return Center of mass of the spherical triangle

Ref: Brock, J. E., "The centroid and inertia tensor for a spherical triangle"
*/
SPECTRUM_DEVICE_FUNC GeoVector SphTriCenter(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3)
{
// Compute vectors normal to the sides of the triangle on the sphere
   GeoVector v1xv2 = (vect1 ^ vect2).Normalize();
   GeoVector v2xv3 = (vect2 ^ vect3).Normalize();
   GeoVector v3xv1 = (vect3 ^ vect1).Normalize();

   return (acos(vect1 * vect2) * v1xv2 + acos(vect2 * vect3) * v2xv3 + acos(vect3 * vect1) * v3xv1)
          / (2.0 * SphTriArea(vect1, vect2, vect3));
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] vect1 Unit vector to point 1 on the first circle
\param[in] vect2 Unit vector to point 2 on the first circle
\param[in] vect3 Unit vector to point 1 on the second circle
\param[in] vect4 Unit vector to point 2 on the second circle
\return Unit vector to the intersection of circles 1 and 2
*/
SPECTRUM_DEVICE_FUNC GeoVector GreatCircleInt(const GeoVector& vect1, const GeoVector& vect2, const GeoVector& vect3, const GeoVector& vect4)
{
// Vector product of the normals gives the common line
   GeoVector vect5 = (vect1 ^ vect2) ^ (vect3 ^ vect4);

// Choose the direction that is closest to point 1 (angle < 90 deg).
   if(vect1 * vect5 < 0.0) vect5 = -vect5;
   return vect5.Normalize();
};

/*!
\author Vladimir Florinski
\date 09/12/2019
\param[in] first First unit vector (the new z-direction)
\return New unit vector normal to "first"
*/
SPECTRUM_DEVICE_FUNC GeoVector GetSecondUnitVec(const GeoVector& first)
{
   int max_ang;

// Find the unit vector making the largest angle with "first" (smallest component of "first")
   if(fabs(first[0]) < fabs(first[1])) {
      if(fabs(first[0]) < fabs(first[2])) max_ang = 0;
      else max_ang = 2;
   }
   else {
      if(fabs(first[1]) < fabs(first[2])) max_ang = 1;
      else max_ang = 2;
   };

   return (first ^ cart_unit_vec[max_ang]).Normalize();
};

/*!
\author Vladimir Florinski
\date 11/02/2020
\param[in] vect_l First vector
\param[in] vect_r Second vector
\return Unit vector normal to "vect_l" in the plane of "vect_l"-"vect_r"
*/
SPECTRUM_DEVICE_FUNC GeoVector GetSecondUnitVec(const GeoVector& vect_l, const GeoVector& vect_r)
{
   GeoVector normal;
   normal = vect_r - ((vect_l * vect_r) / vect_l.Norm2()) * vect_l;
   return normal.Normalize();
};

};
