/*!
\file boundary_space.cc
\brief Implements several classes representing spatial boundaries
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "boundary_space.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlane methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
template <typename HConfig>
BoundaryPlane<HConfig>::BoundaryPlane(void)
             : BoundaryBase(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Vladimir Florinski
\date 01/25/2021
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryPlane<HConfig>::BoundaryPlane(const std::string& name_in, uint16_t status_in)
             : BoundaryBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryPlane<HConfig>::BoundaryPlane(const BoundaryPlane& other)
             : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryPlane<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin);
   container.Read(norm);
   norm.Normalize();
};

/*!
\author Vladimir Florinski
\date 01/27/2021
*/
template <typename HConfig>
void BoundaryPlane<HConfig>::EvaluateBoundary(void)
{
//         |
// delta<0 | delta>0  
//         |
//         |  norm
//      r0 *---->
//         |
//         |
//         |
//         |

   _delta = (_coords.Pos() - origin) * norm;
   _normal = norm;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlaneAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
template <typename HConfig>
BoundaryPlaneAbsorb<HConfig>::BoundaryPlaneAbsorb(void)
                   : BoundaryPlane(bdy_name, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryPlaneAbsorb<HConfig>::BoundaryPlaneAbsorb(const BoundaryPlaneAbsorb& other)
                   : BoundaryPlane(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryPlaneAbsorb<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryPlane::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlaneReflect methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/25/2021
*/
template <typename HConfig>
BoundaryPlaneReflect<HConfig>::BoundaryPlaneReflect(void)
                    : BoundaryPlane(bdy_name, BOUNDARY_SPACE | BOUNDARY_REFLECT)
{
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryPlaneReflect<HConfig>::BoundaryPlaneReflect(const BoundaryPlaneReflect& other)
                    : BoundaryPlane(other)
{
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryPlaneReflect<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryPlane::SetupBoundary(false);
   if (max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlanePass methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/25/2021
*/
template <typename HConfig>
BoundaryPlanePass<HConfig>::BoundaryPlanePass(void)
                 : BoundaryPlane(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryPlanePass<HConfig>::BoundaryPlanePass(const BoundaryPlanePass& other)
                 : BoundaryPlane(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryPlanePass<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryPlane::SetupBoundary(false);
   if (max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBox methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
*/
template <typename HConfig>
BoundaryBox<HConfig>::BoundaryBox(void)
           : BoundaryBase(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryBox<HConfig>::BoundaryBox(const std::string& name_in, uint16_t status_in)
           : BoundaryBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryBox<HConfig>::BoundaryBox(const BoundaryBox& other)
           : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryBox<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(corners[0]);
   container.Read(normals[0]);
   container.Read(normals[1]);
   container.Read(normals[2]);

// The vectors in the array "normals" are initially given by user as any three edges defining the box, which means their magnitudes are relevant for computing the corner opposite to the one given as an input.
   corners[1] = corners[0] + normals[0] + normals[1] + normals[2];

// After computing the corner opposite to the one given as an input, the edge vectors are normalized to serve as the normal vectors of the sides to which they are each perpendicular.
   normals[0].Normalize();
   normals[1].Normalize();
   normals[2].Normalize();
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
*/
template <typename HConfig>
void BoundaryBox<HConfig>::EvaluateBoundary(void)
{
//                          delta>0
//   +---------------------+
//   |             delta<0 | 
//   ^                     |
//   |                     |
// sides[1]                |
//   |                     |
//   *-- sides[0] -----> --+
// corner[0]

   double delta_tmp;
   _delta = -sp_large;

   for (auto s = 0; s < 3; s++) {
      delta_tmp = -(_coords.Pos() - corners[0]) * normals[s];
      if (delta_tmp > _delta) {
         _delta = delta_tmp;
         _normal = -normals[s];
      };
      delta_tmp = (_coords.Pos() - corners[1]) * normals[s];
      if (delta_tmp > _delta) {
         _delta = delta_tmp;
         _normal = normals[s];
      };
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBoxReflect methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
*/
template <typename HConfig>
BoundaryBoxReflect<HConfig>::BoundaryBoxReflect(void)
                  : BoundaryBox(bdy_name, BOUNDARY_SPACE | BOUNDARY_REFLECT)
{
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryBoxReflect<HConfig>::BoundaryBoxReflect(const BoundaryBoxReflect& other)
                  : BoundaryBox(other)
{
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryBoxReflect<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBox::SetupBoundary(false);
   if (max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphere methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/27/2021
*/
template <typename HConfig>
BoundarySphere<HConfig>::BoundarySphere(void)
              : BoundaryBase(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Vladimir Florinski
\date 01/27/2021
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundarySphere<HConfig>::BoundarySphere(const std::string& name_in, uint16_t status_in)
              : BoundaryBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundarySphere<HConfig>::BoundarySphere(const BoundarySphere& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundarySphere<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin);
   container.Read(radius);
};

/*!
\author Vladimir Florinski
\date 01/25/2022
*/
template <typename HConfig>
void BoundarySphere<HConfig>::EvaluateBoundary(void)
{
//        . -- ~~~ -- .
//    .-~               ~-.
//   /                     \
//  /    delta<0            \
// |                         |  norm
// |            * r0         |---->
// |            |            |
//  \           |           /
//   \        R |          /   delta>0
//    `-.       |       .-'
//        ~- . _V_ . -~

   _delta = (_coords.Pos() - origin).Norm() - radius;
   _normal = UnitVec(_coords.Pos() - origin);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphereAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
template <typename HConfig>
BoundarySphereAbsorb<HConfig>::BoundarySphereAbsorb(void)
                    : BoundarySphere(bdy_name, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundarySphereAbsorb<HConfig>::BoundarySphereAbsorb(const BoundarySphereAbsorb& other)
                    : BoundarySphere(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundarySphereAbsorb<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundarySphere::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphereReflect methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/25/2021
*/
template <typename HConfig>
BoundarySphereReflect<HConfig>::BoundarySphereReflect(void)
                     : BoundarySphere(bdy_name, BOUNDARY_SPACE | BOUNDARY_REFLECT)
{
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundarySphereReflect<HConfig>::BoundarySphereReflect(const BoundarySphereReflect& other)
                     : BoundarySphere(other)
{
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundarySphereReflect<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundarySphere::SetupBoundary(false);
   if (max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRankine methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
template <typename HConfig>
BoundaryRankine<HConfig>::BoundaryRankine(void)
               : BoundaryBase(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryRankine<HConfig>::BoundaryRankine(const std::string& name_in, uint16_t status_in)
               : BoundaryBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/25/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryRankine<HConfig>::BoundaryRankine(const BoundaryRankine& other)
               : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryRankine<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin);
   container.Read(axis);
   container.Read(z_nose);
};

/*!
\author Juan G Alonso Guzman
\date 01/25/2022
*/
template <typename HConfig>
void BoundaryRankine<HConfig>::EvaluateBoundary(void)
{
   GeoVector pos_rel = _coords.Pos() - origin;
   double z = pos_rel * axis;
   double r = pos_rel.Norm();

// TODO Develop a better method for boundary proximity test
   _delta = r - z_nose * sqrt(2.0 / (1.0 + z / r));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRankineAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
template <typename HConfig>
BoundaryRankineAbsorb<HConfig>::BoundaryRankineAbsorb(void)
                     : BoundaryRankine(bdy_name, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryRankineAbsorb<HConfig>::BoundaryRankineAbsorb(const BoundaryRankineAbsorb& other)
                     : BoundaryRankine(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryRankineAbsorb<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryRankine::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryCylinder methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
template <typename HConfig>
BoundaryCylinder<HConfig>::BoundaryCylinder(void)
                : BoundaryBase(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryCylinder<HConfig>::BoundaryCylinder(const std::string& name_in, uint16_t status_in)
                : BoundaryBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryCylinder<HConfig>::BoundaryCylinder(const BoundaryCylinder& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryCylinder<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin);
   container.Read(fa_basis[2]);
   container.Read(radius);

// Find cylinder reference frame where z || symmetry axis
   fa_basis[2].Normalize();
   fa_basis[0] = GetSecondUnitVec(fa_basis[2]);
   fa_basis[1] = fa_basis[2] ^ fa_basis[0];
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
template <typename HConfig>
void BoundaryCylinder<HConfig>::EvaluateBoundary(void)
{
//        . -- ~~~ -- .
//    .-~               ~-.
//   /                     \
//  /    delta<0            \
// |                         |  norm
// |            * r0         |---->
// |            |            |
//  \           |           /
//   \        R |          /   delta>0
//    `-.       |       .-'
//        ~- . _V_ . -~

   GeoVector pos_rel = _coords.Pos() - origin;
   pos_rel.ChangeToBasis(fa_basis);
   _delta = sqrt(Sqr(pos_rel[0]) + Sqr(pos_rel[1])) - radius;
   _normal = pos_rel;
   _normal[2] = 0.0;
   _normal = UnitVec(_normal);
   _normal.ChangeFromBasis(fa_basis);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryCylinderAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
template <typename HConfig>
BoundaryCylinderAbsorb<HConfig>::BoundaryCylinderAbsorb(void)
                      : BoundaryCylinder(bdy_name, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryCylinderAbsorb<HConfig>::BoundaryCylinderAbsorb(const BoundaryCylinderAbsorb& other)
                      : BoundaryCylinder(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryCylinderAbsorb<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryCylinder::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRegion methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
*/
template <typename HConfig>
BoundaryRegion<HConfig>::BoundaryRegion(void)
              : BoundaryBase(bdy_name, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryRegion<HConfig>::BoundaryRegion(const std::string& name_in, uint16_t status_in)
              : BoundaryBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryRegion<HConfig>::BoundaryRegion(const BoundaryRegion& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryRegion<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   // TODO: allow the indicator variable type to vary
//   container.Read(region_ind);
   container.Read(region_val);
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
*/
template <typename HConfig>
void BoundaryRegion<HConfig>::EvaluateBoundary(void)
{
//        
//    _ .. -- ~~~~ -- . _
//   |                    \
//    \                    \
//     \                    \
//      )        delta<0     ) delta>0
//     /                    /
//    /                    /
//   |                    / 
//    ` ~ - . ____ . - ~ *
//

   // TODO: allow the indicator variable type to vary
   _delta = _fields.Iv0() - region_val;
//   _delta = region[region_ind] - region_val;

//"_normal" is impossible to predict for an arbitrary geometry. This means that reflections across this type of boundary are unreliable.
   _normal = gv_zeros;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRegionAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
*/
template <typename HConfig>
BoundaryRegionAbsorb<HConfig>::BoundaryRegionAbsorb(void)
                    : BoundaryRegion(bdy_name, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryRegionAbsorb<HConfig>::BoundaryRegionAbsorb(const BoundaryRegionAbsorb& other)
                    : BoundaryRegion(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryRegionAbsorb<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryRegion::SetupBoundary(false);
   max_crossings = 1;
};

};
