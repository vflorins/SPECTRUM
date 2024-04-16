/*!
\file boundary_space.cc
\brief Implements several classes representing spatial boundaries
\author Vladimir Florinski
\author Juan G Alonso Guzman

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
BoundaryPlane::BoundaryPlane(void)
             : BoundaryBase("", 0, BOUNDARY_SPACE)
{
};

/*!
\author Vladimir Florinski
\date 01/25/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryPlane::BoundaryPlane(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
             : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] other Object to initialize from
*/
BoundaryPlane::BoundaryPlane(const BoundaryPlane& other)
             : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryPlane::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin.Data());
   container.Read(norm.Data());
   norm.Normalize();
};

/*!
\author Vladimir Florinski
\date 01/27/2021
*/
void BoundaryPlane::EvaluateBoundary(void)
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

   _delta = (_pos - origin) * norm;
   _normal = norm;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlaneAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
BoundaryPlaneAbsorb::BoundaryPlaneAbsorb(void)
                   : BoundaryPlane(bnd_name_plane_absorb, 0, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from
*/
BoundaryPlaneAbsorb::BoundaryPlaneAbsorb(const BoundaryPlaneAbsorb& other)
                   : BoundaryPlane(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryPlaneAbsorb::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryPlane::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlaneReflect methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/25/2021
*/
BoundaryPlaneReflect::BoundaryPlaneReflect(void)
                    : BoundaryPlane(bnd_name_plane_reflect, 0, BOUNDARY_SPACE | BOUNDARY_REFLECT)
{
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from
*/
BoundaryPlaneReflect::BoundaryPlaneReflect(const BoundaryPlaneReflect& other)
                    : BoundaryPlane(other)
{
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryPlaneReflect::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryPlane::SetupBoundary(false);
   if(max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlanePass methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/25/2021
*/
BoundaryPlanePass::BoundaryPlanePass(void)
                 : BoundaryPlane(bnd_name_plane_pass, 0, BOUNDARY_SPACE)
{
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param[in] other Object to initialize from
*/
BoundaryPlanePass::BoundaryPlanePass(const BoundaryPlanePass& other)
                 : BoundaryPlane(other)
{
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/21/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryPlanePass::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryPlane::SetupBoundary(false);
   if(max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBox methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
*/
BoundaryBox::BoundaryBox(void)
           : BoundaryBase("", 0, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryBox::BoundaryBox(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
           : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param[in] other Object to initialize from
*/
BoundaryBox::BoundaryBox(const BoundaryBox& other)
           : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryBox::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(corners[0].Data());
   container.Read(normals[0].Data());
   container.Read(normals[1].Data());
   container.Read(normals[2].Data());
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
void BoundaryBox::EvaluateBoundary(void)
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
   _delta = -large;
   for(int s = 0; s < 3; s++) {
      delta_tmp = -(_pos - corners[0]) * normals[s];
      if(delta_tmp > _delta) {
         _delta = delta_tmp;
         _normal = -normals[s];
      };
      delta_tmp = (_pos - corners[1]) * normals[s];
      if(delta_tmp > _delta) {
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
BoundaryBoxReflect::BoundaryBoxReflect(void)
                  : BoundaryBox(bnd_name_box_reflect, 0, BOUNDARY_SPACE | BOUNDARY_REFLECT)
{
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param[in] other Object to initialize from
*/
BoundaryBoxReflect::BoundaryBoxReflect(const BoundaryBoxReflect& other)
                  : BoundaryBox(other)
{
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 11/02/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryBoxReflect::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBox::SetupBoundary(false);
   if(max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphere methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/27/2021
*/
BoundarySphere::BoundarySphere(void)
              : BoundaryBase("", 0, BOUNDARY_SPACE)
{
};

/*!
\author Vladimir Florinski
\date 01/27/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundarySphere::BoundarySphere(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
              : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] other Object to initialize from
*/
BoundarySphere::BoundarySphere(const BoundarySphere& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundarySphere::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin.Data());
   container.Read(&radius);
};

/*!
\author Vladimir Florinski
\date 01/25/2022
*/
void BoundarySphere::EvaluateBoundary(void)
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

   _delta = (_pos - origin).Norm() - radius;
   _normal = UnitVec(_pos - origin);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphereAbsorb methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/17/2020
*/
BoundarySphereAbsorb::BoundarySphereAbsorb(void)
                    : BoundarySphere(bnd_name_sphere_absorb, 0, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param[in] other Object to initialize from
*/
BoundarySphereAbsorb::BoundarySphereAbsorb(const BoundarySphereAbsorb& other)
                    : BoundarySphere(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundarySphereAbsorb::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundarySphere::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphereReflect methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/25/2021
*/
BoundarySphereReflect::BoundarySphereReflect(void)
                     : BoundarySphere(bnd_name_sphere_reflect, 0, BOUNDARY_SPACE | BOUNDARY_REFLECT)
{
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param[in] other Object to initialize from
*/
BoundarySphereReflect::BoundarySphereReflect(const BoundarySphereReflect& other)
                     : BoundarySphere(other)
{
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundarySphereReflect::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundarySphere::SetupBoundary(false);
   if(max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRankine methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
BoundaryRankine::BoundaryRankine(void)
               : BoundaryBase("", 0, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryRankine::BoundaryRankine(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
               : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/25/2022
\param[in] other Object to initialize from
*/
BoundaryRankine::BoundaryRankine(const BoundaryRankine& other)
               : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 01/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryRankine::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin.Data());
   container.Read(axis.Data());
   container.Read(&z_nose);
};

/*!
\author Juan G Alonso Guzman
\date 01/25/2022
*/
void BoundaryRankine::EvaluateBoundary(void)
{
   GeoVector pos_rel = _pos - origin;
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
BoundaryRankineAbsorb::BoundaryRankineAbsorb(void)
                     : BoundaryRankine(bnd_name_rankine_absorb, 0, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] other Object to initialize from
*/
BoundaryRankineAbsorb::BoundaryRankineAbsorb(const BoundaryRankineAbsorb& other)
                     : BoundaryRankine(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryRankineAbsorb::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryRankine::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryCylinder methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
BoundaryCylinder::BoundaryCylinder(void)
                : BoundaryBase("", 0, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryCylinder::BoundaryCylinder(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] other Object to initialize from
*/
BoundaryCylinder::BoundaryCylinder(const BoundaryCylinder& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryCylinder::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(origin.Data());
   container.Read(fa_basis[2].Data());
   container.Read(&radius);

// Find cylinder reference frame where z || symmetry axis
   fa_basis[2].Normalize();
   fa_basis[0] = GetSecondUnitVec(fa_basis[2]);
   fa_basis[1] = fa_basis[2] ^ fa_basis[0];
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
void BoundaryCylinder::EvaluateBoundary(void)
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

   GeoVector pos_rel = _pos - origin;
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
BoundaryCylinderAbsorb::BoundaryCylinderAbsorb(void)
                      : BoundaryCylinder(bnd_name_cylinder_absorb, 0, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] other Object to initialize from
*/
BoundaryCylinderAbsorb::BoundaryCylinderAbsorb(const BoundaryCylinderAbsorb& other)
                      : BoundaryCylinder(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryCylinderAbsorb::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryCylinder::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRegion methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
*/
BoundaryRegion::BoundaryRegion(void)
              : BoundaryBase("", 0, BOUNDARY_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryRegion::BoundaryRegion(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
              : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param[in] other Object to initialize from
*/
BoundaryRegion::BoundaryRegion(const BoundaryRegion& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryRegion::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(&region_ind);
   container.Read(&region_val);
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
*/
void BoundaryRegion::EvaluateBoundary(void)
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

   _delta = region[region_ind] - region_val;
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
BoundaryRegionAbsorb::BoundaryRegionAbsorb(void)
                    : BoundaryRegion(bnd_name_region_absorb, 0, BOUNDARY_SPACE | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param[in] other Object to initialize from
*/
BoundaryRegionAbsorb::BoundaryRegionAbsorb(const BoundaryRegionAbsorb& other)
                    : BoundaryRegion(other)
{
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 12/07/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BoundaryRegionAbsorb::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BoundaryRegion::SetupBoundary(false);
   max_crossings = 1;
};

};
