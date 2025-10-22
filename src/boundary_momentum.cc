/*!
\file boundary_momentum.cc
\brief Implements several classes representing momentum boundaries
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "boundary_momentum.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentum methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
template <typename HConfig>
BoundaryMomentum<HConfig>::BoundaryMomentum(void)
                : BoundaryBase(bdy_name, BOUNDARY_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryMomentum<HConfig>::BoundaryMomentum(const std::string& name_in, uint16_t status_in)
                : BoundaryBase(name_in,  status_in)
{
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryMomentum<HConfig>::BoundaryMomentum(const BoundaryMomentum& other)
                : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryMomentum<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   container.Read(momentum);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 02/23/2024
*/
template <typename HConfig>
void BoundaryMomentum<HConfig>::EvaluateBoundary(void)
{
   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline){
// TODO
   }
   else {
      _delta = _coords.AbsMom() - momentum;
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInject methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
template <typename HConfig>
BoundaryMomentumInject<HConfig>::BoundaryMomentumInject(void)
                      : BoundaryMomentum(bdy_name, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 06/26/2025
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BoundaryMomentumInject<HConfig>::BoundaryMomentumInject(const std::string& name_in, uint16_t status_in)
                      : BoundaryMomentum(name_in, status_in)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryMomentumInject<HConfig>::BoundaryMomentumInject(const BoundaryMomentumInject& other)
                      : BoundaryMomentum(other)
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
void BoundaryMomentumInject<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryMomentum::SetupBoundary(false);
   max_crossings = 1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumPass methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/25/2025
*/
template <typename HConfig>
BoundaryMomentumPass<HConfig>::BoundaryMomentumPass(void)
                    : BoundaryMomentum(bdy_name, BOUNDARY_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/25/2025
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryMomentumPass<HConfig>::BoundaryMomentumPass(const BoundaryMomentumPass& other)
                    : BoundaryMomentum(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Juan G Alonso Guzman
\date 06/25/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryMomentumPass<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryMomentum::SetupBoundary(false);
   if (max_crossings == 1) RAISE_BITS(_status, BOUNDARY_TERMINAL);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInjectRestrictSlab methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/26/2025
*/
template <typename HConfig>
BoundaryMomentumInjectRestrictSlab<HConfig>::BoundaryMomentumInjectRestrictSlab(void)
                                  : BoundaryMomentumInject(bdy_name, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
};

/*!
\author Vladimir Florinski
\date 06/26/2025
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryMomentumInjectRestrictSlab<HConfig>::BoundaryMomentumInjectRestrictSlab(const BoundaryMomentumInjectRestrictSlab& other)
                                  : BoundaryMomentumInject(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 06/26/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryMomentumInjectRestrictSlab<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryMomentumInject::SetupBoundary(false);
   container.Read(r0);
   container.Read(r1);
   container.Read(normal);
   normal.Normalize();
};

/*!
\author Vladimir Florinski
\date 06/26/2025
*/
template <typename HConfig>
void BoundaryMomentumInjectRestrictSlab<HConfig>::EvaluateBoundary(void)
{
   BoundaryMomentumInject::EvaluateBoundary();

// If the momentum boundary crossing occurred outside the slab, "_delta_old" is reset to have the same sign as "_delta" to avoid triggering the crossing event.
   if (_delta * _delta_old < 0.0) {
      if (((_coords.Pos() - r0) * normal) * ((_coords.Pos() - r1) * normal) > 0.0) _delta_old = _delta;
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInjectRestrictShell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/27/2025
*/
template <typename HConfig>
BoundaryMomentumInjectRestrictShell<HConfig>::BoundaryMomentumInjectRestrictShell(void)
                                   : BoundaryMomentumInject(bdy_name, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
};

/*!
\author Vladimir Florinski
\date 06/27/2025
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryMomentumInjectRestrictShell<HConfig>::BoundaryMomentumInjectRestrictShell(const BoundaryMomentumInjectRestrictShell& other)
                                   : BoundaryMomentumInject(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
};

/*!
\author Vladimir Florinski
\date 06/27/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryMomentumInjectRestrictShell<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryMomentumInject::SetupBoundary(false);
   container.Read(r0);
   container.Read(r1);
   container.Read(r2);
};

/*!
\author Vladimir Florinski
\date 06/27/2025
*/
template <typename HConfig>
void BoundaryMomentumInjectRestrictShell<HConfig>::EvaluateBoundary(void)
{
   BoundaryMomentumInject::EvaluateBoundary();

// If the momentum boundary crossing occurred outside the slab, "_delta_old" is reset to have the same sign as "_delta" to avoid triggering the crossing event.
   if (_delta * _delta_old < 0.0) {
      auto rdist = (_coords.Pos() - r0).Norm();
      if ((rdist < r1) || (rdist > r2)) _delta_old = _delta;
   };
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMirror methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
template <typename HConfig>
BoundaryMirror<HConfig>::BoundaryMirror(void)
              : BoundaryBase(bdy_name, 0, BOUNDARY_MOMENTUM | BOUNDARY_REFLECT)
{
   max_crossings = -1;
};

/*!
\author Vladimir Florinski
\date 03/07/2022
\param[in] other Object to initialize from
*/
template <typename HConfig>
BoundaryMirror<HConfig>::BoundaryMirror(const BoundaryMirror& other)
              : BoundaryBase(other)
{
   RAISE_BITS(_status, BOUNDARY_MOMENTUM);
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
   max_crossings = -1;
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BoundaryMirror<HConfig>::SetupBoundary(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BoundaryBase::SetupBoundary(false);
   max_crossings = -1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/07/2022
*/
template <typename HConfig>
void BoundaryMirror<HConfig>::EvaluateBoundary(void)
{
// Delta is the parallel momentum component
   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding) {
      _delta = _coords.MomPara();
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Focused) {
      _delta = _coords.AbsMom() * _coords.MomMu();
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz) {
      _delta = _coords.Mom() * _fields.AbsMag();
   }
};


};
