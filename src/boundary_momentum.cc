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
BoundaryMomentum::BoundaryMomentum(void)
                : BoundaryBase("", 0, BOUNDARY_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryMomentum::BoundaryMomentum(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param[in] other Object to initialize from
*/
BoundaryMomentum::BoundaryMomentum(const BoundaryMomentum& other)
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
void BoundaryMomentum::SetupBoundary(bool construct)
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
void BoundaryMomentum::EvaluateBoundary(void)
{
#if (TRAJ_TYPE == TRAJ_LORENTZ) || (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   _delta = _mom.Norm() - momentum;
#elif (TRAJ_TYPE == TRAJ_FOCUSED) || (TRAJ_TYPE == TRAJ_PARKER)
   _delta = _mom[0] - momentum;
#endif
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInject methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
BoundaryMomentumInject::BoundaryMomentumInject(void)
                      : BoundaryMomentum(bnd_name_momentum_inject, 0, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 06/26/2025
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
BoundaryMomentumInject::BoundaryMomentumInject(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                      : BoundaryMomentum(name_in, specie_in, status_in)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] other Object to initialize from
*/
BoundaryMomentumInject::BoundaryMomentumInject(const BoundaryMomentumInject& other)
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
void BoundaryMomentumInject::SetupBoundary(bool construct)
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
BoundaryMomentumPass::BoundaryMomentumPass(void)
                    : BoundaryMomentum(bnd_name_momentum_pass, 0, BOUNDARY_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/25/2025
\param[in] other Object to initialize from
*/
BoundaryMomentumPass::BoundaryMomentumPass(const BoundaryMomentumPass& other)
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
void BoundaryMomentumPass::SetupBoundary(bool construct)
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
BoundaryMomentumInjectRestrictSlab::BoundaryMomentumInjectRestrictSlab(void)
                                  : BoundaryMomentumInject(bnd_name_momentum_inject_restrict_slab, 0, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
};

/*!
\author Vladimir Florinski
\date 06/26/2025
\param[in] other Object to initialize from
*/
BoundaryMomentumInjectRestrictSlab::BoundaryMomentumInjectRestrictSlab(const BoundaryMomentumInjectRestrictSlab& other)
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
void BoundaryMomentumInjectRestrictSlab::SetupBoundary(bool construct)
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
void BoundaryMomentumInjectRestrictSlab::EvaluateBoundary(void)
{
   BoundaryMomentumInject::EvaluateBoundary();

// If the momentum boundary crossing occurred outside the slab, "_delta_old" is reset to have the same sign as "_delta" to avoid triggering the crossing event.
   if (_delta * _delta_old < 0.0) {
      if (((_pos - r0) * normal) * ((_pos - r1) * normal) > 0.0) _delta_old = _delta;
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInjectRestrictShell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/27/2025
*/
BoundaryMomentumInjectRestrictShell::BoundaryMomentumInjectRestrictShell(void)
                                   : BoundaryMomentumInject(bnd_name_momentum_inject_restrict_shell, 0, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
};

/*!
\author Vladimir Florinski
\date 06/27/2025
\param[in] other Object to initialize from
*/
BoundaryMomentumInjectRestrictShell::BoundaryMomentumInjectRestrictShell(const BoundaryMomentumInjectRestrictShell& other)
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
void BoundaryMomentumInjectRestrictShell::SetupBoundary(bool construct)
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
void BoundaryMomentumInjectRestrictShell::EvaluateBoundary(void)
{
   BoundaryMomentumInject::EvaluateBoundary();

// If the momentum boundary crossing occurred outside the slab, "_delta_old" is reset to have the same sign as "_delta" to avoid triggering the crossing event.
   if (_delta * _delta_old < 0.0) {
      auto rdist = (_pos - r0).Norm(); 
      if ((rdist < r1) || (rdist > r2)) _delta_old = _delta;
   };
};

#if (TRAJ_TYPE != TRAJ_PARKER) && (TRAJ_TYPE != TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMirror methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
BoundaryMirror::BoundaryMirror(void)
              : BoundaryBase(bnd_name_mirror, 0, BOUNDARY_MOMENTUM | BOUNDARY_REFLECT)
{
   max_crossings = -1;
};

/*!
\author Vladimir Florinski
\date 03/07/2022
\param[in] other Object to initialize from
*/
BoundaryMirror::BoundaryMirror(const BoundaryMirror& other)
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
void BoundaryMirror::SetupBoundary(bool construct)
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
void BoundaryMirror::EvaluateBoundary(void)
{
// Delta is the parallel momentum component
#if (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   _delta = _mom[2];
#elif TRAJ_TYPE == TRAJ_FOCUSED
   _delta = _mom[0] * _mom[1];
#elif TRAJ_TYPE == TRAJ_LORENTZ
   _delta = _mom * bhat;
#endif
};

#endif

};
