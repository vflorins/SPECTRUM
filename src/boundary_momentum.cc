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
template <typename Trajectory>
BoundaryMomentum<Trajectory>::BoundaryMomentum(void)
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
template <typename Trajectory>
BoundaryMomentum<Trajectory>::BoundaryMomentum(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                : BoundaryBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 03/25/2022
\param[in] other Object to initialize from
*/
template <typename Trajectory>
BoundaryMomentum<Trajectory>::BoundaryMomentum(const BoundaryMomentum& other)
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
template <typename Trajectory>
void BoundaryMomentum<Trajectory>::SetupBoundary(bool construct)
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
template <typename Trajectory>
void BoundaryMomentum<Trajectory>::EvaluateBoundary(void)
{
   if constexpr (std::same_as<Trajectory, TrajectoryLorentz<Fields>> || std::derived_from<Trajectory, TrajectoryGuidingBase<Trajectory, Fields>>) {
      _delta = _mom.Norm() - momentum;
   }
   else if constexpr (std::same_as<Trajectory, TrajectoryFocused<Fields>> || std::same_as<Trajectory, TrajectoryParker<Fields>>) {
      _delta = _mom[0] - momentum;
   }
   else if constexpr (std::derived_from<Trajectory, TrajectoryFieldlineBase<Trajectory, Fields>>){
// TODO
      ;
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInject methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
*/
template <typename Trajectory>
BoundaryMomentumInject<Trajectory>::BoundaryMomentumInject(void)
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
template <typename Trajectory>
BoundaryMomentumInject<Trajectory>::BoundaryMomentumInject(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                      : BoundaryMomentum(name_in, specie_in, status_in)
{
   max_crossings = 1;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2024
\param[in] other Object to initialize from
*/
template <typename Trajectory>
BoundaryMomentumInject<Trajectory>::BoundaryMomentumInject(const BoundaryMomentumInject& other)
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
template <typename Trajectory>
void BoundaryMomentumInject<Trajectory>::SetupBoundary(bool construct)
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
template <typename Trajectory>
BoundaryMomentumPass<Trajectory>::BoundaryMomentumPass(void)
                    : BoundaryMomentum(bnd_name_momentum_pass, 0, BOUNDARY_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/25/2025
\param[in] other Object to initialize from
*/
template <typename Trajectory>
BoundaryMomentumPass<Trajectory>::BoundaryMomentumPass(const BoundaryMomentumPass& other)
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
template <typename Trajectory>
void BoundaryMomentumPass<Trajectory>::SetupBoundary(bool construct)
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
template <typename Trajectory>
BoundaryMomentumInjectRestrictSlab<Trajectory>::BoundaryMomentumInjectRestrictSlab(void)
                                  : BoundaryMomentumInject(bnd_name_momentum_inject_restrict_slab, 0, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
};

/*!
\author Vladimir Florinski
\date 06/26/2025
\param[in] other Object to initialize from
*/
template <typename Trajectory>
BoundaryMomentumInjectRestrictSlab<Trajectory>::BoundaryMomentumInjectRestrictSlab(const BoundaryMomentumInjectRestrictSlab& other)
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
template <typename Trajectory>
void BoundaryMomentumInjectRestrictSlab<Trajectory>::SetupBoundary(bool construct)
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
template <typename Trajectory>
void BoundaryMomentumInjectRestrictSlab<Trajectory>::EvaluateBoundary(void)
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
template <typename Trajectory>
BoundaryMomentumInjectRestrictShell<Trajectory>::BoundaryMomentumInjectRestrictShell(void)
                                   : BoundaryMomentumInject(bnd_name_momentum_inject_restrict_shell, 0, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
};

/*!
\author Vladimir Florinski
\date 06/27/2025
\param[in] other Object to initialize from
*/
template <typename Trajectory>
BoundaryMomentumInjectRestrictShell<Trajectory>::BoundaryMomentumInjectRestrictShell(const BoundaryMomentumInjectRestrictShell& other)
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
template <typename Trajectory>
void BoundaryMomentumInjectRestrictShell<Trajectory>::SetupBoundary(bool construct)
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
template <typename Trajectory>
void BoundaryMomentumInjectRestrictShell<Trajectory>::EvaluateBoundary(void)
{
   BoundaryMomentumInject::EvaluateBoundary();

// If the momentum boundary crossing occurred outside the slab, "_delta_old" is reset to have the same sign as "_delta" to avoid triggering the crossing event.
   if (_delta * _delta_old < 0.0) {
      auto rdist = (_pos - r0).Norm(); 
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
template <typename Trajectory>
BoundaryMirror<Trajectory>::BoundaryMirror(void)
              : BoundaryBase(bnd_name_mirror, 0, BOUNDARY_MOMENTUM | BOUNDARY_REFLECT)
{
   max_crossings = -1;
};

/*!
\author Vladimir Florinski
\date 03/07/2022
\param[in] other Object to initialize from
*/
template <typename Trajectory>
BoundaryMirror<Trajectory>::BoundaryMirror(const BoundaryMirror& other)
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
template <typename Trajectory>
void BoundaryMirror<Trajectory>::SetupBoundary(bool construct)
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
template <typename Trajectory>
void BoundaryMirror<Trajectory>::EvaluateBoundary(void)
{
// Delta is the parallel momentum component
   if constexpr (std::derived_from<Trajectory, TrajectoryGuidingBase<Trajectory, Fields>>) {
      _delta = _mom[2];
   }
   else if constexpr (std::same_as<Trajectory, TrajectoryFocused<Fields>>) {
      _delta = _mom[0] * _mom[1];
   }
   else if constexpr (std::same_as<Trajectory, TrajectoryLorentz<Fields>>) {
      _delta = _mom * _fields.AbsMag();
   }
};


};
