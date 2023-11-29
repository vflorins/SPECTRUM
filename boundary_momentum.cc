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
\author Vladimir Florinski
\date 03/25/2022
*/
BoundaryMomentum::BoundaryMomentum(void)
                : BoundaryBase(bnd_name_momentum, 0, BOUNDARY_MOMENTUM | BOUNDARY_TERMINAL)
{
   max_crossings = 1;
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
   max_crossings = 1;
   RAISE_BITS(_status, BOUNDARY_TERMINAL);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
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
   if(!construct) BoundaryBase::SetupBoundary(false);
   container.Read(&momentum);
};

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
void BoundaryMomentum::EvaluateBoundary(void)
{
   _delta = _mom.Norm() - momentum;
};

#if TRAJ_TYPE != TRAJ_PARKER

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
   max_crossings = -1;
   RAISE_BITS(_status, BOUNDARY_REFLECT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBoundary(true);
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
   if(!construct) BoundaryBase::SetupBoundary(false);
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
