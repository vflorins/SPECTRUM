/*!
\file boundary_momentum.hh
\brief Declares several classes representing momentum boundaries
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BOUNDARY_MOMENTUM_HH
#define SPECTRUM_BOUNDARY_MOMENTUM_HH

#include "boundary_base.hh"

#ifndef TRAJ_TYPE
#error Trajectory type is undefined!
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentum class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryMomentum class
const std::string bnd_name_momentum = "BoundaryMomentum";

/*!
\brief Energy boundary (e.g. injection)
\author Vladimir Florinski

Parameters: (BoundaryBase), double momentum
*/
class BoundaryMomentum : public BoundaryBase {

protected:

//! Momentum value
   double momentum;

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Default constructor
   BoundaryMomentum(void);

//! Copy constructor
   BoundaryMomentum(const BoundaryMomentum& other);

//! Destructor
   ~BoundaryMomentum() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryMomentum);
};

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMirror class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryMirror class
const std::string bnd_name_mirror = "BoundaryMirror";

/*!
\brief Mirroring event counter.
\author Vladimir Florinski

Parameters: (BoundaryBase)
*/
class BoundaryMirror : public BoundaryBase {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Default constructor
   BoundaryMirror(void);

//! Copy constructor
   BoundaryMirror(const BoundaryMirror& other);

//! Destructor
   ~BoundaryMirror() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryMirror);
};

#endif

};

#endif
