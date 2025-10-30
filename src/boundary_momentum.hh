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

/*!
\brief Momentum boundary
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

//! Default constructor
   BoundaryMomentum(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryMomentum(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BoundaryMomentum(const BoundaryMomentum& other);

public:

//! Destructor
   ~BoundaryMomentum() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInject class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryMomentumInject class
const std::string bnd_name_momentum_inject = "BoundaryMomentumInject";

/*!
\brief Absorbing (injection) momentum boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryMomentum)
*/
class BoundaryMomentumInject : public BoundaryMomentum {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryMomentumInject(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryMomentumInject(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BoundaryMomentumInject(const BoundaryMomentumInject& other);

//! Destructor
   ~BoundaryMomentumInject() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryMomentumInject);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumPass class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryMomentumPass class
const std::string bnd_name_momentum_pass = "BoundaryMomentumPass";

/*!
\brief Passive momentum boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryMomentum)
*/
class BoundaryMomentumPass : public BoundaryMomentum {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryMomentumPass(void);

//! Copy constructor
   BoundaryMomentumPass(const BoundaryMomentumPass& other);

//! Destructor
   ~BoundaryMomentumPass() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryMomentumPass);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInjectRestrictSlab class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryMomentumInjectRestrictSlab class
const std::string bnd_name_momentum_inject_restrict_slab = "BoundaryMomentumInjectRestrictSlab";

/*!
\brief Absorbing (injection) momentum boundary restricted to space between two parallel planes
\author Vladimir Florinski

Parameters: (BoundaryMomentumInject), GeoVector r0, GeoVector r1, GeoVector normal 
*/
class BoundaryMomentumInjectRestrictSlab : public BoundaryMomentumInject {

protected:

//! Position of the leading face of the slab
   GeoVector r0;

//! Position of the trailing face of the slab
   GeoVector r1;

//! Normal to the slab
   GeoVector normal;

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Default constructor
   BoundaryMomentumInjectRestrictSlab(void);

//! Copy constructor
   BoundaryMomentumInjectRestrictSlab(const BoundaryMomentumInjectRestrictSlab& other);

//! Destructor
   ~BoundaryMomentumInjectRestrictSlab() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryMomentumInjectRestrictSlab);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInjectRestrictShell class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryMomentumInjectRestrictShell class
const std::string bnd_name_momentum_inject_restrict_shell = "BoundaryMomentumInjectRestrictShell";

/*!
\brief Absorbing (injection) momentum boundary restricted to space between two concentric spheres
\author Vladimir Florinski

Parameters: (BoundaryMomentumInject), GeoVector r0, double r1, double r2 
*/
class BoundaryMomentumInjectRestrictShell : public BoundaryMomentumInject {

protected:

//! Center of the shperes
   GeoVector r0;

//! Radius of the inner sphere
   double r1;

//! Radius of the outer sphere
   double r2;

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Default constructor
   BoundaryMomentumInjectRestrictShell(void);

//! Copy constructor
   BoundaryMomentumInjectRestrictShell(const BoundaryMomentumInjectRestrictShell& other);

//! Destructor
   ~BoundaryMomentumInjectRestrictShell() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryMomentumInjectRestrictShell);
};

#if (TRAJ_TYPE != TRAJ_PARKER) && (TRAJ_TYPE != TRAJ_PARKER_SOURCE) && (TRAJ_TYPE != TRAJ_FIELDLINE)

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
