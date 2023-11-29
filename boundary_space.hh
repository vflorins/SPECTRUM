/*!
\file boundary_space.hh
\brief Declares several classes representing spatial boundaries
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BOUNDARY_SPACE_HH
#define SPECTRUM_BOUNDARY_SPACE_HH

#include "boundary_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlane class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plane boundary base class
\author Vladimir Florinski

Parameters: (BoundaryBase), GeoVector origin, GeoVector normal
*/
class BoundaryPlane : public BoundaryBase {

protected:

//! Origin (persistent)
   GeoVector origin;

//! Normal vector, not to be confused with "_normal" which could have two possible directions (persistent)
   GeoVector norm;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryPlane(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryPlane(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryPlane(const BoundaryPlane& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundaryPlane() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlaneAbsorb class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryPlaneAbsorb class
const std::string bnd_name_plane_absorb = "BoundaryPlaneAbsorb";

/*!
\brief Absorbing plane boundary
\author Vladimir Florinski

Parameters: (BoundaryPlane)
*/
class BoundaryPlaneAbsorb : public BoundaryPlane {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryPlaneAbsorb(void);

//! Copy constructor
   BoundaryPlaneAbsorb(const BoundaryPlaneAbsorb& other);

//! Destructor
   ~BoundaryPlaneAbsorb() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryPlaneAbsorb);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlaneReflect class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryPlaneReflect class
const std::string bnd_name_plane_reflect = "BoundaryPlaneReflect";

/*!
\brief Reflecting plane boundary
\author Vladimir Florinski

Parameters: (BoundaryPlane)
*/
class BoundaryPlaneReflect : public BoundaryPlane {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryPlaneReflect(void);

//! Copy constructor
   BoundaryPlaneReflect(const BoundaryPlaneReflect& other);

//! Destructor
   ~BoundaryPlaneReflect() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryPlaneReflect);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlanePass class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryPlanePass class
const std::string bnd_name_plane_pass = "BoundaryPlanePass";

/*!
\brief Plane crossing recording boundary (event)
\author Vladimir Florinski

Parameters: (BoundaryPlane), int max_crossings
*/
class BoundaryPlanePass : public BoundaryPlane {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryPlanePass(void);

//! Copy constructor
   BoundaryPlanePass(const BoundaryPlanePass& other);

//! Destructor
   ~BoundaryPlanePass() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryPlanePass);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBox class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Box boundary base class
\author Juan G Alonso Guzman

Parameters: (BoundaryBase), GeoVector corners[0], GeoVector[3] normals
*/
class BoundaryBox : public BoundaryBase {

protected:

//! Opposite corners bounding the box (persistent)
   GeoVector corners[2];

//! Normal vectors to three non-parallel sides (i.e. three sides with a common vertex) (persistent)
   GeoVector normals[3];

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryBox(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryBox(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryBox(const BoundaryBox& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundaryBox() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryBoxReflect class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryBoxReflect class
const std::string bnd_name_box_reflect = "BoundaryBoxReflect";

/*!
\brief Reflecting box boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryBox)
*/
class BoundaryBoxReflect : public BoundaryBox {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryBoxReflect(void);

//! Copy constructor
   BoundaryBoxReflect(const BoundaryBoxReflect& other);

//! Destructor
   ~BoundaryBoxReflect() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryBoxReflect);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphere class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Spherical boundary base class
\author Vladimir Florinski

Parameters: (BoundaryBase), GeoVector origin, double radius
*/
class BoundarySphere : public BoundaryBase {

protected:

//! Origin (persistent)
   GeoVector origin;

//! Radius of the sphere (persistent)
   double radius;

//! Default constructor (protected, class not designed to be instantiated)
   BoundarySphere(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundarySphere(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundarySphere(const BoundarySphere& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundarySphere() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphereAbsorb class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundarySphereAbsorb class
const std::string bnd_name_sphere_absorb = "BoundarySphereAbsorb";

/*!
\brief Absorbing spherical boundary
\author Vladimir Florinski

Parameters: (BoundarySphere)
*/
class BoundarySphereAbsorb : public BoundarySphere {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundarySphereAbsorb(void);

//! Copy constructor
   BoundarySphereAbsorb(const BoundarySphereAbsorb& other);

//! Destructor
   ~BoundarySphereAbsorb() override = default;

//! Clone function
   CloneFunctionBoundary(BoundarySphereAbsorb);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundarySphereReflect class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundarySphereReflect class
const std::string bnd_name_sphere_reflect = "BoundarySphereReflect";

/*!
\brief Reflecting spherical boundary
\author Vladimir Florinski

Parameters: (BoundarySphere)
*/
class BoundarySphereReflect : public BoundarySphere {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundarySphereReflect(void);

//! Copy constructor
   BoundarySphereReflect(const BoundarySphereReflect& other);

//! Destructor
   ~BoundarySphereReflect() override = default;

//! Clone function
   CloneFunctionBoundary(BoundarySphereReflect);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRankineAbsorb class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryRankineAbsorb class
const std::string bnd_name_rankine_absorb = "BoundaryRankineAbsorb";

/*!
\brief Absorbing boundary in the shape of the Rankine half-body
\author Vladimir Florinski

Parameters: (BoundaryBase), GeoVector origin, GeoVector axis, double z_nose
*/
class BoundaryRankineAbsorb : public BoundaryBase {

protected:

//! Origin
   GeoVector origin;

//! Symmetry axis
   GeoVector axis;

//! Distance to the nose
   double z_nose;

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Default constructor
   BoundaryRankineAbsorb(void);

//! Copy constructor
   BoundaryRankineAbsorb(const BoundaryRankineAbsorb& other);

//! Destructor
   ~BoundaryRankineAbsorb() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryRankineAbsorb);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryCylinder class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Cylindrical boundary base class
\author Juan G Alonso Guzman

Parameters: (BoundaryBase), GeoVector origin, GeoVector fa_basis[2], double radius
*/
class BoundaryCylinder : public BoundaryBase {

protected:

//! Point on cylinder axis (persisitent)
   GeoVector origin;

//! Cylinder coordinate reference frame with z || fa_basis[2] (persitent)
   GeoVector fa_basis[3];

//! Radius of the cylinder (persistent)
   double radius;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryCylinder(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryCylinder(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryCylinder(const BoundaryCylinder& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundaryCylinder() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryCylinderAbsorb class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BoundaryCylinderAbsorb class
const std::string bnd_name_cylinder_absorb = "BoundaryCylinderAbsorb";

/*!
\brief Absorbing cylindrical boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryCylinder)
*/
class BoundaryCylinderAbsorb : public BoundaryCylinder {

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryCylinderAbsorb(void);

//! Copy constructor
   BoundaryCylinderAbsorb(const BoundaryCylinderAbsorb& other);

//! Destructor
   ~BoundaryCylinderAbsorb() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryCylinderAbsorb);
};

};

#endif
