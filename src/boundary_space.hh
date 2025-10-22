/*!
\file boundary_space.hh
\brief Declares several classes representing spatial boundaries
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BOUNDARY_SPACE_HH
#define SPECTRUM_BOUNDARY_SPACE_HH

#include "boundary_base.hh"
#include "common/fields/generated/field_lists.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryPlane class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plane boundary base class
\author Vladimir Florinski

Parameters: (BoundaryBase), GeoVector origin, GeoVector normal
*/
template <typename HConfig_>
class BoundaryPlane : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryPlane";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::_delta;
   using BoundaryBase::_normal;

protected:

//! Origin (persistent)
   GeoVector origin;

//! Normal vector, not to be confused with "_normal" which could have two possible directions (persistent)
   GeoVector norm;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryPlane(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryPlane(const std::string& name_in, uint16_t status_in);

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

/*!
\brief Absorbing plane boundary
\author Vladimir Florinski

Parameters: (BoundaryPlane)
*/
template <typename HConfig_>
class BoundaryPlaneAbsorb : public BoundaryPlane<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryPlaneAbsorb";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryPlane = BoundaryPlane<HConfig>;

   using BoundaryBase::max_crossings;

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

/*!
\brief Reflecting plane boundary
\author Vladimir Florinski

Parameters: (BoundaryPlane)
*/
template <typename HConfig_>
class BoundaryPlaneReflect : public BoundaryPlane<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryPlaneReflect";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryPlane = BoundaryPlane<HConfig>;

   using BoundaryBase::max_crossings;

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

/*!
\brief Plane crossing recording boundary (event)
\author Vladimir Florinski

Parameters: (BoundaryPlane)
*/
template <typename HConfig_>
class BoundaryPlanePass : public BoundaryPlane<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryPlanePass";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryPlane = BoundaryPlane<HConfig>;

   using BoundaryBase::max_crossings;

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
template <typename HConfig_>
class BoundaryBox : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryBox";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::max_crossings;
   using BoundaryBase::_delta;
   using BoundaryBase::_normal;

protected:

//! Opposite corners bounding the box (persistent)
   GeoVector corners[2];

//! Normal vectors to three non-parallel sides (i.e. three sides with a common vertex) (persistent)
   GeoVector normals[3];

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryBox(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryBox(const std::string& name_in, uint16_t status_in);

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

/*!
\brief Reflecting box boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryBox)
*/
template <typename HConfig_>
class BoundaryBoxReflect : public BoundaryBox<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryBoxReflect";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryBox = BoundaryBox<HConfig>;

   using BoundaryBase::max_crossings;

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
template <typename HConfig_>
class BoundarySphere : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundarySphere";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::_delta;
   using BoundaryBase::_normal;

protected:

//! Origin (persistent)
   GeoVector origin;

//! Radius of the sphere (persistent)
   double radius;

//! Default constructor (protected, class not designed to be instantiated)
   BoundarySphere(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundarySphere(const std::string& name_in, uint16_t status_in);

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

/*!
\brief Absorbing spherical boundary
\author Vladimir Florinski

Parameters: (BoundarySphere)
*/
template <typename HConfig_>
class BoundarySphereAbsorb : public BoundarySphere<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundarySphereAbsorb";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundarySphere = BoundarySphere<HConfig>;

   using BoundaryBase::max_crossings;

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

/*!
\brief Reflecting spherical boundary
\author Vladimir Florinski

Parameters: (BoundarySphere)
*/
template <typename HConfig_>
class BoundarySphereReflect : public BoundarySphere<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundarySphereReflect";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundarySphere = BoundarySphere<HConfig>;

   using BoundaryBase::max_crossings;

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
// BoundaryRankine class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Boundary in the shape of the Rankine half-body
\author Juan G Alonso Guzman

Parameters: (BoundaryBase), GeoVector origin, double radius
*/
template <typename HConfig_>
class BoundaryRankine : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryRankine";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::_delta;

protected:

//! Origin (persistent)
   GeoVector origin;

//! Symmetry axis (persistent)
   GeoVector axis;

//! Distance to the nose (persistent)
   double z_nose;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryRankine(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryRankine(const std::string& name_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryRankine(const BoundaryRankine& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundaryRankine() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRankineAbsorb class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Absorbing boundary in the shape of the Rankine half-body
\author Juan G Alonso Guzman

Parameters: (BoundaryRankine)
*/
template <typename HConfig_>
class BoundaryRankineAbsorb : public BoundaryRankine<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryRankineAbsorb";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryRankine = BoundaryRankine<HConfig>;

   using BoundaryBase::max_crossings;

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

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
template <typename HConfig_>
class BoundaryCylinder : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryCylinder";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::_delta;
   using BoundaryBase::_normal;

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
   BoundaryCylinder(const std::string& name_in, uint16_t status_in);

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

/*!
\brief Absorbing cylindrical boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryCylinder)
*/
template <typename HConfig_>
class BoundaryCylinderAbsorb : public BoundaryCylinder<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryCylinderAbsorb";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryCylinder = BoundaryCylinder<HConfig>;

   using BoundaryBase::max_crossings;

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

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRegion class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Region boundary base class
\author Juan G Alonso Guzman

Parameters: (BoundaryBase), int region_ind, double region_val
*/
// TODO: allow the region indicator variable to vary
template <typename HConfig_>
class BoundaryRegion : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryRegion";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::_delta;
   using BoundaryBase::_normal;

protected:

//! Region index (persistent)
//   int region_ind;

//! Region threshold value (persistent)
   double region_val;

//! Default constructor (protected, class not designed to be instantiated)
   BoundaryRegion(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryRegion(const std::string& name_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BoundaryRegion(const BoundaryRegion& other);

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

//! Compute the distance to the boundary
   void EvaluateBoundary(void) override;

public:

//! Destructor
   ~BoundaryRegion() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryRegionAbsorb class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Absorbing Region boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryRegion)
*/
template <typename HConfig_>
class BoundaryRegionAbsorb : public BoundaryRegion<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryRegionAbsorb";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryRegion = BoundaryRegion<HConfig>;

   using BoundaryBase::max_crossings;

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryRegionAbsorb(void);

//! Copy constructor
   BoundaryRegionAbsorb(const BoundaryRegionAbsorb& other);

//! Destructor
   ~BoundaryRegionAbsorb() override = default;

//! Clone function
   CloneFunctionBoundary(BoundaryRegionAbsorb);
};

};

// Something like this is needed for templated classes
#include "boundary_space.cc"

#endif
