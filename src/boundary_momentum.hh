/*!
\file boundary_momentum.hh
\brief Declares several classes representing momentum boundaries
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BOUNDARY_MOMENTUM_HH
#define SPECTRUM_BOUNDARY_MOMENTUM_HH

#include "boundary_base.hh"


namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentum class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Momentum boundary
\author Vladimir Florinski

Parameters: (BoundaryBase), double momentum
*/
template <typename HConfig_>
class BoundaryMomentum : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryMomentum";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryBase::_delta;
   using BoundaryBase::max_crossings;

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
   BoundaryMomentum(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BoundaryMomentum(const BoundaryMomentum& other);

public:

//! Destructor
   ~BoundaryMomentum() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMomentumInject class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Absorbing (injection) momentum boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryMomentum)
*/
template <typename HConfig_>
class BoundaryMomentumInject : public BoundaryMomentum<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryMomentumInject";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryMomentum = BoundaryMomentum<HConfig>;
   using BoundaryMomentum::_status;
   using BoundaryMomentum::container;
   using BoundaryMomentum::_coords;
   using BoundaryMomentum::_fields;
   using BoundaryMomentum::max_crossings;

protected:

//! Set up the boundary evaluator based on "params"
   void SetupBoundary(bool construct) override;

public:

//! Default constructor
   BoundaryMomentumInject(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BoundaryMomentumInject(const std::string& name_in, uint16_t status_in);

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

/*!
\brief Passive momentum boundary
\author Juan G Alonso Guzman

Parameters: (BoundaryMomentum)
*/
template <typename HConfig_>
class BoundaryMomentumPass : public BoundaryMomentum<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryMomentumPass";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryMomentum = BoundaryMomentum<HConfig>;

   using BoundaryBase::max_crossings;

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

/*!
\brief Absorbing (injection) momentum boundary restricted to space between two parallel planes
\author Vladimir Florinski

Parameters: (BoundaryMomentumInject), GeoVector r0, GeoVector r1, GeoVector normal 
*/
template <typename HConfig_>
class BoundaryMomentumInjectRestrictSlab : public BoundaryMomentumInject<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryMomentumInjectRestrictSlab";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;
   using BoundaryMomentumInject = BoundaryMomentumInject<HConfig>;

   using BoundaryBase::_delta;
   using BoundaryBase::_delta_old;

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

/*!
\brief Absorbing (injection) momentum boundary restricted to space between two concentric spheres
\author Vladimir Florinski

Parameters: (BoundaryMomentumInject), GeoVector r0, double r1, double r2 
*/
template <typename HConfig_>
class BoundaryMomentumInjectRestrictShell : public BoundaryMomentumInject<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryMomentumInjectRestrictShell";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryMomentumInject = BoundaryMomentumInject<HConfig>;
   using BoundaryMomentumInject::_status;
   using BoundaryMomentumInject::container;
   using BoundaryMomentumInject::_coords;
   using BoundaryMomentumInject::_fields;
   using BoundaryMomentumInject::_delta;
   using BoundaryMomentumInject::_delta_old;

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


//----------------------------------------------------------------------------------------------------------------------------------------------------
// BoundaryMirror class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Mirroring event counter.
\author Vladimir Florinski

Parameters: (BoundaryBase)
*/
template <typename HConfig_>
class BoundaryMirror : public BoundaryBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bdy_name = "BoundaryMirror";

public:

   using HConfig = HConfig_;
   using BoundaryBase = BoundaryBase<HConfig>;
   using BoundaryBase::_status;
   using BoundaryBase::container;
   using BoundaryBase::_coords;
   using BoundaryBase::_fields;

   using BoundaryBase::max_crossings;
   using BoundaryBase::_delta;

   static_assert(!(HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Parker || HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline), "BoundaryMirror boundary type cannot be applied to the selected Trajectory type.");

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


};

// Something like this is needed for templated classes
#include "boundary_momentum.cc"

#endif
