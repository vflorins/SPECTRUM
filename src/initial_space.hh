/*!
\file initial_space.hh
\brief Declares several classes to specify spatial initial conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INITIAL_SPACE_HH
#define SPECTRUM_INITIAL_SPACE_HH

#include "initial_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceFixed class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Starting points at a fixed position
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector initpos
*/
template <typename HConfig_>
class InitialSpaceFixed : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceFixed";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;

protected:

//! Position (persistent)
   GeoVector initpos;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceFixed(void);

//! Copy constructor
   InitialSpaceFixed(const InitialSpaceFixed& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceFixed);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceLine class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Uniformly distributed starting points on a line segment
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector startpos, GeoVector endpos, int n_intervals
*/
template <typename HConfig_>
class InitialSpaceLine : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceLine";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;

   using InitialBase::rng;

protected:

//! Random or evenly distributed
   bool randompos;

//! First point (persistent)
   GeoVector startpos;

//! Second point (persistent)
   GeoVector endpos;

//! Increment (persistent)
   GeoVector increment;

//! Constructor with arguments (to speed up construction of derived classes)
   InitialSpaceLine(const std::string& name_in, uint16_t status_in);

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceLine(void);

//! Copy constructor
   InitialSpaceLine(const InitialSpaceLine& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceLine);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCircle class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Uniformly distributed starting points on a circle
\author Juan G Alonso Guzman

Parameters: (InitialBase), GeoVector center, GeoVector normal, double radius
*/
template <typename HConfig_>
class InitialSpaceCircle : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceCircle";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;

   using InitialBase::rng;

protected:

//! Center of circle (persistent)
   GeoVector origin;

//! Radius vector of circle in the local x direction (persistent)
   GeoVector radius_x;

//! Radius vector of circle in the local y direction (persistent)
   GeoVector radius_y;

//! Constructor with arguments (to speed up construction of derived classes)
   InitialSpaceCircle(const std::string& name_in, uint16_t status_in);

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceCircle(void);

//! Copy constructor
   InitialSpaceCircle(const InitialSpaceCircle& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceCircle);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceBox class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Starting points uniformly distributed over the volume of a box
\author Vladimir Florinski

Parameters: (InitialSpaceLine)
*/
template <typename HConfig_>
class InitialSpaceBox : public InitialSpaceLine<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceBox";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialSpaceLine = InitialSpaceLine<HConfig>;

   using InitialBase::rng;

   using InitialSpaceLine::startpos;
   using InitialSpaceLine::endpos;
   using InitialSpaceLine::SetupInitial;


protected:

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceBox(void);

//! Copy constructor
   InitialSpaceBox(const InitialSpaceBox& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceBox);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphere class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Uniformly distributed starting points on a sphere
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector origin, double radius
*/
template <typename HConfig_>
class InitialSpaceSphere : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceSphere";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;

   using InitialBase::rng;

protected:

//! Origin (persistent)
   GeoVector origin;

//! Radius of the sphere (persistent)
   double radius;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceSphere(void);

//! Constructor with arguments (to speed up construction of derived classes)
   InitialSpaceSphere(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   InitialSpaceSphere(const InitialSpaceSphere& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceSphere);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphereSector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Uniformly distributed starting points on a sphere
\author Vladimir Florinski

Parameters: (InitialSpaceSphere), theta1, theta2, phi1, phi2
*/
template <typename HConfig_>
class InitialSpaceSphereSector : public InitialSpaceSphere<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceSphereSector";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialSpaceSphere = InitialSpaceSphere<HConfig>;

   using InitialBase::rng;

   using InitialSpaceSphere::origin;
   using InitialSpaceSphere::radius;

protected:

//! Theta limits for sector
   double theta1, theta2;

//! Cosine of theta limits for sector
   double costheta1, costheta2;

//! Phi limits for sector
   double phi1, phi2;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceSphereSector(void);

//! Copy constructor
   InitialSpaceSphereSector(const InitialSpaceSphereSector& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceSphereSector);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceRankineHalfBody class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Starting points on Rankine Half Body
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (InitialBase), GeoVector origin, GeoVector axis, double z_nose, double radius
*/
template <typename HConfig_>
class InitialSpaceRankineHalfBody : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceRankineHalfBody";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;

   using InitialBase::rng;

protected:

//! Origin
   GeoVector origin;

//! Symmetry axis
   GeoVector axis;

//! Distance to the nose
   double z_nose;

//! Limiting radius
   double radius;

//! Limiting cos(theta)
   double cos_lim;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceRankineHalfBody(void);

//! Copy constructor
   InitialSpaceRankineHalfBody(const InitialSpaceRankineHalfBody& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceRankineHalfBody);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceTable class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Starting points from a table
\author Juan G Alonso Guzman

Parameters: (InitialTable)
*/
template <typename HConfig_>
class InitialSpaceTable : public InitialTable<HConfig_, GeoVector> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceTable";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialTable = InitialTable<HConfig, GeoVector>;

   using InitialBase::rng;

   using InitialTable::random;
   using InitialTable::table_counter;
   using InitialTable::initquant;
   using InitialTable::SetupInitial;

protected:

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceTable(void);

//! Copy constructor
   InitialSpaceTable(const InitialSpaceTable& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceTable);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinder class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Uniformly distributed starting points on a cylinder (just the curved portion, not the top and bottom circular surfaces)
\author Juan G Alonso Guzman

Parameters: (InitialBase), GeoVector origin, GeoVector height, GeoVector radius_x
*/
template <typename HConfig_>
class InitialSpaceCylinder : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceCylinder";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;

   using InitialBase::rng;

protected:

//! Origin (persistent)
   GeoVector origin;

//! Height (persistent)
   GeoVector height;

//! Radius vector of cylinder in the local x direction (persistent)
   GeoVector radius_x;

//! Radius vector of cylinder in the local y direction (persistent)
   GeoVector radius_y;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceCylinder(void);

//! Constructor with arguments (to speed up construction of derived classes)
   InitialSpaceCylinder(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   InitialSpaceCylinder(const InitialSpaceCylinder& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceCylinder);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinderSector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Uniformly distributed starting points on a cylindrical sector (just the curved portion, not the top and bottom wedge surfaces)
\author Juan G Alonso Guzman

Parameters: (InitialSpaceCylinder), phi1, phi2
*/
template <typename HConfig_>
class InitialSpaceCylinderSector : public InitialSpaceCylinder<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialSpaceCylinderSector";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;
   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialSpaceCylinder = InitialSpaceCylinder<HConfig>;

   using InitialBase::rng;

   using InitialSpaceCylinder::origin;
   using InitialSpaceCylinder::radius_x;
   using InitialSpaceCylinder::radius_y;
   using InitialSpaceCylinder::height;

protected:

//! Phi limits for sector
   double phi1, phi2;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceCylinderSector(void);

//! Copy constructor
   InitialSpaceCylinderSector(const InitialSpaceCylinderSector& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceCylinderSector);
};

};

// Something like this is needed for templated classes
#include "initial_space.cc"

#endif


