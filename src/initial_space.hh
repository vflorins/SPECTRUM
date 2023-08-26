/*!
\file initial_space.hh
\brief Declares several classes to specify spatial initial conditions
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INITIAL_SPACE_HH
#define SPECTRUM_INITIAL_SPACE_HH

#include "initial_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceFixed class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialSpaceFixed class
const std::string init_name_space_fixed = "InitialSpaceFixed";

/*!
\brief Starting points at a fixed position
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector initpos
*/
class InitialSpaceFixed : public InitialBase {

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

//! Readable name of the InitialSpaceLine class
const std::string init_name_space_line = "InitialSpaceLine";

/*!
\brief Uniformly distributed starting points on a line segment
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector startpos, GeoVector endpos, int n_intervals
*/
class InitialSpaceLine : public InitialBase {

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
   InitialSpaceLine(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

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
// InitialSpaceCube class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialSpaceCube class
const std::string init_name_space_cube = "InitialSpaceCube";

/*!
\brief Starting points uniformly distributed over the volume of a cube
\author Vladimir Florinski

Parameters: (InitialSpaceLine)
*/
class InitialSpaceCube : public InitialSpaceLine {

protected:

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialSpaceCube(void);

//! Copy constructor
   InitialSpaceCube(const InitialSpaceCube& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceCube);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphere class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialSpaceSphere class
const std::string init_name_space_sphere = "InitialSpaceSphere";

/*!
\brief Uniformly distributed starting points on a sphere
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector origin, double radius
*/
class InitialSpaceSphere : public InitialBase {

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
   InitialSpaceSphere(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   InitialSpaceSphere(const InitialSpaceSphere& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceSphere);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphereSector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialSpaceSphereSector class
const std::string init_name_space_sphere_sector = "InitialSpaceSphereSector";

/*!
\brief Uniformly distributed starting points on a sphere
\author Vladimir Florinski

Parameters: (InitialSpaceSphere), theta1, theta2, phi1, phi2
*/
class InitialSpaceSphereSector : public InitialSpaceSphere {

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

//! Readable name of the InitialSpaceRankineHalfBody class
const std::string init_name_space_rankine_half_body = "InitialSpaceRankineHalfBody";

/*!
\brief Starting points on Rankine Half Body
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (InitialBase), GeoVector origin, GeoVector axis, double z_nose, double radius
*/
class InitialSpaceRankineHalfBody : public InitialBase {

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

//! Readable name of the InitialSpaceTable class
const std::string init_name_space_table = "InitialSpaceTable";

/*!
\brief Starting points from a table
\author Juan G Alonso Guzman

Parameters: (InitialBase), std::string
*/
class InitialSpaceTable : public InitialBase {

protected:

//! Array with initial positions
   std::vector <GeoVector> initpos;

//! Table entry counter
   int table_counter;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

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

//! Readable name of the InitialSpaceCylinder class
const std::string init_name_space_cylinder = "InitialSpaceCylinder";

/*!
\brief Uniformly distributed starting points on a cylinder (just the curved portion, not the top and bottom circular surfaces)
\author Juan G Alonso Guzman

Parameters: (InitialBase), GeoVector origin, GeoVector height, GeoVector radius_x
*/
class InitialSpaceCylinder : public InitialBase {

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
   InitialSpaceCylinder(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   InitialSpaceCylinder(const InitialSpaceCylinder& other);

//! Clone function
   CloneFunctionInitial(InitialSpaceCylinder);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinderSector class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialSpaceCylinderSector class
const std::string init_name_space_cylinder_sector = "InitialSpaceCylinderSector";

/*!
\brief Uniformly distributed starting points on a cylindrical sector (just the curved portion, not the top and bottom wedge surfaces)
\author Juan G Alonso Guzman

Parameters: (InitialSpaceCylinder), phi1, phi2
*/
class InitialSpaceCylinderSector : public InitialSpaceCylinder {

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

#endif
