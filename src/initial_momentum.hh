/*!
\file initial_momentum.hh
\brief Declares several classes to specify momentum initial conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INITIAL_MOMENTUM_HH
#define SPECTRUM_INITIAL_MOMENTUM_HH

#include "initial_base.hh"

#ifndef TRAJ_TYPE
#error Trajectory type is undefined!
#endif

namespace Spectrum {

#if (TRAJ_TYPE == TRAJ_LORENTZ) || (TRAJ_TYPE == TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumFixed class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumFixed class
const std::string init_name_momentum_fixed = "InitialMomentumFixed";

//! Flag to inidicate coordinates for initial momentum: 0 = Cartesian, 1 = spherical
#define INITIAL_MOM_FIXED_COORD 0

/*!
\brief Fixed initial momentum vector (relative to a preferred direction)
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector initmom, or double p0, double theta0, double phi0
*/
class InitialMomentumFixed : public InitialBase {

protected:

#if INITIAL_MOM_FIXED_COORD == 0
//! Momentum (persistent)
   GeoVector initmom;
#else
//! Magnitude of momentum (persistent)
   double p0;

//! Pitch angle cosine (persistent)
   double mu0;

//! Sine of the pitch angle (persistent)
   double st0;

//! Sine of the azimuthal angle (persistent)
   double sp0;

//! Cosine of the azimuthal angle (persistent)
   double cp0;
#endif

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumFixed(void);

//! Copy constructor
   InitialMomentumFixed(const InitialMomentumFixed& other);

//! Destructor
   ~InitialMomentumFixed() override = default;

//! Clone function
   CloneFunctionInitial(InitialMomentumFixed);
};

#endif

#if (TRAJ_TYPE != TRAJ_PARKER) && (TRAJ_TYPE != TRAJ_PARKER_SOURCE) && (TRAJ_TYPE != TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumBeam class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumBeam class
const std::string init_name_momentum_beam = "InitialMomentumBeam";

/*!
\brief Cold beam initial distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p0
*/
class InitialMomentumBeam : public InitialBase {

protected:

//! Magnitude of momentum (persistent)
   double p0;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumBeam(void);

//! Copy constructor
   InitialMomentumBeam(const InitialMomentumBeam& other);

//! Destructor
   ~InitialMomentumBeam() override = default;

//! Clone function
   CloneFunctionInitial(InitialMomentumBeam);
};

#endif

#if (TRAJ_TYPE != TRAJ_PARKER) && (TRAJ_TYPE != TRAJ_PARKER_SOURCE) && (TRAJ_TYPE != TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumRing class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumRing class
const std::string init_name_momentum_ring = "InitialMomentumRing";

/*!
\brief Ring initial distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p0, double theta0
*/
class InitialMomentumRing : public InitialBase {

protected:

//! Magnitude of momentum (persistent)
   double p0;

//! Pitch angle cosine (persistent)
   double mu0;

//! Sine of the pitch angle (persistent)
   double st0;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumRing(void);

//! Copy constructor
   InitialMomentumRing(const InitialMomentumRing& other);

//! Destructor
   ~InitialMomentumRing() override = default;

//! Clone function
   CloneFunctionInitial(InitialMomentumRing);
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumShell class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumShell class
const std::string init_name_momentum_shell = "InitialMomentumShell";

/*!
\brief Isotropic non-drifting shell distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p0
*/
class InitialMomentumShell : public InitialBase {

protected:

//! Radius of the shell (persistent)
   double p0;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumShell(void);

//! Copy constructor
   InitialMomentumShell(const InitialMomentumShell& other);

//! Destructor
   ~InitialMomentumShell() override = default;

//! Clone function
   CloneFunctionInitial(InitialMomentumShell);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumThickShell class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumThickShell class
const std::string init_name_momentum_thickshell = "InitialMomentumThickShell";

/*!
\brief Uniform and isotropic thick shell distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p1, double p2, bool log_bias
*/
class InitialMomentumThickShell : public InitialBase {

protected:

//! Smallest momentum (persistent)
   double p1;

//! Largest momentum (persistent)
   double p2;

//! Use logarithmic bias (persistent)
   bool log_bias;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumThickShell(void);

//! Copy constructor
   InitialMomentumThickShell(const InitialMomentumThickShell& other);

//! Destructor
   ~InitialMomentumThickShell() override = default;

//! Clone function
   CloneFunctionInitial(InitialMomentumThickShell);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumTable class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumTable class
const std::string init_name_momentum_table = "InitialMomentumTable";

/*!
\brief Starting points from a table
\author Juan G Alonso Guzman

Parameters: (InitialTable)
*/
class InitialMomentumTable : public InitialTable<GeoVector> {

protected:

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumTable(void);

//! Copy constructor
   InitialMomentumTable(const InitialMomentumTable& other);

//! Clone function
   CloneFunctionInitial(InitialMomentumTable);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumMaxwell class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the InitialMomentumMaxwell class
const std::string init_name_momentum_maxwell = "InitialMomentumMaxwell";

/*!
\brief Non-relativistic bi-Maxwellian
\author Vladimir Florinski

Parameters: (InitialBase), double p0, double dp_para, double dp_perp
*/
class InitialMomentumMaxwell : public InitialBase {

protected:

//! Drift momentum (persistent)
   double p0;

//! Parallel spread (persistent)
   double dp_para;

//! Parallel spread (persistent)
   double dp_perp;

//! Set up the initial condition generator based on "params"
   void SetupInitial(bool construct) override;

//! Compute the internal position or momentum
   void EvaluateInitial(void) override;

public:

//! Default constructor
   InitialMomentumMaxwell(void);

//! Copy constructor
   InitialMomentumMaxwell(const InitialMomentumMaxwell& other);

//! Destructor
   ~InitialMomentumMaxwell() override = default;

//! Clone function
   CloneFunctionInitial(InitialMomentumMaxwell);
};

};

#endif
