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


namespace Spectrum {

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
template <typename Trajectory_>
class InitialMomentumFixed : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;

   static_assert(std::same_as<Trajectory, TrajectoryLorentz<Fields>> || std::derived_from<Trajectory, TrajectoryFieldlineBase<Trajectory, HConfig>>, "InitialMomentumFixed initial type cannot be applied to the selected Trajectory type.");

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
template <typename Trajectory_>
class InitialMomentumBeam : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;
   // methods:
   using InitialBase::GetMomSample;

   static_assert(!std::same_as<Trajectory, TrajectoryParker<Fields>> && !std::derived_from<Trajectory, TrajectoryFieldlineBase<Trajectory, HConfig>>, "InitialMomentumBeam initial type cannot be applied to the selected Trajectory type.");

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
template <typename Trajectory_>
class InitialMomentumRing : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;
   using InitialBase::rng;

   static_assert(!std::same_as<Trajectory, TrajectoryParker<Fields>> && !std::derived_from<Trajectory, TrajectoryFieldlineBase<Trajectory, HConfig>>, "InitialMomentumRing initial type cannot be applied to the selected Trajectory type.");

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
template <typename Trajectory_>
class InitialMomentumShell : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;
   using InitialBase::rng;

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
template <typename Trajectory_>
class InitialMomentumThickShell : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;
   using InitialBase::rng;

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
template <typename Trajectory_>
class InitialMomentumTable : public InitialTable<Trajectory_, GeoVector> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;
   using InitialTable = InitialTable<Trajectory, GeoVector>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;
   using InitialBase::rng;

   using InitialTable::table_counter;
   using InitialTable::initquant;
   using InitialTable::random;
   // methods:
   using InitialTable::SetupInitial;

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
template <typename Trajectory_>
class InitialMomentumMaxwell : public InitialBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using InitialBase = InitialBase<Trajectory>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_mom;
   using InitialBase::axis;
   using InitialBase::rng;

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

// Something like this is needed for templated classes
#include "initial_momentum.cc"

#endif
