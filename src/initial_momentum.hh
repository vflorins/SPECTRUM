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

/*!
\brief Fixed initial momentum vector (relative to a preferred direction)
\author Vladimir Florinski

Parameters: (InitialBase), GeoVector initmom, or double p0, double theta0, double phi0
*/
template <typename HConfig_>
class InitialMomentumFixed : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumFixed";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialBase::axis;

   static_assert(HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz || HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline, "InitialMomentumFixed initial type cannot be applied to the selected Trajectory type.");

protected:

//! Momentum (persistent)
   GeoVector initmom;

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

/*!
\brief Cold beam initial distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p0
*/
template <typename HConfig_>
class InitialMomentumBeam : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumBeam";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialBase::axis;

   // methods:
   using InitialBase::GetMomSample;

   static_assert(!(HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Parker) && !(HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline), "InitialMomentumRing initial type cannot be applied to the selected Trajectory type.");

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

/*!
\brief Ring initial distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p0, double theta0
*/
template <typename HConfig_>
class InitialMomentumRing : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumRing";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialBase::axis;

   using InitialBase::rng;

   static_assert(!(HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Parker) && !(HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline), "InitialMomentumRing initial type cannot be applied to the selected Trajectory type.");

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

/*!
\brief Isotropic non-drifting shell distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p0
*/
template <typename HConfig_>
class InitialMomentumShell : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumShell";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
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

/*!
\brief Uniform and isotropic thick shell distribution
\author Vladimir Florinski

Parameters: (InitialBase), double p1, double p2, bool log_bias
*/
template <typename HConfig_>
class InitialMomentumThickShell : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumThickShell";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
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

/*!
\brief Starting points from a table
\author Juan G Alonso Guzman

Parameters: (InitialTable)
*/
template <typename HConfig_>
class InitialMomentumTable : public InitialTable<HConfig_, GeoVector> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumTable";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
   using InitialBase::axis;
   using InitialTable = InitialTable<HConfig, GeoVector>;

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

/*!
\brief Non-relativistic bi-Maxwellian
\author Vladimir Florinski

Parameters: (InitialBase), double p0, double dp_para, double dp_perp
*/
template <typename HConfig_>
class InitialMomentumMaxwell : public InitialBase<HConfig_> {
private:

   //! Readable name of the class
   static constexpr std::string_view init_name = "InitialMomentumMaxwell";

public:

   using HConfig = HConfig_;
   using InitialBase = InitialBase<HConfig>;

   using InitialBase::_status;
   using InitialBase::container;
   using InitialBase::_coords;
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
