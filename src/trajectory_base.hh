/*!
\file trajectory_base.hh
\brief Declares a base class for trajectory tracing
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_BASE_HH
#define SPECTRUM_TRAJECTORY_BASE_HH

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, params, physics, random, spatial_data, vectors

#include "background.hh"
#include "diffusion.hh"

#include "distribution_base.hh"
#include "boundary_base.hh"
#include "initial_base.hh"
#include "common/rk_lists.hh"
#include "utils_records.hh"

#include "common/fields2/field_ops.hh"


namespace Spectrum {

//! Clone function pattern
#define CloneFunctionTrajectory(T) std::unique_ptr<TrajectoryBase> Clone(void) const override {return std::make_unique<T>();};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A base class to integrate a trajectory
\author Vladimir Florinski
\author Lucius Schoenbaum

A trajectory should be thought of as a self-contained simulation. There are four ingredients to a trajectory: (a) the transport physics that is built into the class itself, (b) the background u, E, B fields provided through "background", (c) the boundary conditions contained in "bcond_t", "bcond_s", and "bcond_m", and (d) the initial conditions provided by "icond_s" and "icond_m". The trajectory is responsible for computing the derived transport coefficients. This is done for efficiency purposes because a separate transport class hierarchy would have to interact with the other components and exchanging different kinds of transport parameters must be done through a container, hence require extra load/store operations.

A trajectory object is considered initialized if (a) background is assigned, (b) at least one time boundary is assigned, (c) the space initial condition is assigned, and (d) the momentum initial condition is assigned.
*/
template <typename HConfig_>
class TrajectoryBase : public Params {

public:

   using HConfig = HConfig_;
   using Config = HConfig::TrajectoryConfig;
   using Coordinates = Config::Coordinates;
   using RecordCoordinates = Config::RecordCoordinates;
   using Fields = FieldOps::Set<typename HConfig::TrajectoryConfig::Fields>;

   using Background = Background<HConfig>;
   using NumericalDerivatives = NumericalDerivatives<Background, HConfig::numeric_derivatives()>;
   using Diffusion = Diffusion<HConfig>;

   using DiffusionConfig = HConfig::DiffusionConfig;
   using DiffusionCoordinates = typename DiffusionConfig::Coordinates;
   using DiffusionFields = typename FieldOps::Set<typename DiffusionConfig::Fields>;
   using DiffusionFieldsRemainder = typename FieldOps::Difference<DiffusionFields, Fields>;

   using ButcherTable = ButcherTable<HConfig::TrajectoryConfig::rk_integrator>;
   using Records = Records<RecordCoordinates, HConfig::specie, HConfig::TrajectoryConfig::record_mag_extrema, HConfig::TrajectoryConfig::record_trajectory>;

// Polymorphic base types
   using DistributionBase = DistributionBase<HConfig>;
   using BoundaryBase = BoundaryBase<HConfig>;
   using InitialBase = InitialBase<HConfig>;

protected:

//! Background object (persistent)
   Background background = Background();

//! Numerical Derivatives manager (persistent)
   NumericalDerivatives numericalderivatives = NumericalDerivatives();

//! Diffusion object (persistent)
   Diffusion diffusion = Diffusion();

//! Array of distribution objects (persistent)
   std::vector<std::shared_ptr<DistributionBase>> distributions;

//! Array of time boundary condition objects (persistent)
   std::vector<std::unique_ptr<BoundaryBase>> bcond_t;

//! Array of spatial boundary condition objects (persistent)
   std::vector<std::unique_ptr<BoundaryBase>> bcond_s;

//! Array of momentum boundary condition objects (persistent)
   std::vector<std::unique_ptr<BoundaryBase>> bcond_m;

//! Initial condition in time (persistent)
   std::unique_ptr<InitialBase> icond_t = nullptr;

//! Initial condition in space (persistent)
   std::unique_ptr<InitialBase> icond_s = nullptr;

//! Initial condition in momentum (persistent)
   std::unique_ptr<InitialBase> icond_m = nullptr;

   Coordinates _coords;

   ButcherTable butcher_table;

   Coordinates local_coords;

//! Background-dependent dmax (transient)
   double _dmax;

//! Spatial data (transient)
   Fields _fields;

   //! Slopes for position in RK step (transient)
   GeoVector slope_pos[ButcherTable::data.rk_stages];

//! Slopes for momentum in RK step (transient)
   GeoVector slope_mom[ButcherTable::data.rk_stages];

//! Actual time step (transient)
   double dt;

//! Physics based time step (transient)
   double dt_physical;

//! Time step from the adaptive scheme (transient)
   double dt_adaptive;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Number of reflections (transient)
   int n_refl;

//! Number of mirrorings (transient)
   int n_mirr;

//! Active time boundary (transient)
   int bactive_t;

//! Active spatial boundary (transient)
   int bactive_s;

//! Active momentum boundary (transient)
   int bactive_m;

//! Action for current boundary (transient)
   int action;

//! Distance to the nearest boundary (transient)
   double nearest_bnd_dist;

//! Coordinates at the start of the trajectory (needed for distributions)
   Coordinates coords0;

//! Field values at the start of the trajectory (needed for distributions)
   Fields fields0;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Default constructor (protected, class not designed to be instantiated)
   TrajectoryBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryBase(const std::string_view& name_in, status_t status_in);

//! Reset all boundary objects
   void ResetAllBoundaries(void);

//! Evaluate all boundary objects
   void ComputeAllBoundaries(void);

//! Update all boundary objects
   void UpdateAllBoundaries(void);

//! Load the last trajectory point
//   virtual void Load(void);

//! Add the current postion and momentum to the end of the trajectory (typically at the end of a time step)
   void Store(void);

//! Load the local trajectory point
   virtual void LoadLocal(void);

//! Save the current postion and momentum in local variables (during advance step)
   void StoreLocal(void);

//! Conversion of momentum from "native" to (p,mu,phi) coordinates
//   virtual GeoVector ConvertMomentum(void) const;

//! Predict whether any time boundaries will be crossed during the currect step
   void TimeBoundaryProximityCheck(void);

//! Fast test for any spatial boundary overshot
   bool SpaceTerminateCheck(void);

//! CommonFields for time and position, and output fields
   template <typename Coordinates, typename Fields, typename RequestedFields = Fields>
   void CommonFields(Coordinates&, Fields&);

//! Compute the RK slopes
   virtual void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) = 0;

//! Compute the physical time step
   virtual void PhysicalStep(void) = 0;

//! Computes RK slopes
   bool RKSlopes(void);

//! Take a step using precomputed RK slopes
   bool RKStep(void);

//! Handle boundaries at the end of time step
   void HandleBoundaries(void);

//! Advance trajectory using RK method
   bool RKAdvance(void);

//! Advance trajectory
   virtual bool Advance(void) = 0;

//! Function for adjusting/computing relevant momentum components in field aligned frame
   virtual void MomentumCorrection(void);

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   virtual bool IsSimulationReady(void) const;

public:

//! Recording and printing/exporting data helper class
   Records records;

//! Copy constructor (class not copyable)
   TrajectoryBase(const TrajectoryBase& other) = delete;

//! Destructor
   virtual ~TrajectoryBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<TrajectoryBase> Clone(void) const = 0;

//! Assign background model parameters
   void SetupBackground(const DataContainer& container_in);

   //! Assign diffusion model parameters
   void SetupDiffusion(const DataContainer& container_in);

//! Connect to an existing distribution object
   void ConnectDistribution(const std::shared_ptr<DistributionBase> distribution_in);

//! Disconnect an existing distribution object
   void DisconnectDistribution(int distro);

//! Replace an existing distribution object with another
   void ReplaceDistribution(int distro, const std::shared_ptr<DistributionBase> distribution_in);

//! Add a boundary condition
   void AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in);

//! Add an initial condition
   void AddInitial(const InitialBase& initial_in, const DataContainer& container_in);

//! Return the number of reflections (at boundaries)
   int Reflections(void) const {return n_refl;};

//! Return the number of mirror events in the trajectory
   int Mirrorings(void) const {return n_mirr;};

//! Return the time elapsed
   double ElapsedTime(void) const;

//! Signals the background that its services can start
   void StartBackground(void);

//! Signals the background that its services are no longer needed
   void StopBackground(void);

//! Clear the trajectory and start a new one
   virtual void SetStart(void);

//! Integrate the entrire trajectory
   void Integrate(void);

//! Return number of boundary crossing for a single boundary
   int Crossings(unsigned int output, unsigned int bnd) const;

////! Print various quantities along the trajectory
//   void PrintTrajectory(const std::string traj_name, bool phys_units, unsigned int output, unsigned int stride = 1, double dt_out = 0.0) const;
//
////! Print a trajectory as CSV
//   void PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride = 1) const;

//! Print the human-readable status
   void InterpretStatus(void) const;

//! Print the information about the object
   void PrintInfo(void) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryBase inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/08/2022
*/
template <typename HConfig>
inline void TrajectoryBase<HConfig>::ConnectDistribution(const std::shared_ptr<DistributionBase> distribution_in)
{
   distributions.push_back(distribution_in);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2022
\param[in] distro index of which distribution to replace
*/
template <typename HConfig>
inline void TrajectoryBase<HConfig>::ReplaceDistribution(int distro, const std::shared_ptr<DistributionBase> distribution_in)
{
   distributions[distro] = distribution_in;
};

/*!
\author Juan G Alonso Guzman
\date 07/15/2022
\param[in] distro index of which distribution to reset
*/
template <typename HConfig>
inline void TrajectoryBase<HConfig>::DisconnectDistribution(int distro)
{
   distributions[distro].reset();
};

///*!
//\author Vladimir Florinski
//\date 09/25/2020
//*/
//template <typename HConfig>
//inline void TrajectoryBase<HConfig>::Load(void)
//{
//   _coords.Vel() = Vel<Config::specie>(_coords.Mom());
//};

// todo deprecated
///*!
//\author Vladimir Florinski
//\date 09/25/2020
//*/
//template <typename HConfig>
//inline void TrajectoryBase<HConfig>::Store(void)
//{
//   records.Store(_coords);
//};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/10/2022
*/
template <typename HConfig>
inline void TrajectoryBase<HConfig>::LoadLocal(void)
{
   _coords = local_coords;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/10/2022
*/
template <typename HConfig>
inline void TrajectoryBase<HConfig>::StoreLocal(void)
{
   local_coords = _coords;
};

///*!
//\author Juan G Alonso Guzman
//\date 06/12/2024
//\return Momentum in (p,mu,phi) coordinates
//*/
//template <typename HConfig>
//inline GeoVector TrajectoryBase<HConfig>::ConvertMomentum(void) const
//{
//   return _coords.Mom();
//};

/*!
\author Vladimir Florinski
\date 01/13/2021
\return Total time spanned by the trajectory
*/
template <typename HConfig>
inline double TrajectoryBase<HConfig>::ElapsedTime(void) const
{
   return _coords.Time() - coords0.Time();
};

};

// Something like this is needed for templated classes
#include "trajectory_base.cc"

#endif
