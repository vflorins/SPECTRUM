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
#include "distribution_base.hh"
#include "background_base.hh"
#include "diffusion_base.hh"
#include "boundary_base.hh"
#include "initial_base.hh"
#include <common/rk_config.hh>


namespace Spectrum {

//! Record trajectory flag
#define RECORD_TRAJECTORY

//! Record |B| extrema flag
// #define RECORD_BMAG_EXTREMA

//! Trajectory advance safety level: 0 means no checks, 1 means check dt only, 2 means check dt, number of segments, and time adaptations per step.
#define TRAJ_ADV_SAFETY_LEVEL 2

#if TRAJ_ADV_SAFETY_LEVEL == 2
//! Largest length for single trajectory
constexpr int max_trajectory_steps = 10000000;

//! Largest number of time step adaptations for a single time step
constexpr int max_time_adaptations = 100;
#endif

//! The trajecory will end after the step is completed
constexpr uint16_t TRAJ_FINISH = 0x0010;

//! Time boundary was crossed
constexpr uint16_t TRAJ_TIME_CROSSED = 0x0020;

//! Spatial boundary was crossed
constexpr uint16_t TRAJ_SPATIAL_CROSSED = 0x0040;

//! Momentum boundary was crossed
constexpr uint16_t TRAJ_MOMENTUM_CROSSED = 0x0080;

//! Trajectory is invalid and must be discarded
constexpr uint16_t TRAJ_DISCARD = 0x0100;

//! Clone function pattern
#define CloneFunctionTrajectory(T) std::unique_ptr<TrajectoryBase<Fields>> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if maximum number of steps is reached in trajectory
\author Juan G Alonso Guzman
*/
class ExMaxStepsReached : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2022
\return Text describing the error
*/
inline const char* ExMaxStepsReached::what(void) const noexcept
{
   return "Maximum number of steps in a trajectory reached.";
};

/*!
\brief Exception if maximum number of time step adaptations in a single step is reached
\author Juan G Alonso Guzman
*/
class ExMaxTimeAdaptsReached : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2022
\return Text describing the error
*/
inline const char* ExMaxTimeAdaptsReached::what(void) const noexcept
{
   return "Maximum number of time adaptations in a single step reached.";
};

/*!
\brief Exception if time step becomes too small
\author Juan G Alonso Guzman
*/
class ExTimeStepTooSmall : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2022
\return Text describing the error
*/
inline const char* ExTimeStepTooSmall::what(void) const noexcept
{
   return "Time step became too small.";
};

/*!
\brief Exception if time step becomes nan
\author Juan G Alonso Guzman
*/
class ExTimeStepNan : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2022
\return Text describing the error
*/
inline const char* ExTimeStepNan::what(void) const noexcept
{
   return "Time step became nan.";
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A base class to integrate a trajectory
\author Vladimir Florinski

A trajectory should be thought of as a self-contained simulation. There are four ingredients to a trajectory: (a) the transport physics that is built into the class itself, (b) the background u, E, B fields provided through "background", (c) the boundary conditions contained in "bcond_t", "bcond_s", and "bcond_m", and (d) the initial conditions provided by "icond_s" and "icond_m". The trajectory is responsible for computing the derived transport coefficients. This is done for efficiency purposes because a separate transport class hierarchy would have to interact with the other components and exchanging different kinds of transport parameters must be done through a container, hence require extra load/store operations.

A trajectory object is considered initialized if (a) background is assigned, (b) at least one time boundary is assigned, (c) the space initial condition is assigned, and (d) the momentum initial condition is assigned.
*/
template <typename Fields_>
class TrajectoryBase : public Params {

public:

   using Fields = Fields_;
   using DistributionBase = DistributionBase<Fields>;
   using BackgroundBase = BackgroundBase<Fields>;
   using DiffusionBase = DiffusionBase<Fields>;

protected:

//! Initial length of trajectory containers (persistent)
   unsigned int presize = 1;

//! Particle's charge to mass ratio (persistent)
   double q;

//! Array of distribution objects (persistent)
   std::vector<std::shared_ptr<DistributionBase>> distributions;

//! Background object (persistent)
   std::unique_ptr<BackgroundBase> background = nullptr;

//! Diffusion object (persistent)
   std::unique_ptr<DiffusionBase> diffusion = nullptr;

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

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Time along trajectory, includes segment counter (transient)
   std::vector<double> traj_t;

//! Position along trajectory (transient)
   std::vector<GeoVector> traj_pos;

//! Momentum along trajectory (transient)
   std::vector<GeoVector> traj_mom;

#ifndef RECORD_TRAJECTORY
//! Number of trajectory segments (transient)
   int n_segs;
#endif

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

//! Local time for Advance function (transient)
   double local_t;

//! Local position for Advance function (transient)
   GeoVector local_pos;

//! Local momentum for Advance function (transient)
   GeoVector local_mom;

//! Spatial data (transient)
   Fields _fields;

//! Spatial data at the start of the trajectory (transient)
   Fields fields0;

//! Slopes for position in RK step (transient)
   GeoVector slope_pos[MAX_RK_STAGES];

//! Slopes for momentum in RK step (transient)
   GeoVector slope_mom[MAX_RK_STAGES];

//! Actual time step (transient)
   double dt;

//! Physics based time step (transient)
   double dt_physical;

//! Time step from the adaptive scheme (transient)
   double dt_adaptive;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Default constructor (protected, class not designed to be instantiated)
   TrajectoryBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Clear the content and set the initial container capacity to improve performance
   void PreSize(int init_cap);

//! Find the nearest time point
   void GetIdx(double t_in, int& pt, double& weight) const;

//! Reset all boundary objects
   void ResetAllBoundaries(void);

//! Evaluate all boundary objects
   void ComputeAllBoundaries(void);

//! Update all boundary objects
   void UpdateAllBoundaries(void);

//! Load the last trajectory point
   virtual void Load(void);

//! Add the current postion and momentum to the end of the trajectory (typically at the end of a time step)
   void Store(void);

//! Load the local trajectory point
   virtual void LoadLocal(void);

//! Save the current postion and momentum in local variables (during advance step)
   void StoreLocal(void);

//! Conversion of momentum from "native" to (p,mu,phi) coordinates
   virtual GeoVector ConvertMomentum(void) const;

//! Momentum transformation on reflection at a boundary
   virtual void ReverseMomentum(void);

//! Predict whether any time boundaries will be crossed during the currect step
   void TimeBoundaryProximityCheck(void);

//! Fast test for any spatial boundary overshot
   bool SpaceTerminateCheck(void);

//! Compute the common fields ("Uvec", etc.) and their dependencies ("Bmag", etc.)
   void CommonFields(void);

//! Overloaded CommonFields for custom time and position, and output field
   void CommonFields(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, SpatialData& spdata);

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

#ifdef RECORD_BMAG_EXTREMA
//! Update |B| maximum and minimum values along trajectory
   void UpdateBmagExtrema(void);
#endif

//! Advance trajectory
   virtual bool Advance(void) = 0;

//! Function for adjusting/computing relevant momentum components in field aligned frame
   virtual void MomentumCorrection(void);

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   virtual bool IsSimmulationReady(void) const;

public:

//! Copy constructor (class not copyable)
   TrajectoryBase(const TrajectoryBase& other) = delete;

//! Destructor
   virtual ~TrajectoryBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<TrajectoryBase> Clone(void) const = 0;

//! Set the particle specie
   void SetSpecie(unsigned int specie_in);

//! Connect to an existing distribution object 
   void ConnectDistribution(const std::shared_ptr<DistributionBase> distribution_in);

//! Disconnect an existing distribution object
   void DisconnectDistribution(int distro);

//! Replace an existing distribution object with another
   void ReplaceDistribution(int distro, const std::shared_ptr<DistributionBase> distribution_in);

//! Add a background object
   void AddBackground(const BackgroundBase& background_in, const DataContainer& container_in);

//! Assign diffusion model parameters
   void AddDiffusion(const DiffusionBase& diffusion_in, const DataContainer& container_in);

//! Add a boundary condition
   void AddBoundary(const BoundaryBase& boundary_in, const DataContainer& container_in);

//! Add an initial condition
   void AddInitial(const InitialBase& initial_in, const DataContainer& container_in);

#ifdef RECORD_BMAG_EXTREMA
//! Get |B| minimum
   double GetBmagMin(void) const;

//! Get |B| maximum
   double GetBmagMax(void) const;
#endif

//! Return the position at a given time
   GeoVector GetPosition(double t_in) const;

//! Return the velocity at a given time
   GeoVector GetVelocity(double t_in) const;

//! Return the kinetic energy at a given time
   double GetEnergy(double t_in) const;

//! Return the distance along the trajectory at a given time
   double GetDistance(double t_in) const;

//! Return the number of segments in the trajectory
   int Segments(void) const;

//! Return the number of reflections (at boundaries)
   int Reflections(void) const {return n_refl;};

//! Return the number of mirror events in the trajectory
   int Mirrorings(void) const {return n_mirr;};

//! Return the time elapsed
   double ElapsedTime(void) const;

//! Signals the background that its services are no longer needed
   void StopBackground(void);

//! Clear the trajectory and start a new one
   virtual void SetStart(void);

//! Integrate the entrire trajectory
   void Integrate(void);

//! Return number of boundary crossing for a single boundary
   int Crossings(unsigned int output, unsigned int bnd) const;

//! Print various quantities along the trajectory
   void PrintTrajectory(const std::string traj_name, bool phys_units, unsigned int output, unsigned int stride = 1, double dt_out = 0.0) const;

//! Print a trajectory as CSV
   void PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride = 1) const;

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
template <typename Fields>
inline void TrajectoryBase<Fields>::ConnectDistribution(const std::shared_ptr<DistributionBase> distribution_in)
{
   distributions.push_back(distribution_in);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2022
\param[in] distro index of which distribution to replace
*/
template <typename Fields>
inline void TrajectoryBase<Fields>::ReplaceDistribution(int distro, const std::shared_ptr<DistributionBase> distribution_in)
{
   distributions[distro] = distribution_in;
};

/*!
\author Juan G Alonso Guzman
\date 07/15/2022
\param[in] distro index of which distribution to reset
*/
template <typename Fields>
inline void TrajectoryBase<Fields>::DisconnectDistribution(int distro)
{
   distributions[distro].reset();
};

/*!
\author Vladimir Florinski
\date 09/25/2020
*/
template <typename Fields>
inline void TrajectoryBase<Fields>::Load(void)
{
#ifdef RECORD_TRAJECTORY
   _t = traj_t.back();
   _pos = traj_pos.back();
   _mom = traj_mom.back();
#endif
   _vel = Vel(_mom, specie);
};

/*!
\author Vladimir Florinski
\date 09/25/2020
*/
template <typename Fields>
inline void TrajectoryBase<Fields>::Store(void)
{
#ifdef RECORD_TRAJECTORY
   traj_t.push_back(_t);
   traj_pos.push_back(_pos);
   traj_mom.push_back(_mom);
#else
   n_segs++;
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/10/2022
*/
template <typename Fields>
inline void TrajectoryBase<Fields>::LoadLocal(void)
{
   _t = local_t;
   _pos = local_pos;
   _mom = local_mom;
   _vel = Vel(_mom, specie);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/10/2022
*/
template <typename Fields>
inline void TrajectoryBase<Fields>::StoreLocal(void)
{
   local_t = _t;
   local_pos = _pos;
   local_mom = _mom;
};

/*!
\author Juan G Alonso Guzman
\date 06/12/2024
\return Momentum in (p,mu,phi) coordinates
*/
template <typename Fields>
inline GeoVector TrajectoryBase<Fields>::ConvertMomentum(void) const
{
   return _mom;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\return Largest index in the trajectory arrays
*/
template <typename Fields>
inline int TrajectoryBase<Fields>::Segments(void) const
{
#ifdef RECORD_TRAJECTORY
   return traj_t.size() - 1;
#else
   return n_segs;
#endif
};

/*!
\author Vladimir Florinski
\date 01/13/2021
\return Total time spanned by the trajectory
*/
template <typename Fields>
inline double TrajectoryBase<Fields>::ElapsedTime(void) const
{
#ifdef RECORD_TRAJECTORY
   return traj_t.back() - traj_t.front();
#else
   return _t - traj_t[0];
#endif
};

};

#endif
