/*!
\file trajectory_guiding_scatt.hh
\brief Declares a class for trajectory based on guiding center equations with pitch angle (elastic) scattering
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_SCATT_HH
#define SPECTRUM_TRAJECTORY_GUIDING_SCATT_HH

#include "trajectory_guiding.hh"

namespace Spectrum {

//! Whether to split the diffusive advance into two (one before and one after the advection).
// #define SPLIT_SCATT

//! Whether to use constant dmumax or constant dthetamax, 0 = constant dthetamax, 1 = constant dmumax
#define CONST_DMUMAX 0

//! Which stochastic method to use for PA scattering, 0 = Euler, 1 = Milstein, 2 = RK2
#define STOCHASTIC_METHOD_MU 0

//! Readable name of the TrajectoryGuidingScatt class
const std::string traj_name_guidingscatt = "TrajectoryGuidingScatt";

//! Default initial size
const unsigned int defsize_guidingscatt = 10000;

#ifdef SPLIT_SCATT
const double alpha = 0.5;
#endif

#if CONST_DMUMAX == 1
//! Desired accuracy in pitch angle cosine
const double dmumax = 0.02;
#else 
//! Desired accuracy in pitch angle (deg x [deg to rad conversion factor])
const double dthetamax = 2.0 * M_PI / 180.0;
#endif

//! CFL condition for pitch angle scattering
const double cfl_pa_gs = 0.1;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Trajectory tracer for the relativistic guiding center equations plus pitch angle scattering
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
class TrajectoryGuidingScatt : virtual public TrajectoryGuiding {

protected:

//! Pitch angle scattering coefficient (transient)
   double Dmumu;

//! Rate of change of Dmumu with mu (transient)
   double Vmu;

#ifdef GEO_DEBUG
//! Number of times |mu| > 1
   int Nabsmugt1 = 0;
#endif

//! Computes the pitch angle diffusion coefficient
   virtual void DiffusionCoeff(void);

//! Performs Euler pitch angle scattering (weak order 1, strong order 1/2)
   void EulerPitchAngleScatt(bool second);

//! Performs Milstein pitch angle scattering (weak order 1, strong order 1)
   void MilsteinPitchAngleScatt(bool second);

//! Performs RK2 pitch angle scattering (weak order 2, strong order 1?)
   void RK2PitchAngleScatt(bool second);

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   bool IsSimmulationReady(void) const override;

public:

//! Default constructor
   TrajectoryGuidingScatt(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuidingScatt(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Destructor
   ~TrajectoryGuidingScatt() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuidingScatt);

#ifdef GEO_DEBUG

//! Return number of |mu| = 1 crossings
   int fabsmugt1(void) const;

//! Reset number of |mu| = 1 crossings
   void ResetNabsmugt1(void);

#endif
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingScatt inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef GEO_DEBUG

/*!
\author Juan G Alonso Guzman
\date 05/03/2022
\return Number of |mu| = 1 crossings
*/
inline int TrajectoryGuidingScatt::fabsmugt1(void) const
{
   return Nabsmugt1;
};

/*!
\author Juan G Alonso Guzman
\date 05/03/2022
*/
inline void TrajectoryGuidingScatt::ResetNabsmugt1(void)
{
   Nabsmugt1 = 0;
};

#endif

//! Trajectory type
#if TRAJ_TYPE == TRAJ_GUIDING_SCATT
typedef TrajectoryGuidingScatt TrajectoryType;
#endif

};

#endif
