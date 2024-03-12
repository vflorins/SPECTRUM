/*!
\file trajectory_guiding_diff.hh
\brief Declares a class for trajectory based on guiding center equations with perpendicular diffusion
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_GUIDING_DIFF_HH
#define SPECTRUM_TRAJECTORY_GUIDING_DIFF_HH

#include "trajectory_guiding.hh"
#include "diffusion_base.hh"

namespace Spectrum {

//! Which stochastic method to use for perpendicular diffusion, 0 = Euler, 1 = Milstein, 2 = RK2
#define STOCHASTIC_METHOD_PERP 0

//! Readable name of the TrajectoryGuidingDiff class
const std::string traj_name_guidingdiff = "TrajectoryGuidingDiff";

//! Default initial size
const unsigned int defsize_guidingdiff = 10000;

//! CFL condition for perpendicular diffusion
const double cfl_dif_gd = 0.5;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TrajectoryGuidingDiff class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Trajectory tracer for the relativistic guiding center equations plus perpendicular diffusion
\author Juan G Alonso Guzman
\author Vladimir Florinski
*/
class TrajectoryGuidingDiff : virtual public TrajectoryGuiding {

protected:

//! The field aligned coordinate system unit vectors (transient)
   GeoVector fa_basis[3];

//! Perpendicular diffusion coefficient (transient)
   double Dperp;

//! Gradient of Dperp (transient)
   GeoVector Vperp;

//! Stochastic increment in guiding center
   GeoVector dr_perp;

//! Compute the basis vectors of the field aligned frame
   void FieldAlignedFrame(void);

//! Computes the perpendicular diffusion coefficient
   virtual void DiffusionCoeff(void);

//! Computes perpendicular diffusion slopes using Euler method
   void EulerPerpDiffSlopes(void);

//! Computes perpendicular diffusion slopes using Milstein method
   void MilsteinPerpDiffSlopes(void);

//! Computes perpendicular diffusion slopes using RK2 method
   bool RK2PerpDiffSlopes(void);

//! Compute the RK(?) slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage) override;

//! Compute the physical time step
   void PhysicalStep(void) override;

//! Take a step
   bool Advance(void) override;

//! Perform all checks to see if a trajectory is ready to be used in a simulation
   bool IsSimmulationReady(void) const override;

public:

//! Default constructor
   TrajectoryGuidingDiff(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryGuidingDiff(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Destructor
   ~TrajectoryGuidingDiff() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryGuidingDiff);

};

//! Trajectory type
#if TRAJ_TYPE == TRAJ_GUIDING_DIFF
typedef TrajectoryGuidingDiff TrajectoryType;
#endif

};

#endif
