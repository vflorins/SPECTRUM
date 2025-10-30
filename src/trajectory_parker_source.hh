/*!
\file trajectory_parker_source.hh
\brief Declares a class for trajectory based on the Parker Transport Equation
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef TRAJECTORY_PARKER_SOURCE_HH
#define TRAJECTORY_PARKER_SOURCE_HH

#include "trajectory_parker.hh"

namespace Spectrum {

//! Readable name of the TrajectoryParkerSource class
const std::string traj_name_parker_source = "TrajectoryParkerSource";

/*!
\brief A derived class for Parker equation (diffusive simulation) with a source term
\author Juan G Alonso Guzman

Components of "traj_mom" are: p_mag (x), unused (y), unused (z)
*/
class TrajectoryParkerSource : public TrajectoryParker {

//! Compute the RK slopes
   void Slopes(GeoVector& slope_pos_istage, GeoVector& slope_mom_istage, double& slope_amp_istage, double& slope_wgt_istage) override;

public:

//! Default constructor
   TrajectoryParkerSource(void);

//! Constructor with arguments (to speed up construction of derived classes)
   TrajectoryParkerSource(const std::string& name_in, unsigned int specie_in, uint16_t status_in, bool presize_in);

//! Copy constructor (class not copyable)
   TrajectoryParkerSource(const TrajectoryParkerSource& other) = delete;

//! Destructor
   ~TrajectoryParkerSource() override = default;

//! Clone function
   CloneFunctionTrajectory(TrajectoryParkerSource);

//! Clear the trajectory and start a new one with specified position and momentum
   void SetStart(void) override;

};

//! Trajectory type
#if TRAJ_TYPE == TRAJ_PARKER_SOURCE
typedef TrajectoryParkerSource TrajectoryType;
#endif

};

#endif
