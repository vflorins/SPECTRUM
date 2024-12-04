/*!
\file background_base.hh
\brief Declares a base class to compute the plasma background
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_BASE_HH
#define SPECTRUM_BACKGROUND_BASE_HH

#include "config.h"

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "common/params.hh"
#include "common/physics.hh"
#include "common/spatial_data.hh"
#include "src/server_config.hh"

#include <memory>
#ifdef USE_SILO
#include <silo.h>
#endif

namespace Spectrum {

//! Number of local coordinate systems in which to evaluate and average numerical derivatives
#define BACKGROUND_NUM_GRAD_EVALS 1

#ifdef USE_SILO

//! Name of the 2D mesh
const std::string mesh2d_name = "square_mesh";

//! Name of the 3D mesh
const std::string mesh3d_name = "cube_mesh";

#endif

//! What fraction of "_dmax" to use to calculate the field increment
const double incr_dmax_ratio = 0.0001;

//! Rotation angle of local (x,y) plane in numerical derivative evaluation/averaging
const double local_rot_ang = M_PI / BACKGROUND_NUM_GRAD_EVALS;

//! Sine and cosine of local rotation angle
const double sin_lra = sin(local_rot_ang);
const double cos_lra = cos(local_rot_ang);

//! Clone function pattern
#define CloneFunctionBackground(T) std::unique_ptr<BackgroundBase> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if field evaluation failed
\author Vladimir Florinski
*/
class ExFieldError : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Vladimir Florinski
\date 12/09/2021
\return Text describing the error
*/
inline const char* ExFieldError::what(void) const noexcept
{
   return "Field evaluation error";
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A base class to generate plasma background
\author Vladimir Florinski

The BackgroundXXXXX" classes describe plasma condition (u, E, B, dB/dt, gradB) at an arbitrary location. In the simplest scenario these are computed analytically from a known formula. A more complex situation is when the background is provided externally, typically from a file produced from a numerical simulation. Derived classes may contain an interface to read mesh based output files and calculate field values using interpolation.

Parameters: double t0, GeoVector r0, GeoVector u0, GeoVector B0, double dmax0
*/
class BackgroundBase : public Params {

protected:

//! Initial time (persistent)
   double t0;

//! Coordinate system origin (persistent)
   GeoVector r0;

//! Reference flow vector (persistent)
   GeoVector u0;

//! Reference magnetic field (persistent)
   GeoVector B0;

//! Distance reference value, how far the trajectory is allowed to move in one step (persistent)
   double dmax0;

//! Spatial data object (transient)
   SpatialData _spdata;

//! Auxiliary spatial data objects for derivative computations (transient)
   SpatialData _spdata_tmp, _spdata_forw, _spdata_back;

//! Gyro-radius for derivative computations (transient)
   double r_g;

//! Gyro-frequency for derivative computations (transient)
   double w_g;

//! Field-aligned basis (transient)
   GeoMatrix fa_basis;

//! Rotation matrix (transient)
   GeoMatrix rot_mat;

#ifdef USE_SILO

//! A handle to a SILO database (transient)
   DBfile* silofile;

//! First corner of the output box (transient)
   GeoVector xyz_min;

//! Second corner of the output box (transient)
   GeoVector xyz_max;

//! Number of output zones in each direction (transient)
   MultiIndex dims_z;

//! Direction of normal for plane cut (transient)
   GeoVector plane;

//! Direction of the x axis in the cut plane (transient)
   GeoVector x_dir;

#endif

//! Default constructor (protected, class not designed to be instantiated)
   BackgroundBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BackgroundBase(const BackgroundBase& other);

//! Compute the fields at an incremented position or time
   void DirectionalDerivative(int xyz);

//! Compute the field derivatives based on the mask
   void NumericalDerivatives(void);

//! Set up the field evaluator based on "params"
   virtual void SetupBackground(bool construct);

//! Compute the internal u, B, and E fields
   virtual void EvaluateBackground(void);

//! Compute the internal derivatives of the fields
   virtual void EvaluateBackgroundDerivatives(void);

//! Calculate the maximum distance allowed per time step
   virtual void EvaluateDmax(void);

//! Calculate magnetic field magnitude
   virtual void EvaluateBmag(void);

public:

//! Destructor
   virtual ~BackgroundBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<BackgroundBase> Clone(void) const = 0;

//! Set up the object's persistent class data members
   void SetupObject(const DataContainer& cont_in);

//! Signal the backend this client no longer needs its service
   virtual void StopServerFront(void);

//! Find "safe" increment in a given direction
   double GetSafeIncr(const GeoVector& dir);

//! Return maximum distance allowed per time step
   double GetDmax(void) const;

//! Return fields at the internal position, evaluated or previously stored
   void GetFields(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, SpatialData& spdata);

#ifdef USE_SILO

//! Set up the plot limits
   void SetBox(const GeoVector& xyz_min_in, const GeoVector& xyz_max_in, const MultiIndex& dims_z_in,
               const GeoVector& normal, const GeoVector& right);

//! Generate a 3D mesh
   void BoxPlot3DMesh(const std::string box_fname, bool phys_units);

//! Generate a 3D box scalar plot
   void BoxPlot3DScalar(const std::string var_name, bool phys_units, double t = 0.0);

//! Generate a 3D box vector plot
   void BoxPlot3DVector(const std::string var_name, bool phys_units, double t = 0.0);

//! Generate a 2D mesh
   void BoxPlot2DMesh(const std::string box_fname, bool phys_units);

//! Generate a 2D box scalar plot
   void BoxPlot2DScalar(const std::string var_name, bool phys_units, double t = 0.0);

//! Generate a 2D box vector plot
   void BoxPlot2DVector(const std::string var_name, bool phys_units, double t = 0.0);

//! Finalize the output
   void BoxPlotFinalize(void);

#endif

};

};

#endif
