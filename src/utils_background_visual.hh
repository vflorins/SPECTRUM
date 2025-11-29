/*!
\file background_base.hh
\brief Declares a base class to compute the plasma background
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_UTILS_BACKGROUND_VISUAL_HH
#define SPECTRUM_UTILS_BACKGROUND_VISUAL_HH


// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "common/params.hh"
#include "common/physics.hh"
#include "common/matrix.hh"
#include "common/derivativedata.hh"
#include "src/server_config.hh"

#include <memory>
#ifdef USE_SILO
#include <silo.h>
#else
#define DBfile nullptr_t
#define DB_CLOBBER 0
#define DB_LOCAL 0
#define DB_PDB 0
#define DB_DOUBLE 0
#define DB_COLLINEAR 0
#define DB_ZONECENT 0
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundVisual class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\brief A class to generate plasma background visualization
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

 */
template <typename HConfig_>
class BackgroundVisual {
public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

private:

   //! Name of the 2D mesh
   static constexpr std::string_view mesh2d_name = "square_mesh";

//! Name of the 3D mesh
   static constexpr std::string_view mesh3d_name = "cube_mesh";

protected:

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

//! Default constructor (protected, class not designed to be instantiated)
   BackgroundVisual(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundVisual(const std::string_view& name_in, status_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BackgroundVisual(const BackgroundVisual& other);


public:

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

};

};

// Something like this is needed for templated classes
#include "utils_background_visual.cc"

#endif
