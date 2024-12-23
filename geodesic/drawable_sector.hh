/*!
\file drawable_sector.hh
\brief Declares PostScript graphic routines to draw a sector and all its elements
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DRAWABLE_SECTOR_HH
#define SPECTRUM_DRAWABLE_SECTOR_HH

#include <string>
#include <fstream>
#include "geodesic/geodesic_sector.hh"

namespace Spectrum {

//! Default PostScript resolution
constexpr double ps_dpi_def = 72.0;

//! Resolution of PostScript drawings
constexpr double ps_dpi = 600.0;

//! Scale factor of PostScript drawings
constexpr double ps_scale_factor = ps_dpi_def / ps_dpi;

//! File name prefix for the face plot
const std::string face_file_pref = "face";

//! File name prefix for the edge plot
const std::string edge_file_pref = "edge";

//! File name prefix for the vertex plot
const std::string vert_file_pref = "vert";

/*!
\brief A class to graphically visualize the element arrangement within a secrtor
\author Vladimir Florinski
*/
template <int verts_per_face>
class DrawableSector : public GeodesicSector<verts_per_face>
{
protected:

   using PolygonalAddressing<verts_per_face>::cardinal_directions;
   using PolygonalAddressing<verts_per_face>::square_fill;
   using PolygonalAddressing<verts_per_face>::edge_dimax;
   using PolygonalAddressing<verts_per_face>::edge_djmax;
   using GeodesicSector<verts_per_face>::side_length;
   using GeodesicSector<verts_per_face>::total_length;
   using GeodesicSector<verts_per_face>::ghost_width;
   using GeodesicSector<verts_per_face>::MaxVertJ;
   using GeodesicSector<verts_per_face>::MaxFaceJ;

//! Plot border width
   static constexpr double border = 30.0 / ps_scale_factor;

//! Width of image
   static constexpr double v_plot = 8.5 * ps_dpi - 2 * border;

//! Left side of the bounding box
   static constexpr int bbox_l = 15;

//! Right side of the bounding box
   static constexpr int bbox_r = 597;

//! Scaling factor for the line thickness and font sizes
   double annotation_scale_factor;

//! Character size as a fraction of the length of the edge
   double char_scale;

//! Thick line width as a fraction of the length of the edge
   double thick_line_frac;

//! Vertical shift of thick lines as a fraction of the length of the edge
   double line_shift;

//! Bottom side of the bounding box
   int bbox_b;

//! Top side of the bounding box
   int bbox_t;

//! Vertical rise of the image to center
   double rise;

//! Array of lengths for use by the drawing routines
   double len[7];

//! PostScript file where all output is sent
   std::ofstream mapfile;

//! Obtain a vertex drawing position from its TAS/QAS indices
   void GetCoord(int i, int j, double& x, double& y);

//! Draw a closed convex polygon with thick lines
   void DrawConvexPolygon(int sides, const double* x, const double* y);

//! Return a horizontal shift for string of numbers depending on the number of digits
   double GetShift(int num);

//! Print the PostScript header
   void MakeHeader(void);

//! Draw all faces including shading for TAS
   void DrawFaces(void);

//! Draw the borders of the neighbor exchange regions
   void DrawRegions(void);

//! Print a number for each vertex
   void PrintVertContent(double* data);

//! Print a number for each edge
   void PrintEdgeContent(double* data);

//! Print a number for each face
   void PrintFaceContent(double* data);

public:

//! Default constructor
   DrawableSector(void) = default;

//! Constructor with arguments
   DrawableSector(int width, int wghost);

//! Allocate memory and set up connectivity between mesh elements
   void SetDimensions(int width, int wghost, bool construct);

//! Generate a PostScript drawing of a sector with data
   void Draw(int request, double* data = nullptr);
};

};

#endif
