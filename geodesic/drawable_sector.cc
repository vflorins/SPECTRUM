/*!
\file drawable_sector.cc
\brief Defines PostScript graphic routines to draw a sector and all its elements
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/vectors.hh"
#include "drawable_sector.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DrawableSector public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/10/2024
*/
template <int verts_per_face>
DrawableSector<verts_per_face>::DrawableSector(int width, int wghost)
                              : GeodesicSector<verts_per_face>::GeodesicSector(width, wghost)
{
   SetDimensions(width, wghost, true);
};

/*!
\author Vladimir Florinski
\date 06/10/2024
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::SetDimensions(int width, int wghost, bool construct)
{
// Call base method
   if (!construct) GeodesicSector<verts_per_face>::SetDimensions(width, wghost, false);

// We would like to have the _height_ of the face to be identical for triangular and quad sectors, so the unit length "len[0]" will be different. To ensure that the font size and line thickness are the same, "annotation_scale_factor" is set to compensate for the difference in "len[0]".
   if (verts_per_face == 3) {
      bbox_b = 140;
      bbox_t = 652;
      rise = 128.0 / ps_scale_factor;
      annotation_scale_factor = 1.0;
      len[0] = v_plot / total_length / annotation_scale_factor;
      len[1] = len[0] * sin(M_PI / 3.0);
      len[2] = len[0] * cos(M_PI / 3.0);
      len[3] = len[0] / sin(M_PI / 3.0);
      len[4] = len[0] / cos(M_PI / 3.0);
      len[5] = len[0] / tan(M_PI / 3.0);
      len[6] = 0.0;
   }
   else if (verts_per_face == 4) {
      bbox_b = 70;
      bbox_t = 652;
      rise = 58.0 / ps_scale_factor;
      annotation_scale_factor = 2.0 / sqrt(3.0) * (side_length + 3 * ghost_width) / (side_length + 2 * ghost_width);
      len[0] = v_plot / total_length / annotation_scale_factor;
      len[1] = len[0];
      len[2] = 0.0;
      len[3] = 0.0;
      len[4] = 0.0;
      len[5] = 0.0;
      len[6] = len[0];
   };

   char_scale = 0.25 * annotation_scale_factor;
   thick_line_frac = 0.05 * annotation_scale_factor;
   line_shift = thick_line_frac / 2.0;
};

/*!
\author Vladimir Florinski
\date 06/10/2024
\param[in] request Request type: 1 for vertices, 2 for edges, 3 for faces
\param[in] data    Optional data associated with this sector
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::Draw(int request, double* data)
{
   std::string filename;

   if (request == 1) filename = vert_file_pref;
   else if (request == 2) filename = edge_file_pref;
   else if (request == 3) filename = face_file_pref;
   else return;

// Open the PS file
   filename += "_" + std::to_string(verts_per_face) + "_" + std::to_string(side_length) + "_" + std::to_string(ghost_width) + ".eps";
   mapfile.open(filename.c_str());
   
// Perform draw operations
   MakeHeader();
   DrawFaces();
   DrawRegions();

// Print the data
   if (request == 1) PrintVertContent(data);
   else if (request == 2) PrintEdgeContent(data);
   else if (request == 3) PrintFaceContent(data);

// Close the file
   mapfile << "showpage\n";
   mapfile << "%%EOF\n";
   mapfile.close();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DrawableSector protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/10/2024
\param[in]  i First index
\param[in]  j Second index
\param[out] x Abscissa coordinate
\param[out] y Ordinate coordinate
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::GetCoord(int i, int j, double& x, double& y)
{
   x = len[0] * i - len[2] * j;
   y = len[1] * j;
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] sides   Number of sides
\param[in] x       x coordinates of points
\param[in] y       y coordinates of points
\param[in] len     Array of lengths for use by the drawing routines
\param[in] mapfile PostScript file
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::DrawConvexPolygon(int sides, const double* x, const double* y)
{
   int vert;
   GeoVector v1, v2, v3;

// Interior points
   double* xint = new double[sides];
   double* yint = new double[sides];

// Compute the interior points. This assumes that the polygon is traversed in the CC direction.
   for (vert = 0; vert < sides; vert++) {
      v1[0] = x[(vert + sides - 1) % sides] - x[vert];
      v1[1] = y[(vert + sides - 1) % sides] - y[vert];
      v1[2] = 0.0;
      v2[0] = x[(vert + 1) % sides] - x[vert];
      v2[1] = y[(vert + 1) % sides] - y[vert];
      v2[2] = 0.0;
      v1.Normalize();
      v2.Normalize();
      v3 = Bisect(v1, v2);
      v3 *= thick_line_frac * len[0] / sqrt(1.0 - Sqr(v2 * v3));
      xint[vert] = v3[0] + x[vert];
      yint[vert] = v3[1] + y[vert];
   };

// Draw the border as a collection of trapezoids
   for (vert = 0; vert < sides; vert++) {
      mapfile << "newpath\n";
      mapfile << x[vert] + border << " " << y[vert] + border + rise << " moveto\n";
      mapfile << x[(vert + 1) % sides] + border << " " << y[(vert + 1) % sides] + border + rise << " lineto\n";
      mapfile << xint[(vert + 1) % sides] + border << " " << yint[(vert + 1) % sides] + border + rise << " lineto\n";
      mapfile << xint[vert] + border << " " << yint[vert] + border + rise << " lineto\n";
      mapfile << "closepath\n";
      mapfile << "fill\n";
   };

   delete[] xint;
   delete[] yint;
};

/*!
\author Vladimir Florinski
\date 06/10/2024
\param[in] num An integer 1-3 digits in length
\return Horizontal shift
*/
template <int verts_per_face>
double DrawableSector<verts_per_face>::GetShift(int num)
{
   if (num > 99) return -0.8 * char_scale * len[0];
   else if (num > 9) return -0.5 * char_scale * len[0];
   else return -0.2 * char_scale * len[0];
};

/*!
\author Vladimir Florinski
\date 06/10/2024
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::MakeHeader(void)
{
// Title and size
   mapfile << "%!PS-Adobe-3.0\n";
   mapfile << "%%BoundingBox: " << bbox_l << " " << bbox_b << " " << bbox_r << " " << bbox_t << std::endl;
   mapfile << "%%Creator: DrawSector\n";
   mapfile << "%%DocumentData: Clean8Bit\n";
   mapfile << "%%Title: Sector mesh and data\n";
   mapfile << "%%Orientation: Portrait\n";
   mapfile << "%%Pages: 1\n";
   mapfile << "%%EndComments\n";
   mapfile << std::endl;

// Scale factor
   mapfile << ps_scale_factor << " " << ps_scale_factor << " scale\n";

// Color palette
   mapfile << "/Color0 {1.0000 1.0000 1.0000} def\n";
   mapfile << "/Color1 {0.0000 0.0000 0.0000} def\n";
   mapfile << "/Color2 {1.0000 0.0000 0.0000} def\n";
   mapfile << "/Color3 {0.0000 1.0000 0.0000} def\n";
   mapfile << "/Color4 {0.0000 0.0000 1.0000} def\n";
   mapfile << "/Color5 {1.0000 1.0000 0.0000} def\n";
   mapfile << "/Color6 {0.7373 0.5608 0.5608} def\n";
   mapfile << "/Color7 {0.8627 0.8627 0.8627} def\n";
   mapfile << "/Color8 {0.5804 0.0000 0.8275} def\n";
   mapfile << "/Color9 {0.0000 1.0000 1.0000} def\n";
   mapfile << "/Color10 {1.0000 0.0000 1.0000} def\n";
   mapfile << "/Color11 {1.0000 0.6471 0.0000} def\n";
   mapfile << "/Color12 {0.4471 0.1294 0.7373} def\n";
   mapfile << "/Color13 {0.4039 0.0275 0.2824} def\n";
   mapfile << "/Color14 {0.2510 0.8784 0.8157} def\n";
   mapfile << "/Color15 {0.0000 0.5451 0.0000} def\n";
   mapfile << "/Color16 {0.7529 0.7529 0.7529} def\n";
   mapfile << "/Color17 {0.5059 0.5059 0.5059} def\n";
   mapfile << "/Color18 {0.2588 0.2588 0.2588} def\n";
   mapfile << "/DeviceRGB setcolorspace\n";

// Font name and size
   mapfile << "/Helvetica findfont " << char_scale * len[0] << " scalefont setfont\n";
};

/*!
\author Vladimir Florinski
\date 06/10/2024
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::DrawFaces(void)
{
   double x0, y0;

   if (verts_per_face == 3) {

// Yellow fill color. Shaded faces must be drawn first and the borders of unshaded faces are drawn on top.
      mapfile << "Color5 setcolor\n";

// Loop over shaded faces (TAS). Drawing starts from the bottom corner.
      for (auto i = 0; i < total_length; i++) {
         for (auto j = 1; j < MaxFaceJ(total_length, i); j += square_fill) {
            GetCoord(i, j / square_fill, x0, y0);
            mapfile << "newpath\n";
            mapfile << border + x0 << " " << border + rise + y0 << " moveto\n";
            mapfile << len[2] << " " << len[1] << " rlineto\n";
            mapfile << -len[0] << " " << 0 << " rlineto\n";
            mapfile << "closepath\n";
            mapfile << "fill\n";
         };
      };
   };

// Black border color
   mapfile << "Color1 setcolor\n";

// Loop over all faces (QAS) or unshaded faces (TAS). Drawing starts from the bottom left corner.
   for (auto i = 0; i < total_length; i++) {
      for (auto j = 0; j <= MaxFaceJ(total_length, i); j += square_fill) {
         GetCoord(i, j / square_fill, x0, y0);
         mapfile << "newpath\n";
         mapfile << border + x0 << " " << border + rise + y0 << " moveto\n";
         mapfile << len[0] << " " << 0 << " rlineto\n";
         mapfile << -len[2] << " " << len[1] << " rlineto\n";

         if (verts_per_face == 4) {
            mapfile << -len[0] << " " << 0 << " rlineto\n";
         };

         mapfile << "closepath\n";
         mapfile << "stroke\n";
      };
   };
};

/*!
\brief Draw the borders of the neighbor exchange regions
\author Vladimir Florinski
\date 06/10/2024
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::DrawRegions(void)
{
   double x[4], y[4];
   std::pair<int, int> base_vert = std::make_pair(square_fill * ghost_width, ghost_width);

// These shifts are always zero for QUAD mesh
   int tri_shift_side = (square_fill - 1) * (total_length - 3 * ghost_width);
   int tri_shift_ghst = (square_fill - 1) * ghost_width;

   mapfile << thick_line_frac * len[0] << " setlinewidth\n";

// Interior
   mapfile << "Color1 setcolor\n";

   GetCoord(square_fill * ghost_width, ghost_width, x[0], y[0]);
   GetCoord(total_length - ghost_width, ghost_width, x[1], y[1]);
   GetCoord(total_length - ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[2], y[2]);
   GetCoord(square_fill * ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[3], y[3]);
   DrawConvexPolygon(verts_per_face, x, y);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Edge neighbors
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Bottom or S edge neighbor
   mapfile << "Color10 setcolor\n";

   GetCoord(total_length - ghost_width, ghost_width, x[0], y[0]);
   GetCoord(square_fill * ghost_width, ghost_width, x[1], y[1]);
   GetCoord(square_fill * ghost_width, 0, x[2], y[2]);
   GetCoord(total_length - ghost_width - tri_shift_ghst, 0, x[3], y[3]);
   DrawConvexPolygon(4, x, y);

// Right or NE edge neighbor
   mapfile << "Color4 setcolor\n";

   GetCoord(total_length - ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[0], y[0]);
   GetCoord(total_length - ghost_width, ghost_width, x[1], y[1]);
   GetCoord(total_length, ghost_width + tri_shift_ghst, x[2], y[2]);
   GetCoord(total_length, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[3], y[3]);
   DrawConvexPolygon(4, x, y);

// Top edge neighbor
   if (verts_per_face == 4) {
      mapfile << "Color9 setcolor\n";

      GetCoord(square_fill * ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, square_fill * ghost_width), x[0], y[0]);
      GetCoord(total_length - ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[1], y[1]);
      GetCoord(total_length - ghost_width, MaxVertJ(total_length, total_length - ghost_width), x[2], y[2]);
      GetCoord(square_fill * ghost_width, MaxVertJ(total_length, square_fill * ghost_width), x[3], y[3]);
      DrawConvexPolygon(4, x, y);
   };

// Left or NW edge neighbor
   mapfile << "Color2 setcolor\n";

   GetCoord(square_fill * ghost_width, ghost_width, x[0], y[0]);
   GetCoord(square_fill * ghost_width + tri_shift_side, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[1], y[1]);
   GetCoord(tri_shift_side + tri_shift_ghst, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[2], y[2]);
   GetCoord(2 * tri_shift_ghst, ghost_width + tri_shift_ghst, x[3], y[3]);
   DrawConvexPolygon(4, x, y);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Vertex neighbors
//----------------------------------------------------------------------------------------------------------------------------------------------------

// SW vertex neighbors (2) - TRI mesh only
   if (verts_per_face == 3) {
      mapfile << "Color3 setcolor\n";

      GetCoord(square_fill * ghost_width, ghost_width, x[0], y[0]);
      GetCoord(square_fill * ghost_width, ghost_width + tri_shift_ghst, x[1], y[1]);
      GetCoord(square_fill * ghost_width - tri_shift_ghst, ghost_width, x[2], y[2]);
      DrawConvexPolygon(verts_per_face, x, y);

      GetCoord(square_fill * ghost_width, ghost_width, x[0], y[0]);
      GetCoord(square_fill * ghost_width - tri_shift_ghst, 0, x[1], y[1]);
      GetCoord(square_fill * ghost_width, 0, x[2], y[2]);
      DrawConvexPolygon(verts_per_face, x, y);
   };

// SW vertex neighbor
   mapfile << "Color15 setcolor\n";

   GetCoord(square_fill * ghost_width, ghost_width, x[0], y[0]);
   GetCoord(tri_shift_ghst, ghost_width, x[1], y[1]);
   GetCoord(tri_shift_ghst, 0, x[2], y[2]);
   GetCoord(square_fill * ghost_width, 0, x[3], y[3]);
   DrawConvexPolygon(verts_per_face, x, y);

// SE vertex neighbors (2) - TRI mesh only
   if (verts_per_face == 3) {
      mapfile << "Color3 setcolor\n";

      GetCoord(total_length - ghost_width, ghost_width, x[0], y[0]);
      GetCoord(total_length - ghost_width - tri_shift_ghst, 0, x[1], y[1]);
      GetCoord(total_length - ghost_width, 0, x[2], y[2]);
      DrawConvexPolygon(verts_per_face, x, y);

      GetCoord(total_length - ghost_width, ghost_width, x[0], y[0]);
      GetCoord(total_length, ghost_width, x[1], y[1]);
      GetCoord(total_length, ghost_width + tri_shift_ghst, x[2], y[2]);
      DrawConvexPolygon(verts_per_face, x, y);
   };

   mapfile << "Color15 setcolor\n";

   GetCoord(total_length - ghost_width, ghost_width, x[0], y[0]);
   GetCoord(total_length - ghost_width, 0, x[1], y[1]);
   GetCoord(total_length, tri_shift_ghst, x[2], y[2]);
   GetCoord(total_length, ghost_width, x[3], y[3]);
   DrawConvexPolygon(verts_per_face, x, y);

// NE or N vertex neighbors (2)
   if (verts_per_face == 3) {
      mapfile << "Color3 setcolor\n";

      GetCoord(total_length - ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[0], y[0]);
      GetCoord(total_length, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[1], y[1]);
      GetCoord(total_length, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width) + tri_shift_ghst, x[2], y[2]);
      DrawConvexPolygon(verts_per_face, x, y);

      GetCoord(total_length - ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[0], y[0]);
      GetCoord(total_length - ghost_width, MaxVertJ(total_length, total_length - ghost_width), x[1], y[1]);
      GetCoord(total_length - ghost_width - tri_shift_ghst, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[2], y[2]);
      DrawConvexPolygon(verts_per_face, x, y);
   };

   mapfile << "Color15 setcolor\n";

   GetCoord(total_length - ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[0], y[0]);
   GetCoord(total_length, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width) + tri_shift_ghst, x[1], y[1]);
   GetCoord(total_length - tri_shift_ghst, MaxVertJ(total_length, total_length) - tri_shift_ghst, x[2], y[2]);
   GetCoord(total_length - ghost_width, MaxVertJ(total_length, total_length), x[3], y[3]);
   DrawConvexPolygon(verts_per_face, x, y);

// NW vertex neighbor
   if (verts_per_face == 4) {
      mapfile << "Color15 setcolor\n";

      GetCoord(square_fill * ghost_width, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[0], y[0]);
      GetCoord(square_fill * ghost_width, MaxVertJ(total_length, total_length), x[1], y[1]);
      GetCoord(0, MaxVertJ(total_length, total_length), x[2], y[2]);
      GetCoord(0, MaxVertJ(base_vert, total_length - (square_fill + 1) * ghost_width, total_length - ghost_width), x[3], y[3]);
      DrawConvexPolygon(verts_per_face, x, y);
   };
};

/*!
\author Vladimir Florinski
\date 06/10/2024
\param[in] data Vertex centered data to print
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::PrintVertContent(double* data)
{
   int data_int, vert;
   double x0, y0, x, y;

// Use black color for text
   mapfile << "Color1 setcolor\n";
   vert = 0;

// Loop over vertices
   for(auto i = 0; i <= total_length; i++) {
      for(auto j = 0; j <= MaxVertJ(total_length, i); j++) {
         GetCoord(i, j, x0, y0);

// Print the data or vertex index
         data_int = (data == nullptr ? vert : (int)data[vert]);
         x = border + x0 + GetShift(data_int);
         y = border + rise + y0 - len[1] * 0.27;

         mapfile << x << " "  << y << " moveto\n";
         mapfile << "(" << data_int << ") show\n";

         vert++;
      };
   };
};

/*!
\author Vladimir Florinski
\date 06/10/2024
\param[in] data Edge centered data to print
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::PrintEdgeContent(double* data)
{
   int data_int, edge;
   double x0, y0, x, y;

// Positioning for edges of different types
   double dx[] = {0.5 * len[0], -0.5 * len[2], 0.5 * len[2]};
   double dy[] = {0.0,           0.4 * len[1], 0.4 * len[1]};

// Use black color for text
   mapfile << "Color1 setcolor\n";
   edge = 0;

   for (auto etype = 0; etype < cardinal_directions; etype++) {
      for (auto i = 0; i <= total_length + edge_dimax[etype]; i++) {
         for (auto j = 0; j <= MaxVertJ(total_length, i) + edge_djmax[etype]; j++) {
            GetCoord(i, j, x0, y0);
      
// Print the data or edge index
            data_int = (data == nullptr ? edge : (int)data[edge]);
            x = border + x0 + GetShift(data_int) + dx[etype];
            y = border + rise + y0 + dy[etype];

            mapfile << x << " "  << y << " moveto\n";
            mapfile << "(" << data_int << ") show\n";

            edge++;
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 06/10/2024
\param[in] data Face centered data to print
*/
template <int verts_per_face>
void DrawableSector<verts_per_face>::PrintFaceContent(double* data)
{
   int data_int, face;
   double x0, y0, x, y;

// Use black color for text
   mapfile << "Color1 setcolor\n";
   face = 0;

// Loop over faces
   for(auto i = 0; i < total_length; i++) {
      for(auto j = 0; j <= MaxFaceJ(total_length, i); j++) {
         GetCoord(i, j / square_fill, x0, y0);

// Print the data or face index
         data_int = (data == nullptr ? face : (int)data[face]);
         
// Shaded TRI face or all QUAD faces
         if ((j % square_fill) || (verts_per_face == 4)) {
            x = border + x0 + GetShift(data_int) + len[6] / 2.0;
            y = border + rise + y0 + len[1] * 0.5;
         }

// Unshaded TRI faces
         else {
            x = border + x0 + GetShift(data_int) + len[0] / 2.0;
            y = border + rise + y0 + len[1] * 0.3;
         };

// Any QAS or shaded TAS face
         mapfile << x << " " << y << " moveto\n";
         mapfile << "(" << data_int << ") show\n";

         face++;
      };
   };
};

template class DrawableSector<3>;
template class DrawableSector<4>;

};
