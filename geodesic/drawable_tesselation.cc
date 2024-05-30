/*!
\file drawable_tesselation.cc
\brief Implements the DrawableTesselation class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <iomanip>
#include <fstream>
#include "geodesic/drawable_tesselation.hh"
#include "common/print_warn.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 02/20/2020
\param[in] div    The division to draw
\param[in] opaque Do not draw lines that are behind the sphere
\param[in] rot_z  First rotation angle about the z-axis
\param[in] rot_x  Second rotation angle about the x-axis
\param[in] smooth Draw curved (true) or straight (false) edges

The output is a series of curves that when visualized on an XY plot represent the great circle arcs that are edges of the tesselation.
*/
template <PolyType poly_type, int max_division>
void DrawableTesselation<poly_type, max_division>::DrawGridArcs(int div, bool opaque, double rot_z, double rot_x, bool smooth) const
{
   int edge, ipt, nint, n_drawn;
   double ddist;

// A drawing would be too busy for division values greated than 6
   if(div > 6) {
      std::cerr << "Tesselation: Division too high for a plot\n";
      return;
   };

// Precompute trig functions of the rotation angles
   double snrz = sin(rot_z);
   double csrz = cos(rot_z);
   double snrx = sin(rot_x);
   double csrx = cos(rot_x);

// Number of segments used to draw each arc
   if(smooth) nint = Pow2(6 - div);
   else nint = 1;

   GeoVector v, v1, v2, v3;
   std::ofstream meshfile;
   meshfile.open(("mesh_" + poly_names[poly_type] + "_" + std::to_string(div) + ".out").c_str(), std::ofstream::out);
   meshfile << "# Drawing edges at division " << div << std::endl;

// Loop over edges
   for(edge = 0; edge < nedges[div]; edge++) {
      v1 = vert_cart[ev_con[div][edge][0]];
      v2 = vert_cart[ev_con[div][edge][1]];
      v3 = (v1 ^ v2) ^ v1;
      v3.Normalize();

// Compute the angle between the points and divide it into equal segments
      ddist = acos(v1 * v2) / (double)nint;

// Draw the points to form an edge
      n_drawn = 0;
      for(ipt = 0; ipt <= nint; ipt++) {

// Interior point vectors are linear combinations of "v1" and "v3"
         v = cos(ipt * ddist) * v1 + sin(ipt * ddist) * v3;

// Rotate the view
         v.Rotate(gv_nz, snrz, csrz);
         v.Rotate(gv_nx, snrx, csrx);

// Print the projection on the xz plane
         if(v[1] < 0.0 || !opaque) {
            meshfile << std::setw(12) << std::setprecision(5) << v[0]
                     << std::setw(12) << std::setprecision(5) << v[2]
                     << std::endl;
            n_drawn++;
         };
      };

// The ampersand is required by ggrace (print only if some points in the line are visible).
      if(n_drawn) meshfile << "&\n";
   };

   meshfile.close();
};

/*!
\author Vladimir Florinski
\date 08/30/2019
*/
template <PolyType poly_type, int max_division>
void DrawableTesselation<poly_type, max_division>::PrintStats(void) const
{
   int div, face, edge, iv;
   double vertex_angle, vertex_angle_min, vertex_angle_max, vertex_angle_mean;
   double edge_length, edge_length_min, edge_length_max, edge_length_mean;
   double face_area, face_area_min, face_area_max, face_area_mean;
   
// Memory for latitude/longitude and 3 Cartesian components
   int memtot = 5 * nverts[max_division] * SZDBL;

   std::cerr << std::endl;
   std::cerr << "╭─────┬─────────┬─────────┬─────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────╮\n";
   std::cerr << "│ Div │ Vertices│   Edges │   Faces │ Min edge │ Avg edge │ Max edge │ Min angl │ Avg angl │ Max angl │ Min area │ Avg area │ Max area │\n";
   std::cerr << "├─────┼─────────┼─────────┼─────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n";

// Compute mesh statistics
   for(div = 0; div <= max_division; div++) {

// Memory for connectivity arrays
      memtot += 3 * nverts[div] * edges_per_vert[div] * SZINT
              + 2 * nedges[div] * 2 * SZINT
              + 3 * nfaces[div] * verts_per_face[div] * SZINT;

// Vertex angles
      vertex_angle_min = M_2PI;
      vertex_angle_max = 0.0;
      vertex_angle_mean = 0.0;
      for(face = 0; face < nfaces[div]; face++) {
         for(iv = 0; iv < verts_per_face[div]; iv++) {
            vertex_angle = acos(VertexAngle(vert_cart[fv_con[div][face][iv]],
                                            vert_cart[fv_con[div][face][(iv + 1) % verts_per_face[div]]],
                                            vert_cart[fv_con[div][face][(iv + 2) % verts_per_face[div]]]));
            vertex_angle_min = fmin(vertex_angle_min, vertex_angle);
            vertex_angle_max = fmax(vertex_angle_max, vertex_angle);
            vertex_angle_mean += vertex_angle;
         };
      };
      vertex_angle_mean /= verts_per_face[div] * nfaces[div];

// Edge lengths
      edge_length_min = M_2PI;
      edge_length_max = 0.0;
      edge_length_mean = 0.0;
      for(edge = 0; edge < nedges[div]; edge++) {
         edge_length = acos(vert_cart[ev_con[div][edge][0]] * vert_cart[ev_con[div][edge][1]]);
         edge_length_min = fmin(edge_length_min, edge_length);
         edge_length_max = fmax(edge_length_max, edge_length);
         edge_length_mean += edge_length;
      };
      edge_length_mean /= nedges[div];

// Face areas
      face_area_min = M_4PI;
      face_area_max = 0.0;
      face_area_mean = 0.0;
      for(face = 0; face < nfaces[div]; face++) {

// Divide the polygon into triangles and sum up their partial areas
         face_area = 0.0;
         for(iv = 2; iv < verts_per_face[div]; iv++) {
            face_area += SphTriArea(vert_cart[fv_con[div][face][0]], vert_cart[fv_con[div][face][iv - 1]], vert_cart[fv_con[div][face][iv]]);
         };
         face_area_min = fmin(face_area_min, face_area);
         face_area_max = fmax(face_area_max, face_area);
         face_area_mean += face_area;
      };

      face_area_mean /= nfaces[div];

// Print the numbers
      std::cerr << "│"
                << std::setw(3) << div << "  │"
                << std::setw(8) << nverts[div] << " │"
                << std::setw(8) << nedges[div] << " │"
                << std::setw(8) << nfaces[div] << " │"
                << std::setprecision(4)
                << std::fixed << std::setw(8) << RadToDeg(edge_length_min)   << "° │"
                << std::fixed << std::setw(8) << RadToDeg(edge_length_mean)  << "° │"
                << std::fixed << std::setw(8) << RadToDeg(edge_length_max)   << "° │"
                << std::fixed << std::setw(8) << RadToDeg(vertex_angle_min)  << "° │"
                << std::fixed << std::setw(8) << RadToDeg(vertex_angle_mean) << "° │"
                << std::fixed << std::setw(8) << RadToDeg(vertex_angle_max)  << "° │"
                << std::setprecision(7)
                << std::fixed << std::setw(10) << face_area_min  * 10.0 << "│"
                << std::fixed << std::setw(10) << face_area_mean * 10.0 << "│"
                << std::fixed << std::setw(10) << face_area_max  * 10.0 << "│"
                << std::endl;
   };
   std::cerr << "╰─────┴─────────┴─────────┴─────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────╯\n";
   std::cerr << "Memory used for tesselation is " << memtot << " bytes\n\n";
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 04/15/2020
\param[in] div Division
*/
template <PolyType poly_type, int max_division>
void DrawableTesselation<poly_type, max_division>::TestConnectivity(int div) const
try {

   static const std::string callerID = "TestConnectivity";
   if((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int vert, vert1, edge, edge1, face, face1, face2, face3, iv, iv1, ie, it;

   for(vert = 0; vert < nverts[div]; vert++) {

// Check whether "vert" is in the vertex neighbor's ("vert1") VV.
      for(iv = 0; iv < NVertNbrs(div, vert); iv++) {
         vert1 = vv_con[div][vert][iv];
         if((vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         iv1 = InList(NVertNbrs(div, vert1), vv_con[div][vert1], vert);
         if(iv1 == -1) throw TessError(callerID, div, TESERR_MISMT, __LINE__);
      };

// Check whether "vert" is in the edge neighbor's ("edge") EV and whether the other vertex of "edge" is the same as VV of "vert" with the same neighbor index.
      for(ie = 0; ie < NVertNbrs(div, vert); ie++) {
         edge = ve_con[div][vert][ie];
         if((edge < 0) || (edge >= nedges[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         iv1 = InList(2, ev_con[div][edge], vert);
         if(iv1 == -1) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         vert1 = ev_con[div][edge][(iv1 + 1) % 2];
         if((vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if(vert1 != vv_con[div][vert][ie]) TessError(callerID, div, TESERR_MISMT, __LINE__);
      };

// Check whether "vert" is in the face neighbor's ("face") FV and whether the next vertex of "face" is the same the previous vertex in VV of "vert".
      for(it = 0; it < NVertNbrs(div, vert); it++) {
         face = vf_con[div][vert][it];
         if((face < 0) || (face >= nfaces[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         iv1 = InList(verts_per_face[div], fv_con[div][face], vert);
         if(iv1 == -1) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         vert1 = fv_con[div][face][(iv1 + 1) % verts_per_face[div]];
         if((vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if(vert1 != vv_con[div][vert][(it - 1 + NVertNbrs(div, vert)) % NVertNbrs(div, vert)]) TessError(callerID, div, TESERR_MISMT, __LINE__);
      };

// Check whether face neighbors of "vert" are in the right place in their edge neighbor's EF, whether edge neighbors of "vert" are in the right place in the face neighbor's FE, and whether face neighbors of "vert" are in the right place in the face neighbor's FF.
      for(ie = 0; ie < NVertNbrs(div, vert); ie++) {
         edge = ve_con[div][vert][ie];
         iv = InList(2, ev_con[div][edge], vert);
         face1 = vf_con[div][vert][ie];
         face2 = vf_con[div][vert][(ie + 1) % NVertNbrs(div, vert)];
         if(ef_con[div][edge][0] != (iv ? face1 : face2)) throw TessError(callerID, div, TESERR_MISMT, __LINE__);
         if(ef_con[div][edge][1] != (iv ? face2 : face1)) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         iv = InList(verts_per_face[div], fv_con[div][face1], vert);
         edge1 = fe_con[div][face1][(iv - 1 + verts_per_face[div]) % verts_per_face[div]];
         if((edge1 < 0) || (edge1 >= nedges[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if(edge1 != edge) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         face3 = ff_con[div][face1][(iv - 1 + verts_per_face[div]) % verts_per_face[div]];
         if((face3 < 0) || (face3 >= nfaces[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if(face3 != face2) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         iv = InList(verts_per_face[div], fv_con[div][face2], vert);
         edge1 = fe_con[div][face2][iv];
         if((edge1 < 0) || (edge1 >= nedges[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if(edge1 != edge) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         face3 = ff_con[div][face2][iv];
         if((face3 < 0) || (face3 >= nfaces[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if(face3 != face1) throw TessError(callerID, div, TESERR_MISMT, __LINE__);
      };
   };
   std::cerr << "Connectivity test successful\n";
}

catch(const TessError& err) {
   err.PrintErrors();
   return;
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] div  The division
\param[in] type The connectivity array to print

The type is one of the folloing: 1 is vertex-vertex, 2 is vertex-edge, 3 is vertex-face, 4 is edge-vertex, 6 is edge-face, 7 is face-vertex, 8 is face-edge,9 is face-face.
*/
template <PolyType poly_type, int max_division>
void DrawableTesselation<poly_type, max_division>::PrintConn(int div, int type) const
{
   switch(type) {

   case 1:
      std::cerr << "Printing vert-vert connectivity at level " << div << std::endl;
      PrintConnectivity(nverts[div], edges_per_vert[div], 0, vv_con[div]);
      break;

   case 2:
      std::cerr << "Printing vert-edge connectivity at level " << div << std::endl;
      PrintConnectivity(nverts[div], edges_per_vert[div], 0, ve_con[div]);
      break;

   case 3:
      std::cerr << "Printing vert-face connectivity at level " << div << std::endl;
      PrintConnectivity(nverts[div], edges_per_vert[div], 0, vf_con[div]);
      break;

   case 4:
      std::cerr << "Printing edge-vert connectivity at level " << div << std::endl;
      PrintConnectivity(nedges[div], 2, 0, ev_con[div]);
      break;

   case 6:
      std::cerr << "Printing edge-face connectivity at level " << div << std::endl;
      PrintConnectivity(nedges[div], 2, 0, ef_con[div]);
      break;

   case 7:
      std::cerr << "Printing face-vert connectivity at level " << div << std::endl;
      PrintConnectivity(nfaces[div], verts_per_face[div], 0, fv_con[div]);
      break;

   case 8:
      std::cerr << "Printing face-edge connectivity at level " << div << std::endl;
      PrintConnectivity(nfaces[div], verts_per_face[div], 0, fe_con[div]);
      break;

   case 9:
      std::cerr << "Printing face-face connectivity at level " << div << std::endl;
      PrintConnectivity(nfaces[div], verts_per_face[div], 0, ff_con[div]);
      break;
   };
};

#endif

template class DrawableTesselation<POLY_TETRAHEDRON, 3>;
template class DrawableTesselation<POLY_HEXAHEDRON, 3>;
template class DrawableTesselation<POLY_OCTAHEDRON, 3>;
template class DrawableTesselation<POLY_DODECAHEDRON, 3>;
template class DrawableTesselation<POLY_ICOSAHEDRON, 3>;

};
