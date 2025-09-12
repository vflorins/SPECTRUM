/*!
\file spherical_tesselation.hh
\brief Declares SphericalTesselation class, a recursive polygonal partitioning of a sphere
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPHERICAL_TESSELATION_HH
#define SPECTRUM_SPHERICAL_TESSELATION_HH

#include <iostream>
#include <cstdint>

#include "common/vectors.hh"
#include "geodesic/polyhedron.hh"

namespace Spectrum {

#define N_TESS_ERRORS 6

using TERR_TYPE = uint8_t;

//! No error
#define TESERR_NOERR 0x00

//! Bad index entry error
#define TESERR_INPUT 0x01

//! Bad index entry error
#define TESERR_INDEX 0x02

//! Overfilled connectivity table error
#define TESERR_OVERF 0x04

//! Underfilled connectivity table error
#define TESERR_UNDER 0x08

//! Connectivity mismatch in two tables
#define TESERR_MISMT 0x10

//! All tesselation errors
constexpr TERR_TYPE tess_errors[] = {TESERR_NOERR, TESERR_INPUT, TESERR_INDEX, TESERR_OVERF, TESERR_UNDER, TESERR_MISMT};

//! All tesselation eror messages
const std::string tess_error_msg[] = {"",
                                      "Bad input",
                                      "Bad index in connectivity table",
                                      "Overfilled connectivity table",
                                      "Underfilled connectivity table",
                                      "Mismatch between connectivity tables"};


//! Computes a general reverse connectivity table
TERR_TYPE BuildReverse(int n_nodes1, int n_nbrs1, int n_sing1, int n_nbrs1s, int start_sing1, const int* const* conn_12,
                       int n_nodes2, int n_nbrs2, int n_sing2, int n_nbrs2s, int start_sing2,       int* const* conn_21);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A simple class to store error information for use in exceptions
\author Vladimir Florinski
*/
struct TessError
{

//! Calling function
   std::string cid;

//! Division
   int div;

//! Error
   TERR_TYPE err;

//! Line in the code where the error occurred
   int lin;

//! Constructor with arguments
   TessError(const std::string& caller, int division, TERR_TYPE error, int line);

//! Decode and print error mesages
   void PrintErrors(void) const;
};

/*!
\author Vladimir Florinski
\date 04/15/2020
\param[in] caller   Name of the function throwing the exception
\param[in] division Division
\param[in] error    Error
*/
inline TessError::TessError(const std::string& caller, int division, TERR_TYPE error, int line)
                : cid(caller),
                  div(division),
                  err(error),
                  lin(line)
{
};

/*!
\author Vladimir Florinski
\date 04/15/2020
*/
inline void TessError::PrintErrors(void) const
{
   for (auto ierr = 0; ierr < N_TESS_ERRORS; ierr++) {
      if (err & tess_errors[ierr]) {
         std::cerr << "Tesselation::" << cid << ": " << tess_error_msg[ierr] << " at division " << div << ", line " << lin << "\n";
      };
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SphericalTesselation class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A class describing a nested polygonal tesselation of a sphere
\author Vladimir Florinski

Tesselation is a partitioning of the surface of a sphere into spherical polygons. This class defines data structures to store corrdinates and methods to compute mesh connectivity.
*/
template <PolyType poly_type, int max_division>
class SphericalTesselation : public Polyhedron<poly_type>
{
protected:

//! The number of vertices in each division
   int nverts[max_division + 1];

//! The number of edges in each division
   int nedges[max_division + 1];

//! The number of t-faces in each division
   int nfaces[max_division + 1];

//! Number of vertices per face; could be different for division 0
   int verts_per_face[max_division + 1];

//! Number of vertex/edge/face neighbors per vertex; could be different for division 0; could be different for singular vertices
   int edges_per_vert[max_division + 1];

//! Into how many faces the parent face is divided; could be different for division 0
   int children_per_face[max_division + 1];

//! New vertices are added for each parent _edge_
   bool newverts_at_edge[max_division + 1];

//! New vertices are added for each parent _face_
   bool newverts_at_face[max_division + 1];

//! Cartesian coordinates of vertices on a unit sphere
   GeoVector* vert_cart = nullptr;

//! Vertex-vertex connectivity array (not ordered)
   int** vv_con[max_division + 1] = {nullptr};

//! Vertex-edge connectivity array (not ordered)
   int** ve_con[max_division + 1] = {nullptr};

//! Vertex-face connectivity array (ordered counter-clockwise)
   int** vf_con[max_division + 1] = {nullptr};

//! Edge-vertex connectivity array (not ordered)
   int** ev_con[max_division + 1] = {nullptr};

//! Edge-face connectivity array (not ordered)
   int** ef_con[max_division + 1] = {nullptr};
   
//! Face-vertex connectivity array (ordered counter-clockwise)
   int** fv_con[max_division + 1] = {nullptr};

//! Face-edge connectivity array (ordered counter-clockwise)
   int** fe_con[max_division + 1] = {nullptr};

//! Face-face connectivity array (ordered counter-clockwise)
   int** ff_con[max_division + 1] = {nullptr};

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Actual number of vertex neighbors
   SPECTRUM_DEVICE_FUNC int NVertNbrs(int div, int vert) const;

//! Return an FV entry CC from the given vertex.
   SPECTRUM_DEVICE_FUNC int VertCC(int div, int face, int vert, int rot) const;

//! Return an FE entry CC from the given edge.
   SPECTRUM_DEVICE_FUNC int EdgeCC(int div, int face, int edge, int rot) const;

//! Find two faces shared by two vertices.
   SPECTRUM_DEVICE_FUNC void Match2vf(int div, int vert1, int vert2, int& face1, int& face2) const;

//! Find the edge connecting two vertices.
   SPECTRUM_DEVICE_FUNC int Match2ve(int div, int vert1, int vert2) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Insert vertices of a daughter division and compute VV (CC-ordered).
   void RefineVert(int div);

//! Connect a set of vertices with edges and compute EV.
   void AddEdges(int div);

//! Insert faces of a daughter division and compute FV (CC-ordered).
   void RefineFace(int div);

//! Compute VE (CC-ordered and synchronized with VV).
   void VertEdgeConn(int div);

//! Compute VF (CC-ordered and synchronized with VV, VE).
   void VertFaceConn(int div);

//! Compute EF (synchronizd with EV).
   void EdgeFaceConn(int div);

//! Compute FE (CC-ordered and synchronized with FV).
   void FaceEdgeConn(int div);

//! Compute FF (CC-ordered and synchronized with FV. FE).
   void FaceFaceConn(int div);

//! Memory allocator
   void AllocateStorage(void);

//! Compute all grid coordinates and connectivity
   void ComputeAll(void);

//! Memory de-allocator
   void FreeStorage(void);

public:

//! Default constructor
   SphericalTesselation(void);

//! Copy constructor
   SphericalTesselation(const SphericalTesselation& other) = delete;

//! Destructor
   ~SphericalTesselation();

//! Return number of vertices in a given division.
   SPECTRUM_DEVICE_FUNC auto NVerts(int div) const {return nverts[div];};

//! Return number of edges in a given division.
   SPECTRUM_DEVICE_FUNC auto NEdges(int div) const {return nedges[div];};

//! Return number of faces in a given division.
   SPECTRUM_DEVICE_FUNC auto NFaces(int div) const {return nfaces[div];};
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SphericalTesselation inlined methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/06/2020
v\param[in] div  Division
\\param[in] vert Vertex
return Number of neighbor vertices (same as edges and faces)
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC inline int SphericalTesselation<poly_type, max_division>::NVertNbrs(int div, int vert) const
{
   if (poly_type == POLY_DODECAHEDRON) {
      if (!div) return edges_per_vert[0];
      else if ((vert >= nverts[0]) && (vert < nverts[0] + nfaces[0])) return verts_per_face[0];
      else return edges_per_vert[div];
   }
   else if (vert < nverts[0]) return edges_per_vert[0];
   else return edges_per_vert[div];
};

/*!
\author Vladimir Florinski
\date 04/05/2020
\param[in] div  Division
\param[in] face Face
\param[in] vert Vertex that belongs to "face"
\param[in] rot  Amount of rotation in the CC direction
\return Another vertex that belongs to "face"
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC inline int SphericalTesselation<poly_type, max_division>::VertCC(int div, int face, int vert, int rot) const
{
   auto iv = InList(verts_per_face[div], fv_con[div][face], vert);
   if (iv == -1) return -1;
   else return fv_con[div][face][(iv + rot + verts_per_face[div]) % verts_per_face[div]];
};

/*!
\author Vladimir Florinski
\date 04/05/2020
\param[in] div  Division
\param[in] face Face
\param[in] vert Edge that belongs to "face"
\param[in] rot  Amount of rotation in the CC direction
\return Another edge that belongs to "face"
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC inline int SphericalTesselation<poly_type, max_division>::EdgeCC(int div, int face, int edge, int rot) const
{
   auto ie = InList(verts_per_face[div], fe_con[div][face], edge);
   if (ie == -1) return -1;
   else return fe_con[div][face][(ie + rot + verts_per_face[div]) % verts_per_face[div]];
};

/*!
\author Vladimir Florinski
\date 04/06/2020
\param[in]  div   Division
\param[in]  vert1 First vertex
\param[in]  vert2 Second vertex
\param[out] face1 First  common face (-1 if vertices are not connected)
\param[out] face2 Second common face (-1 if vertices are not connected)
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC inline void SphericalTesselation<poly_type, max_division>::Match2vf(int div, int vert1, int vert2, int& face1, int& face2) const
{
   bool found_two = false;
   int face3, it1 = 0;
   face1 = face2 = -1;

// If "vert1" is the same as "vert2", this will pick the first two faces that share this vertex.
   while ((it1 < NVertNbrs(div, vert1)) && !found_two) {
      face3 = vf_con[div][vert1][it1];
      if (InList(NVertNbrs(div, vert2), vf_con[div][vert2], face3) != -1) {
         if (face1 == -1) face1 = face3;
         else {
            face2 = face3;
            found_two = true;
         };
      };
      it1++;
   };
};

/*!
\author Vladimir Florinski
\date 04/06/2020
\param[in] div   Division
\param[in] vert1 First vertex
\param[in] vert2 Second vertex
\return Edge connecting "vert1" and "vert2" (-1 if vertices are not connected)
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC inline int SphericalTesselation<poly_type, max_division>::Match2ve(int div, int vert1, int vert2) const
{
   bool found_one = false;
   int edge, ie1 = 0;

// If "vert1" is the same as "vert2", this will pick its first edge.
   while ((ie1 < NVertNbrs(div, vert1)) && !found_one) {
      edge = ve_con[div][vert1][ie1];
      if (InList(NVertNbrs(div, vert2), ve_con[div][vert2], edge) != -1) return edge;
      ie1++;
   };

   return -1;
};

};

#endif