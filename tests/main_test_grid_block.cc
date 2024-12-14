// g++ -I.. -Wall -Wno-comment -std=c++20 -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_test_grid_block main_test_grid_block.cc ../geometry/distance_map.cc ../geodesic/spherical_tesselation.cc ../geodesic/requestable_tesselation.cc ../geodesic/traversable_tesselation.cc ../geodesic/geodesic_sector.cc ../geodesic/grid_block.cc ../common/params.cc ../common/data_container.cc ../common/vectors.cc -lsiloh5  -lgsl -lgslcblas

#include <silo.h>
#include "geodesic/grid_block.hh"
#include "geodesic/traversable_tesselation.hh"

using namespace Spectrum;

#define POLY_TYPE POLY_HEXAHEDRON
//#define POLY_TYPE POLY_ICOSAHEDRON
#define MAX_DIVISION 6
#define verts_per_face 4

const double rmin = 30.0;
const double rmed = 300.0;
const double rmax = 1000.0;
const double power_law = 10.0;
const double fraction = 0.4;
const double logy = 10.0;

const int face_division = 6;
const int sect_division = 0;
const int ghost_width = 2;

const int slab_height = 720;
const int n_slabs = 1;
const int ghost_height = 2;
const int sector = 0;
const int slab = 0;

int main(int argc, char *argv[])
{
   bool corners[verts_per_face];
   bool borders[2];
   borders[0] = (slab == 0 ? true : false);
   borders[1] = (slab == n_slabs - 1 ? true : false);

// Create a Tesselation object
   TraversableTesselation<POLY_TYPE, MAX_DIVISION> tess;

// Set up a radial map
   DataContainer container;
   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
//   container.Insert(power_law);
   container.Insert(rmed);
   container.Insert(fraction);
   container.Insert(logy);

   std::shared_ptr<DistanceBase> dist_map = std::make_shared<DistanceLinearExp>();
   dist_map->SetupObject(container);

// Create one grid block
   GridBlock<verts_per_face> block(Pow2(face_division - sect_division), ghost_width, slab_height, ghost_height);
   block.SetIndex(0);

   int* face_index = new int[block.TotalFaces()];
   int* vert_index = new int[block.TotalVerts()];
   GeoVector* vcart = new GeoVector[block.TotalVerts()];

   tess.GetAllInsideFaceNative(sect_division, sector, face_division, ghost_width, face_index, vert_index, corners);

   tess.FillVertCoordArrays(block.TotalVerts(), vert_index, vcart);
   block.AssociateMesh(slab / (double)n_slabs, (slab + 1.0) / (double)n_slabs, corners, borders, vcart, dist_map);
   block.PrintStats();
//   block.DrawZone(3, 57, DegToRad(0.0), DegToRad(0.0));
//   block.DrawZone(4, 57, DegToRad(0.0), DegToRad(0.0));

   DBfile* silofile = DBCreate("onegridlock.silo", DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
   block.WriteSilo(silofile, false);
   DBClose(silofile);

   delete[] face_index;
   delete[] vert_index;
   delete[] vcart;

   return 0;
};

