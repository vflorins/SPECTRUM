// mpicxx -I.. -Wall -Wno-comment -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_onestenciledblock main_onestenciledblock.cc ../geometry/distance_map.cc ../geodesic/spherical_tesselation.cc ../geodesic/requestable_tesselation.cc ../geodesic/traversable_tesselation.cc ../geodesic/geodesic_sector.cc ../geodesic/grid_block.cc ../geodesic/stenciled_block.cc ../common/params.cc ../common/data_container.cc ../common/vectors.cc ../common/polynomial.cc -lsiloh5

#include <silo.h>
#include "geodesic/stenciled_block.hh"
#include "geodesic/traversable_tesselation.hh"

using namespace Spectrum;

constexpr PolyType poly_type = POLY_ICOSAHEDRON;
constexpr int verts_per_face = 3;

const double rmin = 1.0;
const double rmax = 100.0;
const double pwl = 2.0;

const int face_division = 4;
const int sect_division = 1;
const int ghost_width = 2;

const int slab_height = 8;
const int n_slabs = 30;
const int ghost_height = 2;
const int sector = 1;
const int slab = 1;

int main(int argc, char *argv[])
{
   bool corners[3];
   bool borders[2] = {false, false};

// Create a Tesselation object
   TraversableTesselation<poly_type, face_division> tess;

// Set up a radial map
   DataContainer container;
   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
   container.Insert(pwl);

   std::shared_ptr<DistanceBase> dist_map = std::make_shared<DistancePowerLaw>();
   dist_map->SetupObject(container);

// Create one stenciled block
   StenciledBlock<verts_per_face> block(Pow2(face_division - sect_division), ghost_width, slab_height, ghost_height);

   int *face_index = new int[block.TotalFaces()];
   int *vert_index = new int[block.TotalVerts()];
   GeoVector* vcart = new GeoVector[block.TotalVerts()];

   tess.GetAllInsideFaceNative(sect_division, sector, face_division, ghost_width, face_index, vert_index, corners);
   tess.FillVertCoordArrays(block.TotalVerts(), vert_index, vcart);
   block.AssociateMesh(0, slab / (double)n_slabs, (slab + 1.0) / (double)n_slabs, corners, borders, vcart, dist_map);
   block.PrintStats();
   block.DrawStencil(3, 57, 1, DegToRad(0.0), DegToRad(0.0));

   DBfile* silofile = DBCreate("onegridlock.silo", DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
   block.WriteSilo(silofile, false);
   DBClose(silofile);

   delete[] face_index;
   delete[] vert_index;
   delete[] vcart;

   return 0;
};

