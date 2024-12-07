// mpicxx -I.. -Wall -Wno-comment -std=c++20 -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_test_drawable_tesselation main_test_drawable_tesselation.cc ../geodesic/spherical_tesselation.cc ../geodesic/drawable_tesselation.cc ../common/vectors.cc

#include "geodesic/drawable_tesselation.hh"

using namespace Spectrum;

//#define POLY_TYPE POLY_TETRAHEDRON
//#define POLY_TYPE POLY_HEXAHEDRON
//#define POLY_TYPE POLY_OCTAHEDRON
//#define POLY_TYPE POLY_DODECAHEDRON
#define POLY_TYPE POLY_ICOSAHEDRON
#define MAX_DIVISION 5

int main(int argc, char** argv)
{
   std::cerr << "Testing the tesselation class for type " << POLY_TYPE << " base solid\n";

   DrawableTesselation<POLY_TYPE, MAX_DIVISION> tess;
   tess.PrintStats();
   tess.PrintConn(1, 1);
   tess.TestConnectivity(MAX_DIVISION);
   tess.DrawGridArcs(4, true, DegToRad(-100.0), DegToRad(-30.0), true);
};

