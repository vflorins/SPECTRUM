#include <iostream>
#include <iomanip>
#include "geodesic/drawable_tesselation.hh"

using namespace Spectrum;

constexpr PolyType poly_type = POLY_ICOSAHEDRON;
constexpr int max_division = 3;

int main(int argc, char *argv[])
{
   std::cerr << "Testing the tesselation class for type " << poly_type << " base solid\n";

   DrawableTesselation<poly_type, max_division> tess;
   tess.PrintStats();
   tess.PrintConn(1, 1);
   tess.TestConnectivity(max_division);
   tess.DrawGridArcs(0, true, DegToRad(-100.0), DegToRad(-30.0), true);

   return 0;
};


