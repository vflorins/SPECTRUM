// g++ -I.. -Wall -std=c++20 -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_test_drawable_tesselation main_test_drawable_tesselation.cc ../geodesic/spherical_tesselation.cc ../geodesic/drawable_tesselation.cc ../common/vectors.cc
// ./main_test_drawable_tesselation <conn_division> <draw_division>

#include <common/print_warn.hh>
#include <geodesic/drawable_tesselation.hh>

using namespace Spectrum;

// Edit these to your needs
const PolyType poly_type = PolyType::POLY_HEXAHEDRON;
const int max_division = 6;

int main(int argc, char** argv)
{
   int conn_division = 0;
   int draw_division = 3;
   if (argc > 1) conn_division = atoi(argv[1]);
   if (argc > 2) draw_division = atoi(argv[2]);

   DrawableTesselation<poly_type, max_division> tesselation;
   std::cerr << "Testing the tesselation class for type " << inf_color << tesselation.GetType() << std_color << " base solid\n";

   tesselation.PrintStats();
   for (auto con = 1; con <= 9; con++) {
      if (con != 5) {
         std::cerr << std::endl;
         tesselation.PrintConn(conn_division, con);
      };
   };
   std::cerr << std::endl;
   for (auto div = 0; div <= max_division; div++) {
      std::cerr << "Testing connectivity at division " << div <<": ";
      tesselation.TestConnectivity(div);
   };
   tesselation.DrawGridArcs(draw_division, true, DegToRad(-100.0), DegToRad(-30.0), true);
};
