// g++ -I.. -Wall -Wno-comment -std=c++20 -DGEO_DEBUG -g -O0 -o main_test_drawable_sector main_test_drawable_sector.cc ../geodesic/drawable_sector.cc ../geodesic/geodesic_sector.cc ../common/vectors.cc
// ./main_test_drawable_sector <vert_per_face> <face_division> <sect_division> <ghost_width>

#include <common/definitions.hh>
#include <geodesic/drawable_sector.hh>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int vert_per_face;
   int face_division = 4;
   int sect_division = 1;
   int ghost_width = 3;

   if (argc > 1) vert_per_face = atoi(argv[1]);
   if ((vert_per_face < 3) || (vert_per_face > 4)) vert_per_face = 3;
   if (argc > 2) face_division = atoi(argv[2]);
   if (argc > 3) sect_division = atoi(argv[3]);
   if (argc > 4) ghost_width = atoi(argv[4]);

   DrawableSector<3> sector3(Pow2(face_division - sect_division), ghost_width);
   DrawableSector<4> sector4(Pow2(face_division - sect_division), ghost_width);

   for (auto type = 1; type <= 3; type++) {
      std::cerr << std::endl;
      if (vert_per_face == 3) sector3.PrintAddresses(type);
      else sector4.PrintAddresses(type);
   };
   for (auto type = 1; type <= 9; type++) {
      if (type != 5) {
         std::cerr << std::endl;
         if (vert_per_face == 3) sector3.PrintConn(type);
         else sector4.PrintConn(type);
      };
   };
   for (auto type = 1; type <= 3; type++) {
      std::cerr << std::endl;
      if (vert_per_face == 3) sector3.PrintMask(type);
      else sector4.PrintMask(type);
   };
   for (auto type = 1; type <= 3; type++) {
      std::cerr << std::endl;
      if (vert_per_face == 3) sector3.Draw(type);
      else sector4.Draw(type);
   };
};
