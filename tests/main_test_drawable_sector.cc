// g++ -I.. -Wall -Wno-comment -std=c++20 -DGEO_DEBUG -g -O0 -o main_test_drawable_sector main_test_drawable_sector.cc ../geodesic/drawable_sector.cc ../geodesic/geodesic_sector.cc ../common/vectors.cc

#include "common/definitions.hh"
#include "geodesic/drawable_sector.hh"

using namespace Spectrum;

const int face_division = 4;
const int sect_division = 1;
const int ghost_width = 3;

int main(int argc, char *argv[])
{
   DrawableSector<3> sector3(Pow2(face_division - sect_division), ghost_width);
   sector3.Draw(1);
   sector3.Draw(2);
   sector3.Draw(3);

   DrawableSector<4> sector4(Pow2(face_division - sect_division), ghost_width);
   sector4.Draw(1);
   sector4.Draw(2);
   sector4.Draw(3);

   return 0;
};

