// g++ -I.. -Wall -Wno-comment -std=c++20 -DGEO_DEBUG -g -O0 -o main_test_drawable_sector main_test_drawable_sector.cc ../geodesic/drawable_sector.cc ../geodesic/geodesic_sector.cc ../common/vectors.cc
// ./main_test_drawable_sector -v <vert_per_face> -d <div_dif> -g <ghost_width>

#include <config.h>

#ifndef USE_CXXOPTS
#error This software requires cxxopts library to build
#endif

#include <cxxopts.hpp>
#include <common/print_warn.hh>
#include <geodesic/drawable_sector.hh>

using namespace Spectrum;

/*!
\brief parse the command-line parameters 
\author Vladimir Florinski
\date 06/29/2026
\param[in]  argc        Number of command-line arguments plus one
\param[in]  argv        Command-line arguments
\param[out] vert_face   Vertices per face (3 or 4)
\param[out] div_dif     Difference between face and sector divisions
\param[out] ghost_width Width of ghost zone halo
\return True on success, false if some of the parameters were bad 
*/
bool ParseCLOpts(int argc, char** argv, int& vert_face, int& div_dif, int& ghost_width)
{
// Declare the command line options
   cxxopts::Options options("main_test_drawable_sector", "Visualize element numbering and print neighbor information within a sector");
   options.add_options()("v,vert_face", "Vertices per face (3 or 4)", cxxopts::value<int>()->default_value("3"))
                        ("d,div_dif", "Difference between face and sector divisions", cxxopts::value<int>()->default_value("3"))
                        ("g,ghost_width", "Width of ghost zone halo", cxxopts::value<int>()->default_value("1"));

   auto result = options.parse(argc, argv);

   vert_face = result["vert_face"].as<int>();
   div_dif = result["div_dif"].as<int>();
   ghost_width = result["ghost_width"].as<int>();

   if ((vert_face < 3) || (vert_face > 4)) {
      std::cerr << "The \"vert_face\" argument must be either 3 or 4\n";
      return false;
   };
   if (div_dif < 1) {
      std::cerr << "\"div_dif\" must be greater than 1\n";
      return false;
   };

   int width = Pow2(div_dif);
   if (width < min_width_to_ghost * ghost_width) {
      std::cerr << "\"ghost_width\" is too large\n";
      return false;
   };
   return true;
};

int main(int argc, char** argv)
{
   int vert_per_face, div_dif, ghost_width;
   if (!ParseCLOpts(argc, argv, vert_per_face, div_dif, ghost_width)) return 1;

   DrawableSector<3> sector3(Pow2(div_dif), ghost_width);
   DrawableSector<4> sector4(Pow2(div_dif), ghost_width);

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
