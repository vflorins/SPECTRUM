// g++ -I.. -Wall -std=c++20 -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_test_drawable_tesselation main_test_drawable_tesselation.cc ../geodesic/spherical_tesselation.cc ../geodesic/drawable_tesselation.cc ../common/vectors.cc
// ./main_test_drawable_tesselation -c <conn_division> -d <draw_division>

#include <config.h>

#ifndef USE_CXXOPTS
#error This software requires cxxopts library to build
#endif

#include <cxxopts.hpp>
#include <common/print_warn.hh>
#include <geodesic/drawable_tesselation.hh>

using namespace Spectrum;

// Edit these to your needs
const PolyType poly_type = PolyType::POLY_ICOSAHEDRON;
const int max_division = 6;

/*!
\brief parse the command-line parameters 
\author Vladimir Florinski
\date 06/26/2026
\param[in]  argc     Number of command-line arguments plus one
\param[in]  argv     Command-line arguments
\param[out] conn_div Division for which connectivity will be printed
\param[out] draw_div Division for which images will be generated
\return True on success, false if some of the parameters were bad 
*/
bool ParseCLOpts(int argc, char** argv, int& conn_div, int& draw_div)
{
// Declare the command line options
   cxxopts::Options options("main_test_drawable_tesselation", "Test the tesselation class and print diagnostic information");
   options.add_options()("c,conn_div", "Division for which connectivity will be printed", cxxopts::value<int>()->default_value("0"))
                        ("d,draw_div", "Division for which images will be generated", cxxopts::value<int>()->default_value("0"));

   auto result = options.parse(argc, argv);

   conn_div = result["conn_div"].as<int>();
   draw_div = result["draw_div"].as<int>();

   if ((conn_div < 0) || (conn_div > max_division)) {
      std::cerr << "The \"conn_div\" argument must be between 0 and " << max_division << std::endl;
      return false;
   };
   if ((draw_div < 0) || (draw_div > max_division)) {
      std::cerr << "The \"draw_div\" argument must be between 0 and " << max_division << std::endl;
      return false;
   };

   return true;
};

int main(int argc, char** argv)
{
   int conn_division, draw_division;
   if (!ParseCLOpts(argc, argv, conn_division, draw_division)) return 1;

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
