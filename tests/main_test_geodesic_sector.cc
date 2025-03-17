#include <iostream>
#include <geodesic/geodesic_sector.hh>

using namespace Spectrum;

constexpr int verts_per_face = 3;
constexpr int len = 8;

int main(int argc, char *argv[])
{
   GeodesicSector<verts_per_face> sector;
   sector.SetDimensions(len, 2, false);
   std::cerr << std::endl;
   sector.PrintAddresses(1);
   std::cerr << std::endl;
   sector.PrintAddresses(2);
   std::cerr << std::endl;
   sector.PrintAddresses(3);
   std::cerr << std::endl;
   sector.PrintConn(1);
   std::cerr << std::endl;
   sector.PrintConn(2);
   std::cerr << std::endl;
   sector.PrintConn(3);
   std::cerr << std::endl;
   sector.PrintConn(4);
   std::cerr << std::endl;
   sector.PrintConn(6);
   std::cerr << std::endl;
   sector.PrintConn(7);
   std::cerr << std::endl;
   sector.PrintConn(8);
   std::cerr << std::endl;
   sector.PrintConn(9);
   std::cerr << std::endl;
};

