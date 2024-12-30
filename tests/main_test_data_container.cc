// g++ -I.. -Wall -DGEO_DEBUG -O0 -g -std=c++20 -o main_test_data_container main_test_data_container.cc ../common/data_container.cc

#include <string>
#include <iostream>
#include <iomanip>

#include "common/vectors.hh"
#include "common/data_container.hh"

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i1 = 10, i2;
   double d1 = 12.3, d2;
   MultiIndex mi1(3, 2, 1), mi2;
   GeoVector gv1(3.0, 2.0, 1.0), gv2;
   std::string s1 = "abc", s2;
   std::vector<int> iv1, iv2;
   iv1.push_back(10);
   iv1.push_back(20);
   iv1.push_back(30);
   std::vector<double> dv1, dv2;
   dv1.push_back(1.5);
   dv1.push_back(3.5);
   dv1.push_back(5.5);

   DataContainer container;

// Put the data into the container
   container.Clear();
   container.Insert(i1);
   container.Insert(d1);
   container.Insert(mi1);
   container.Insert(gv1);
   container.Insert(s1);
   container.Insert(iv1);
   container.Insert(dv1);

// Read the data and verify its integrity
   container.Reset();
   container.Read(i2);
   container.Read(d2);
   container.Read(mi2);
   container.Read(gv2);
   container.Read(s2);
   container.Read(iv2);
   container.Read(dv2);

   std::cerr << "Testing int... " << (i2 == i1 ? "PASS" : "FAIL") << std::endl;
   std::cerr << "Testing double... " << (d2 == d1 ? "PASS" : "FAIL") << std::endl;
   std::cerr << "Testing MultiIndex... " << (mi2 == mi1 ? "PASS" : "FAIL") << std::endl;
   std::cerr << "Testing GeoVector... " << (gv2 == gv1 ? "PASS" : "FAIL") << std::endl;
   std::cerr << "Testing std::string... " << (s2 == s1 ? "PASS" : "FAIL") << std::endl;
   std::cerr << "Testing std::vector<int>... " << (iv2 == iv1 ? "PASS" : "FAIL") << std::endl;
};
