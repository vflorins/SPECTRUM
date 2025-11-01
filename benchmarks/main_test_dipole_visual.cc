#include "src/background_dipole.hh"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{
   BackgroundDipole background;

   SpatialData spdata;
   double t = 0.0;
   int i,j,k;
   GeoVector pos, mom = gv_ones;

   spdata._mask = BACKGROUND_ALL;

   DataContainer container;
   container.Clear();

// Initial time
   double t0 = 0.0;
   container.Insert(t0);

// Origin
   container.Insert(gv_zeros);

// Velocity
   container.Insert(gv_zeros);

// Magnetic field
   double Bmag = 0.311 / unit_magnetic_fluid;
   GeoVector B0(0.0, 0.0, Bmag);
   container.Insert(B0);

// Largest absolute step size
   double RE = 6.37e8 / unit_length_fluid;
   double dmax_fraction = 0.1;
   double dmax0 = dmax_fraction * RE;
   container.Insert(dmax0);

// Reference distance
   container.Insert(RE);

// Relative step size
   container.Insert(dmax_fraction);

   background.SetupObject(container);

   std::ofstream outfile_x;
   std::ofstream outfile_y;
   std::ofstream outfile_z;
   std::string background_file = "main_test_dipole_visualization_B";
   outfile_x.open(background_file + "x.dat");
   outfile_y.open(background_file + "y.dat");
   outfile_z.open(background_file + "z.dat");
// Output field
   int N = 1000;
   int M = 100;
   double corner = 3.0 * RE;
   double dx = 2.0 * corner / (double)(N - 1);
   double dz = 2.0 * corner / (double)(N - 1);
   pos[1] = 0.0;
   outfile_x << std::setprecision(8); 
   outfile_y << std::setprecision(8); 
   outfile_z << std::setprecision(8); 
   for(i = 0; i < N; i++) {
      pos[0] = -corner + i * dx;
      for(j = 0; j < N; j++) {
         pos[2] = -corner + j * dz;
         background.GetFields(t, pos, mom, spdata);
         outfile_x << std::setw(16) << spdata.Bvec[0] * unit_magnetic_fluid;
         outfile_y << std::setw(16) << spdata.Bvec[1] * unit_magnetic_fluid;
         outfile_z << std::setw(16) << spdata.Bvec[2] * unit_magnetic_fluid;
      };
      outfile_x << std::endl;
      outfile_y << std::endl;
      outfile_z << std::endl;
      if(i % M == 0) std::cerr << "i = " << i << std::endl;
   };
   std::cerr << "i = " << i << std::endl;
   outfile_x.close();
   outfile_y.close();
   outfile_z.close();

   std::cout << std::endl;
   std::cout << "DIPOLE FIELD VISUALIZATION (y=0)" << std::endl;
   std::cout << "++++++++++++++++++++" << std::endl;
   std::cout << "Dipole axis: " << UnitVec(B0) << std::endl;
   std::cout << "Resolution: " << N << " x " << N << std::endl;
   std::cout << "Domain size (km): " << 2.0*corner*unit_length_fluid/1.0e5 << " x " << 2.0*corner*unit_length_fluid/1.0e5 << std::endl;
   std::cout << "Domain center: " << gv_zeros << std::endl;
   std::cout << "++++++++++++++++++++" << std::endl;
   std::cout << "Background outputed to " << background_file << "*.dat" << std::endl;
   std::cout << std::endl;

   return 0;
};


