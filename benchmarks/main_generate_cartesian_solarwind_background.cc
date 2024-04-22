#include "src/background_solarwind.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace Spectrum;

int main(void)
{
   BackgroundSolarWind background;
   SpatialData spdata;
   int Nb = 20;
   std::string Nbs = std::to_string(Nb);
   std::string fname_pattern = "cartesian_backgrounds/parker_" + Nbs + "_" + Nbs + "_" + Nbs;

// Block configuration
   MultiIndex block_size (4, 4, 4); // Number of zones per block in each dimension
   MultiIndex Nblocks (Nb, Nb, Nb); // Number of blocks in each dimension

// Domain limits
   double corner = 100.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   GeoVector dom_min (-corner, -corner, -corner); // Minimum domain values
   GeoVector dom_max ( corner,  corner,  corner); // Maximum domain values

// Variables to output
   int Nvar = 9;
   RAISE_BITS(spdata._mask, BACKGROUND_U);
   RAISE_BITS(spdata._mask, BACKGROUND_B);
   RAISE_BITS(spdata._mask, BACKGROUND_E);

   DataContainer container;
   container.Clear();

// ========== PARKER SPIRAL ==========

// Origin
   container.Insert(gv_zeros);

// Velocity
   double umag = 4.0e7 / unit_velocity_fluid;
   GeoVector u0(umag, 0.0, 0.0);
   container.Insert(u0);

// Magnetic field
   double RS = 6.957e10 / unit_length_fluid;
   double r_ref = 3.0 * RS;
   double BmagE = 5.0e-5 / unit_magnetic_fluid;
   double Bmag_ref = BmagE * Sqr((GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid) / r_ref);
   GeoVector B0(Bmag_ref, 0.0, 0.0);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax_fraction = 0.1;
   double dmax = dmax_fraction * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

// solar rotation vector
   double w0 = twopi / (25.0 * 24.0 * 3600.0) / unit_frequency_fluid;
   GeoVector Omega(0.0, 0.0, w0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction for distances closer to the Sun
   container.Insert(dmax_fraction);

   background.SetupObject(container);

// Generate header file
   int var;

   std::cerr << "Outputting header file. Progress:     ";

   std::string header_fname = fname_pattern + ".info";
   std::ofstream header_file (header_fname);

   header_file << std::setprecision(12);
   header_file << std::setw(20) << Nblocks[0] 
               << std::setw(20) << Nblocks[1]
               << std::setw(20) << Nblocks[2]
               << std::endl;
   header_file << std::setw(20) << dom_min[0]
               << std::setw(20) << dom_min[1]
               << std::setw(20) << dom_min[2]
               << std::endl;
   header_file << std::setw(20) << dom_max[0]
               << std::setw(20) << dom_max[1]
               << std::setw(20) << dom_max[2]
               << std::endl;
   header_file << std::setw(20) << Nvar
               << std::endl;

   header_file.close();

   std::cerr << "\e[4D100%\n";

// Generate data file
   int ib,jb,kb; // block indices
   int iz,jz,kz; // zone indices
   GeoVector block_length = (dom_max - dom_min) / Nblocks; // length of each block
   GeoVector zone_length = (block_length) / block_size; // length of each zone
   GeoVector delta = 0.5 * zone_length; // half-length of each zone
   double t = 0.0; // time at which to output data
   GeoVector center_min, pos; // position vectors
   MultiIndex block_idx, zone_idx; // block and zone indices

   std::cerr << "Outputting data file. Progress:     ";

   std::string data_fname = fname_pattern + ".out";
   std::ofstream data_file (data_fname, std::ofstream::binary);

// Iterate over blocks
   for(kb = 0; kb < Nblocks[2]; kb++) {
      block_idx.k = kb;
      for(jb = 0; jb < Nblocks[1]; jb++) {
         block_idx.j = jb;
         for(ib = 0; ib < Nblocks[0]; ib++) {
            block_idx.i = ib;
// Find position of block min center cell and iterate over zones
            center_min = dom_min + delta + block_idx * block_length;
            for(kz = 0; kz < block_size[2]; kz++) {
               zone_idx.k = kz;
               for(jz = 0; jz < block_size[1]; jz++) {
                  zone_idx.j = jz;
                  for(iz = 0; iz < block_size[0]; iz++) {
                     zone_idx.i = iz;
                     pos = center_min + zone_idx * zone_length;
                     background.GetFields(t, pos, spdata);
                     for(var = 0; var < 3; var++) data_file.write((char*)(&spdata.Uvec[var]), sizeof(double));
                     for(var = 0; var < 3; var++) data_file.write((char*)(&spdata.Bvec[var]), sizeof(double));
                     for(var = 0; var < 3; var++) data_file.write((char*)(&spdata.Evec[var]), sizeof(double));
                  };
               };
            };
         };
      };

      std::cerr << "\e[4D";
      std::cerr << std::setw(3) << int(double(kb) / double(Nblocks[2]) * 100.0) << "%";
   };

   data_file.close();

   std::cerr << "\e[4D100%\n";

   std::cout << std::endl;
   std::cout << "Cartesian solarwind background generated to " << fname_pattern << std::endl;

   return 0;
};


