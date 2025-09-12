#include "common/physics.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

const int specie = SPECIES_PROTON_BEAM;
const double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
const double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;
const double one_MeV = SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle;
const double t_final = 100.0 * one_day;
const double z_spectrum = 0.1 * one_au;
const double p0 = Mom(1.0 * one_MeV, specie);
const double pf = Mom(100.0 * one_MeV, specie);
const double s = 4.0;
const double U_up = 4.0e7 / unit_velocity_fluid;
const double U_dn = U_up / s;
const double DeltaU = U_up - U_dn;
const double kappa_up = 1.0e20 / unit_diffusion_fluid;
const double kappa_dn = kappa_up * Sqr(U_dn / U_up);
const double beta = 1.5 * (s + 1.0) / (s - 1.0);
const double tau = 4.0 * kappa_up / Sqr(U_up);
const double amp = 3.0 / (M_8PI * DeltaU * Cube(p0));

// Phase-space density downstream
inline double f_dn(double z, double p, double t)
{
   double a = pow(p0 / p, beta);
   double b = 2.0 * z / (tau * U_dn);
   double c = sqrt(t / tau);
   double d = 0.5 * beta / c * log(p / p0) + z / (sqrt(t * tau) * U_dn);
   return amp * sqrt(Cube(p0 / p)) * exp(U_dn * z / (2.0 * kappa_dn)) * (exp(-b) * erfc(d - c) * a + exp(b) * erfc(d + c) / a);
};

int main(int argc, char** argv)
{
   int j; 
   std::string infilename;
   std::string outfilename;
   std::string line;
   int Np = 100, sum_c[Np];
   double coord[Np], distro[Np], sum_w[Np];

// Spectrum vs momentum
   infilename = "main_test_diff_shock_acc_TrajectoryParkerSource_spectrum.dat";
   outfilename = "main_test_diff_shock_acc_TrajectoryParkerSource_spectrum_pp.dat";

// Open input analytic distro file
   std::ifstream input_dsa_file(infilename);

// Read first two lines of distro file
   std::getline(input_dsa_file, line);
   std::getline(input_dsa_file, line);

// Read data
   for(j = 0; j < Np; j++) {
      input_dsa_file >> coord[j];
      input_dsa_file >> distro[j];
      input_dsa_file >> sum_w[j];
      input_dsa_file >> sum_c[j];
   };

// Close input cartesian distro file
   input_dsa_file.close();

// Open output distro file
   std::ofstream output_dsa_file(outfilename);

// Output data
   output_dsa_file << std::setprecision(8);
   for(j = 0; j < Np; j++) {
      output_dsa_file << std::setw(20) << EnrKin(coord[j], specie) / one_MeV
                      << std::setw(20) << f_dn(z_spectrum, coord[j], t_final) * M_4PI * Sqr(coord[j])
                      << std::setw(20) << amp * M_8PI * distro[j] * Sqr(coord[j])
                      << std::endl;
   };

// Close output distro file
   output_dsa_file.close();

   std::cout << "Post-processed distribution files outputed to " << outfilename << std::endl;
   std::cout << std::endl;
   return 0;
};
