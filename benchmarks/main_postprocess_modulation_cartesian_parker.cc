#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

const double J0 = 1.0;
const double mu = 1.8;
const double T0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT;
const double kap0 = 1.5e22;
const double r0 = GSL_CONST_CGSM_ASTRONOMICAL_UNIT;
const double a = 1.0;
const double b = 2.0;
const double delta = 2.0 / (1.0 - b);
const double eta = 1.0 - delta;
const double alpha = 1.5;
const double V = 4.0e7;

// Unmodulated spectrum
inline double unmod_spectrum(double T) {return J0 * pow(T / T0, -mu);};

// Diffusion coefficient
inline double diff_coeff(double T, double r) {return kap0 * pow(T / T0, a) * pow(r / r0, b);};

// Modulated spectrum
inline double mod_spectrum(double T)
{
   double x = (2.0 * (1.0 + alpha * (mu - 1.0) / 3.0) / Sqr(1.0 - b)) * (V * r0 / diff_coeff(T,r0));
   return (2.0 * unmod_spectrum(T) / std::tgamma(eta)) * pow(x, eta / 2.0) * std::cyl_bessel_k(eta, 2.0 * sqrt(x));
};

int main(int argc, char** argv)
{
   std::string infilename1 = "main_test_modulation_cartesian_TrajectoryParker_spectrum.dat";
   std::string outfilename = "main_test_modulation_cartesian_TrajectoryParker_spectrum_pp.dat";
   std::string line;
   int i, N = 100, sum_c1[N], sum_c2[N];
   double energy1[N], distro1[N], sum_w1[N]; 
   double scal, v, p2, energy0 = SPC_CONST_CGSM_MEGA_ELECTRON_VOLT;

// Open input analytic distro file
   std::ifstream input_spectrum_file1(infilename1);

// Read first two lines of distro file
   std::getline(input_spectrum_file1, line);
   std::getline(input_spectrum_file1, line);

// Read data
   for(i = 0; i < N; i++) {
      input_spectrum_file1 >> energy1[i];
      input_spectrum_file1 >> distro1[i];
      input_spectrum_file1 >> sum_w1[i];
      input_spectrum_file1 >> sum_c1[i];
   };

// Close input cartesian distro file
   input_spectrum_file1.close();

// Open output distro file
   std::ofstream output_spectrum_file(outfilename);

// Output data
   output_spectrum_file << std::setprecision(8);
   scal = Vel(Mom(T0 / unit_energy_particle)) * unmod_spectrum(T0);
   for(i = 0; i < N; i++) {
      p2 = Sqr(Mom(energy1[i] / unit_energy_particle));
      v = Vel(Mom(energy1[i] / unit_energy_particle));
      output_spectrum_file << std::setw(20) << energy1[i] / energy0;
      output_spectrum_file << std::setw(20) << v * unmod_spectrum(energy1[i]) / scal;
      output_spectrum_file << std::setw(20) << v * mod_spectrum(energy1[i]) / scal;
      output_spectrum_file << std::setw(20) << p2 * distro1[i] / scal;
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

   std::cout << "Post-processed distribution files outputed to " << outfilename << std::endl;
   std::cout << std::endl;

   return 0;
};
