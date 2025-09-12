/*!
\file physics.cc
\brief Defines some supplementary unit conversion routines
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <iostream>
#include <iomanip>

#include "common/physics.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 06/02/2025
*/
void PrintUnits(void)
{

#ifdef USE_GSL

   int spc = SPECIES_PROTON_CORE;

   std::cout << std::endl;
   std::cout << "CGS speed of light: " << GSL_CONST_CGSM_SPEED_OF_LIGHT << std::endl;
   std::cout << "CGS mass of proton: " << GSL_CONST_CGSM_MASS_PROTON << std::endl;
   std::cout << "CGS electron charge: " << SPC_CONST_CGSM_ELECTRON_CHARGE << std::endl;
   std::cout << "CGS Boltzmann constant: " << GSL_CONST_CGSM_BOLTZMANN << std::endl;
   std::cout << "CGS astronomical unit: " << GSL_CONST_CGSM_ASTRONOMICAL_UNIT << std::endl;
   std::cout << "CGS electron volt: " << GSL_CONST_CGSM_ELECTRON_VOLT << std::endl;
   std::cout << "--------------------------------------------------------------------------------\n";
   std::cout << "Printing physical units used by the code\n";
   std::cout << "-------------------------------fluid--------------------------------------------\n";
   std::cout << "Unit of length: " << unit_length_fluid << " cm\n";
   std::cout << "Unit of time: " << unit_time_fluid << " s\n";
   std::cout << "Unit of frequency: " << unit_frequency_fluid << " s^-1\n";
   std::cout << "Unit of velocity: " << unit_velocity_fluid << " cm s^-1\n";
   std::cout << "Unit of number density: " << unit_number_density_fluid << " cm^-3\n";
   std::cout << "Unit of density: " << unit_density_fluid << " g cm^-3\n";
   std::cout << "Unit of pressure: " << unit_pressure_fluid << " dyn cm^-2\n";
   std::cout << "Unit of magnetic field: " << unit_magnetic_fluid << " G\n";
   std::cout << "Unit of diffusion: " << unit_diffusion_fluid << " cm^2 s^-1\n";
   std::cout << "Unit of temperature: " << unit_temperature_fluid << " K\n";
   std::cout << "------------------------------particle------------------------------------------\n";
   std::cout << "Unit of mass: " << unit_mass_particle << " g\n";
   std::cout << "Unit of charge: " << unit_charge_particle << " CGSq\n";
   std::cout << "Unit of energy: " << unit_energy_particle << " erg\n";
   std::cout << "Unit of momentum: " << unit_momentum_particle << " g cm s^-1\n";
   std::cout << "Unit of rigidity: " << unit_rigidity_particle << " erg/CGSq\n";
   std::cout << "Charge to mass factor: " << charge_mass_particle << "\n";

   std::cout << "--------------------------------------------------------------------------------\n";
   std::cout << "Speed of light in code units: " << c_code;
   std::cout << "     Should be: " << GSL_CONST_CGSM_SPEED_OF_LIGHT / unit_velocity_fluid << "\n";
   
   std::cout << "Boltzmann constant in code units: " << kb_code;
   std::cout << "     Should be: " << GSL_CONST_CGSM_BOLTZMANN / unit_energy_particle * unit_temperature_fluid << "\n";

   std::cout << "Particle mass in code units: " << SpeciesMasses[spc];
   std::cout << "     Should be: " << SpeciesMasses[spc] * unit_mass_particle / unit_mass_particle << "\n";

   std::cout << "Particle charge in code units: " << SpeciesCharges[spc];
   std::cout << "     Should be: " << SpeciesCharges[spc] * unit_charge_particle / unit_charge_particle << "\n";

   std::cout << "Thermal speed in code units: " << ThermalSpeed(1.0, spc);
   std::cout << "     Should be: " << sqrt(2.0 * GSL_CONST_CGSM_BOLTZMANN * unit_temperature_fluid / (SpeciesMasses[spc] * unit_mass_particle))
                                    / unit_velocity_fluid << "\n";

   std::cout << "Particle rigidity in code units: " << Rigidity(1.0, spc);
   std::cout << "     Should be: " << unit_momentum_particle * GSL_CONST_CGSM_SPEED_OF_LIGHT / unit_charge_particle / unit_rigidity_particle << "\n";

   std::cout << "Alfven speed in code units: " << AlfvenSpeed(1.0, 1.0);
   std::cout << "     Should be: " << unit_magnetic_fluid / sqrt(M_4PI * unit_density_fluid) / unit_velocity_fluid << "\n";

   std::cout << "Cyclotron frequency in code units: " << CyclotronFrequency(1.0, 1.0, spc);
   std::cout << "     Should be: " << SpeciesCharges[spc] * unit_charge_particle * unit_magnetic_fluid
                                    / (SpeciesMasses[spc] * unit_mass_particle * GSL_CONST_CGSM_SPEED_OF_LIGHT) / unit_frequency_fluid << "\n";
   
   std::cout << "Larmor radius in code units: " << LarmorRadius(SpeciesMasses[spc] * 1.0, 1.0, spc);
   std::cout << "     Should be: " << SpeciesMasses[spc] * unit_mass_particle * unit_velocity_fluid * GSL_CONST_CGSM_SPEED_OF_LIGHT
                                    / (SpeciesCharges[spc] * unit_charge_particle * unit_magnetic_fluid) / unit_length_fluid << "\n";

   std::cout << "Plasma frequency in code units: " << PlasmaFrequency(1.0, spc);
   std::cout << "     Should be: " << sqrt(M_4PI * unit_density_fluid) * SpeciesCharges[spc] * unit_charge_particle
                                    / (SpeciesMasses[spc] * unit_mass_particle) / unit_frequency_fluid << "\n";

   std::cout << "Collision frequency in code units: " << CollisionFrequency(1.0, 1.0, 20.0, spc);
   std::cout << "     Should be: " << 160.0 * M_SQRT2 * M_SQRTPI * unit_density_fluid * Quad(SpeciesCharges[spc] * unit_charge_particle) / 3.0
                                    / Cube(SpeciesMasses[spc] * unit_mass_particle * sqrt(2.0 * GSL_CONST_CGSM_BOLTZMANN * unit_temperature_fluid
                                    / (SpeciesMasses[spc] * unit_mass_particle))) / unit_frequency_fluid << "\n";

   std::cout << "--------------------------------------------------------------------------------\n";
   std::cout << std::endl;

#endif

};

};
