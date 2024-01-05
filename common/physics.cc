/*!
\file physics.cc
\brief Defines some supplementary unit conversion routines
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "physics.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 07/14/2020
*/
void PrintUnits(void)
{
   std::cout << std::endl;
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
   std::cerr << "------------------------------particle------------------------------------------\n";
   std::cout << "Unit of mass: " << unit_mass_particle << " g\n";
   std::cout << "Unit of charge: " << unit_charge_particle << " CGSq\n";
   std::cout << "Unit of energy: " << unit_energy_particle << " erg\n";
   std::cout << "Unit of momentum: " << unit_momentum_particle << " g cm s^-1\n";
   std::cout << "Charge to mass factor: " << charge_mass_particle << "\n";
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cout << "Speed of light in code units: " << c_code << "\n";
   std::cout << "Boltzmann constant in code units: " << kb_code << "\n";
   std::cout << "Particle mass in code units: " << mass[0] << "\n";
   std::cout << "Particle charge in code units: " << charge[0] << "\n";
   std::cout << "Thermal speed in code units: " << ThermalSpeed(1.0, 0) << "\n";
   std::cout << "Alfven speed in code units: " << AlfvenSpeed(1.0, 1.0) << "\n";
   std::cout << "Cyclotron frequency in code units: " << CyclotronFrequency(1.0, 1.0, 0) << "\n";
   std::cout << "Larmor radius in code units: " << LarmorRadius(RelFactor(1.0) * mass[0], 1.0, 0) << "\n";
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << std::endl;
};

};
