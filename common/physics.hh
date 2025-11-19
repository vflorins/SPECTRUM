/*!
\file physics.hh
\brief Declares physical constants and conversion routines for various physical quantities
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_PHYSICS_HH
#define SPECTRUM_PHYSICS_HH

#include "common/specie.hh"
#include "common/coordinates.hh"

namespace Spectrum {

namespace Fluid {

//-TODO: temporary code start (should be moved out)---------------------------------------------------------------------------------------------------

//! Length: primary (1 au is a good value for the heliosphere)
SPECTRUM_CONSTEXPR double unit_length = SPC_CONST_CGSM_ASTRONOMICAL_UNIT;

//! Velocity: primary (100 km/s is a good value for the heliosphere)
SPECTRUM_CONSTEXPR double unit_velocity = 1.0E7;

//! Mass unit: primary (m_p, m_e, etc.)
SPECTRUM_CONSTEXPR double unit_mass = SPC_CONST_CGSM_MASS_PROTON;

//! Number density: primary (should be 1 cm^-3). We can choose a number density unit independently of the length unit because certain conservation laws (namely MHD) are invariant under the transformation rho->a*rho, p->a*p, B^2->a*B^2. However, this is not the case for other physical laws such as those governing collisions, where corrections must are applied.
SPECTRUM_CONSTEXPR double unit_number_density = 1.0;

//! Temperature: primary (should be 1 eV/k_B~10^4 K)
SPECTRUM_CONSTEXPR double unit_temperature = SPC_CONST_CGSM_ELECTRON_VOLT / SPC_CONST_CGSM_BOLTZMANN;

//-temporary code end---------------------------------------------------------------------------------------------------------------------------------

//! Time: derived from length and velocity (should be hours/days)
SPECTRUM_CONSTEXPR double unit_time = unit_length / unit_velocity;

//! Density: derived from mass and number density (should be 1 m_p*cm^-3)
SPECTRUM_CONSTEXPR double unit_density = unit_mass * unit_number_density;

//! Magnetic field: derived from density and velocity (should be uG)
SPECTRUM_CONSTEXPR double unit_magnetic = unit_velocity * sqrt(unit_density);

//! Pressure and energy density: derived from magnetic field (should be 10^-10)
SPECTRUM_CONSTEXPR double unit_pressure = Sqr(unit_magnetic);

//! Speed of light in fluid units
SPECTRUM_CONSTEXPR double c_code = SPC_CONST_CGSM_SPEED_OF_LIGHT / unit_velocity;

//! Boltzmann constant. Note that for fluids the ideal gas EOS is p=n*k_B*T.
SPECTRUM_CONSTEXPR double kb_code = SPC_CONST_CGSM_BOLTZMANN / unit_pressure * unit_temperature * unit_number_density;

//! Ratio of charge/mass has a conversion factor required to calculate particle's cyclotron frequency and radius.
//SPECTRUM_CONSTEXPR double charge_mass_conversion ((unit_charge_particle * unit_magnetic_fluid / unit_mass_particle / unit_velocity_fluid) / unit_frequency_fluid)

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Fluid-fluid routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

namespace Fluid {

/*!
\brief Calculate gas pressure from temperature
\author Vladimir Florinski
\date 08/08/2025
\param[in] nden Number density in fluid units
\param[in] T    Temperature in fluid units
\return Sound speed in fluid units
*/
SPECTRUM_DEVICE_FUNC inline constexpr double Pressure(double nden, double T)
{
   return nden * kb_code * T;
};

/*!
\brief Calculate the MHD energy density
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density in fluid units
\param[in] u2  Square of velocity in fluid units
\param[in] prs Gas pressure in fluid units
\param[in] B2  Square of the magnetic field (if defined) in fluid units
\return Energy per unit volume in fluid units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double EnergyDensity(double den, double u2, double prs, double B2 = 0.0)
{
   return 0.5 * den * u2 + prs / (specie.pt_idx - 1.0) + B2 / M_8PI;
};

/*!
\brief Calculate gas pressure from MHD energy density
\author Vladimir Florinski
\date 07/31/2019
\param[in] den  Mass density
\param[in] u2   Square of velocity
\param[in] eden Energy density
\param[in] B2   Square of the magnetic field (if defined)
\return gas pressure
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double Pressure(double den, double u2, double eden, double B2 = 0.0)
{
   return (eden - 0.5 * den * u2 - B2 / M_8PI) * (specie.pt_idx - 1.0);
};

/*!
\brief Compute sound speed
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density in fluid units
\param[in] prs Pressure in fluid units
\return Sound speed in fluid units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double SoundSpeed(double den, double prs)
{
   return sqrt(specie.pt_idx * prs / den);
};

/*!
\brief Compute square of the sound speed
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density in fluid units
\param[in] prs 
\param[in] prs 
\return Sound speed squared in fluid units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double Sound2(double den, double prs)
{
   return specie.pt_idx * prs / den;
};

/*!
\brief Compute Alfven speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Density in fluid units
\param[in] B2  Square of the magnetic field in fluid units
\return Alfven speed in fluid units
*/
SPECTRUM_DEVICE_FUNC inline constexpr double AlfvenSpeed(double den, double B2)
{
   return sqrt(B2 / den / M_4PI);
};

/*!
\brief Compute square of the Alfven speed
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density in fluid units
\param[in] B2  Square of the magnetic field in fluid units
\return Alfven speed squared in fluid units
*/
SPECTRUM_DEVICE_FUNC inline constexpr double Alfven2(double den, double B2)
{
   return B2 / den / M_4PI;
};

/*!
\brief Compute fast magnetosonic speed
\author Vladimir Florinski
\date 11/10/2025
\param[in] Cs2   Square of the sound speed in fluid units
\param[in] Va2   Square of the Alfven speed in fluid units
\param[in] cosa2 Square of the cosine of the propagation angle
\return Fast magnetosonic speed in fluid units
*/
SPECTRUM_DEVICE_FUNC inline constexpr double FastMagnetosonicSpeed(double Cs2, double Va2, double cosa2)
{
   return sqrt(Cs2 + Va2 + sqrt(Sqr(Cs2 + Va2) - 4.0 * Cs2 * Va2 * cosa2)) / M_SQRT2;
};

/*!
\brief Compute slow magnetosonic speed
\author Vladimir Florinski
\date 11/10/2025
\param[in] Cs2   Square of the sound speed in fluid units
\param[in] Va2   Square of the Alfven speed in fluid units
\param[in] cosa2 Square of the cosine of the propagation angle
\return Slow magnetosonic speed in fluid units
*/
SPECTRUM_DEVICE_FUNC inline constexpr double SlowMagnetosonicSpeed(double Cs2, double Va2, double cosa2)
{
   return sqrt(Cs2 + Va2 - sqrt(Sqr(Cs2 + Va2) - 4.0 * Cs2 * Va2 * cosa2)) / M_SQRT2;
};

/*!
\brief Compute squares of both fast and slow magnetosonic speeds
\author Vladimir Florinski
\date 11/10/2025
\param[in]  Cs2   Square of the sound speed in fluid units
\param[in]  Va2   Square of the Alfven speed in fluid units
\param[in]  cosa2 Square of the cosine of the propagation angle
\param[out] Vf2   Square of the fast magnetosonic speed in fluid units
\param[out] Vs2   Square of the slow magnetosonic speed in fluid units
*/
SPECTRUM_DEVICE_FUNC inline constexpr void Magnetosonic2(double Cs2, double Va2, double cosa2, double& Vf2, double& Vs2)
{
   double Vd2 = sqrt(Sqr(Cs2 + Va2) - 4.0 * Cs2 * Va2 * cosa2);
   Vf2 = 0.5 * (Cs2 + Va2 + Vd2);
   Vs2 = Vf2 - Vd2;
};

/*!
\brief Effective temperature (Maxwell) based on characteristic velocity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/11/2025
\param[in] v_th Thermal speed
\return Effective temperature
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double EffectiveTemperature(double v_th)
{
   return specie.mass * Sqr(v_th) / (2.0 * kb_code);
};

};

namespace Particle {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle-particle routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Relativistic factor
\author Vladimir Florinski
\date 07/15/2020
\param[in] vmag Velocity magnitude in particle units
\return Lorentz factor
*/
SPECTRUM_DEVICE_FUNC inline constexpr double RelFactor(double vmag)
{
   return 1.0 / sqrt(1.0 - Sqr(vmag / c_code));
};

/*!
\brief Relativistic factor based on momentum
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mmag Momentum magnitude in particle units
\return Lorentz factor
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double RelFactor1(double mmag)
{
   return sqrt(1.0 + Sqr(mmag / specie.mass / c_code));
};

/*!
\brief Calculate particle momentum from its kinetic energy
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] enr Kinetic energy in particle units
\return Relativistic momentum in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double Mom(double enr)
{
   return sqrt(enr * (enr + 2.0 * specie.mass * Sqr(c_code))) / c_code;
};

/*!
\brief Calculate particle total energy from momentum magnitude
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum magnitude in particle units
\return Total relativistic energy in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double EnrTot(double mmag)
{
   return sqrt(Sqr(mmag) + Sqr(specie.mass * c_code)) * c_code;
};

/*!
\brief Calculate particle kinetic energy from momentum magnitude
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum magnitude in particle units
\return Relativistic kinetic energy in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double EnrKin(double mmag)
{
   return EnrTot<specie>(mmag) - specie.mass * Sqr(c_code);
};

/*!
\brief Calculate particle velocity magnitude from momentum magnitude
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum magnitude in particle units
\return Velocity in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double Vel(double mmag)
{
   return mmag * Sqr(c_code) / EnrTot<specie>(mmag);
};

/*!
\brief Calculate particle velocity from momentum
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/10/2025
\param[in] mom Momentum in particle units
\return Velocity in particle units
*/
template <Specie specie, CoordinateSystem vel_sys>
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector Vel(const GeoVector& mom)
{
   double mom_mag = Metric<vel_sys>::Length(mom);
   GeoVector retval = mom;
   Metric<vel_sys>::Rescale(retval, Vel<specie>(mom_mag) / mom_mag);
   return retval;
};

/*!
\brief Calculate particle momentum from velocity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/10/2025
\param[in] vel Velocity in particle units
\return Relativistic momentum in particle units
*/
template <Specie specie, CoordinateSystem vel_sys>
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector Mom(const GeoVector& vel)
{
   GeoVector retval = vel;
   Metric<vel_sys>::Rescale(retval, specie.mass * RelFactor(Metric<vel_sys>::Length(retval)));
   return retval;
};

/*!
\brief Calculate particle rigidity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mmag Momentum magnitude in particle units
\return Rigidity (\f$pc/q\f$) in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double Rigidity(double mmag)
{
   return mmag * c_code / fabs(specie.charge);
};

/*!
\brief Effective temperature (Maxwell) based on characteristic velocity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/11/2025
\param[in] v_th Thermal speed
\return Effective temperature
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double EffectiveTemperature(double v_th)
{
   return specie.mass * Sqr(v_th) / (2.0 * kb_code);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle-fluid routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Electric field from ideal Ohm's law
\author Vladimir Florinski
\date 11/17/2025
\param[in] vel  Velocity vector in particle units
\param[in] Bvec Magnetic field vector in particle units
\return Electric field, \f$\mathbf{E}=-(\mathbf{v}\times\mathbf{B})/c\f$
*/
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector InducedEfield(const GeoVector& vel, const GeoVector& Bvec)
{
   return  -(vel ^ Bvec) / c_code;
};

/*!
\brief Magnetic force on a particle
\author Vladimir Florinski
\date 11/11/2025
\param[in] vel  Velocity in particle units
\param[in] Evec Electric field in particle units
\param[in] Bvec Magnetic field in particle units
\return Magnetic Lorentz force in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr GeoVector LorentzForce(const GeoVector& vel, const GeoVector& Evec, const GeoVector& Bvec)
{
   return specie.charge * (Evec + (vel ^ Bvec) / c_code);
};

/*!
\brief Thermal speed (Maxwell) based on temperature
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/1/2025
\param[in] T Temperature in particle units
\return Thermal speed
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double ThermalSpeed(double T)
{
   return sqrt(2.0 * kb_code * T / specie.mass);
};

/*!
\brief Cyclotron frequency
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/11/2025
\param[in] vmag Velocity magnitude in particle units
\param[in] Bmag Magnetic field magnitude in particle units
\return Cyclotron frequency in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double CyclotronFrequency(double vmag, double Bmag)
{
   return specie.charge * Bmag / (RelFactor(vmag) * specie.mass * c_code);
};

/*!
\brief Larmor radius
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/11/2025
\param[in] mmag Momentum magnitude in particle units
\param[in] Bmag Magnetic field magnitude in particle units
\return Larmor radius
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double LarmorRadius(double mmag, double Bmag)
{
   return mmag * c_code / (fabs(specie.charge) * Bmag);
};

/*!
\brief Magnetic moment
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/11/2025
\param[in] mperp Momentum (perp. component) in particle units
\param[in] Bmag  Magnetic field magnitude in particle units
\return Magnetic moment
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double MagneticMoment(double mperp, double Bmag)
{
   return Sqr(mperp) / (2.0 * specie.mass * Bmag);
};

/*!
\brief Perpendicular momentum from magnetic moment
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mag_mom Magnetic moment in particle units
\param[in] Bmag    Magnetic field magnitude in particle units
\return Momentum (perp. component) in particle units
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double PerpMomentum(double mag_mom, double Bmag)
{
   return sqrt(2.0 * specie.mass * mag_mom * Bmag);
};

/*!
\brief Relativistic factor based on parallel momentum and magnetic moment
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mpara   Parallel component of momentum in particle units
\param[in] mag_mom Magnetic moment in particle units
\param[in] Bmag    Magnetic field magnitude in particle units
\return Lorentz factor
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline constexpr double RelFactor2(double mpara, double mag_mom, double Bmag)
{
   return sqrt(1.0 + (2.0 * specie.mass * mag_mom * Bmag + Sqr(mpara)) / Sqr(specie.mass * c_code));
};

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Plasma physics routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma frequency
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density
\return Plasma frequency
*/
// TODO fix this function
//template <Specie specie>
//SPECTRUM_DEVICE_FUNC inline double PlasmaFrequency(double den)
//{
//   return sqrt(M_4PI * den) * charge_mass_conversion * fabs(specie.charge) / specie.mass;
//};

/*!
\brief Collision frequency
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density
\param[in] T   Temperature
\param[in] Lam Coulomb logarithm
\return Collision frequency
\note Our system of units (MHD centric) requires a coefficient to relate the number density unit to the inverse cube of distance unit.
*/
// TODO fix this function
//template <Specie specie>
//SPECTRUM_DEVICE_FUNC inline double CollisionFrequency(double den, double T, double Lam)
//{
   //return 2.0 / 3.0 * M_SQRT2 / M_SQRTPI * Lam * specie.mass / Cube(ThermalSpeed<specie>(T) * unit_length_fluid) / unit_density_fluid
          //* unit_mass_particle * Sqr(PlasmaFrequency<specie>(den) * charge_mass_conversion * specie.charge / specie.mass);
//};

//! Print all units and test each conversion routine
void TestPhysics(void);

};

#endif
