/*!
\file physics.cc
\brief Defines some supplementary unit conversion routines
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <format>

#include "common/print_warn.hh"
#include "common/physics.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 11/10/2025
*/
void TestPhysics(void)
{
   std::cout << std::endl;

   std::cout << "-------------------------------Physical constants-------------------------------\n";
   std::cout << "CGS speed of light: " << SPC_CONST_CGSM_SPEED_OF_LIGHT << "⨯cm⨯s⁻¹\n";
   std::cout << "CGS mass of electron: " << SPC_CONST_CGSM_MASS_ELECTRON << "⨯g\n";
   std::cout << "CGS mass of proton: " << SPC_CONST_CGSM_MASS_PROTON << "⨯g\n";
   std::cout << "CGS electron charge: " << SPC_CONST_CGSM_ELECTRON_CHARGE << "⨯Fr\n";
   std::cout << "CGS Boltzmann constant: " << SPC_CONST_CGSM_BOLTZMANN << "⨯erg⨯K⁻¹\n";
   std::cout << "CGS astronomical unit: " << SPC_CONST_CGSM_ASTRONOMICAL_UNIT << "⨯cm\n";
   std::cout << "CGS electron volt: " << SPC_CONST_CGSM_ELECTRON_VOLT << "⨯erg\n";
   std::cout << std::endl;

   std::cout << "---------------------------------Particle-units---------------------------------\n";
   std::cout << "Unit of length: " << Particle::unit_length << "⨯cm\n";
   std::cout << "Unit of velocity: " << Particle::unit_velocity << "⨯cm⨯s⁻¹\n";
   std::cout << "Unit of mass: " << Particle::unit_mass << "⨯g\n";
   std::cout << "Unit of charge: " << Particle::unit_charge << "⨯Fr\n";
   std::cout << "Unit of temperature: " << Particle::unit_temperature << "⨯K\n";
   std::cout << "Unit of momentum: " << Particle::unit_momentum << "⨯g⨯cm⨯s⁻¹\n";
   std::cout << "Unit of energy: " << Particle::unit_energy << "⨯erg or " << Particle::unit_energy / SPC_CONST_CGSM_MEGA_ELECTRON_VOLT << "⨯MeV\n";
   std::cout << "Unit of rigidity: " << Particle::unit_rigidity << " erg*Fr^-1\n";
   std::cout << "Unit of time: " << Particle::unit_time << "⨯s\n";
   std::cout << "Unit of force: " << Particle::unit_force << "⨯g⨯cm⨯s⁻²\n";
   std::cout << "Unit of magnetic field: " << Particle::unit_magnetic << "⨯G\n";
   std::cout << "Unit of diffusion: " << Particle::unit_diffusion << "⨯cm²⨯s⁻¹\n";
   std::cout << "Unit of magnetic moment: " << Particle::unit_magnetic_moment << "⨯Fr⨯cm\n";
//   std::cout << "EM to particle conversion unit: " << Particle::charge_conv << "\n";
   std::cout << std::endl;

   std::cout << "----------------------------------Fluid-units-----------------------------------\n";

   std::cout << "Unit of length: " << Fluid::unit_length << "⨯cm or " << Fluid::unit_length / SPC_CONST_CGSM_ASTRONOMICAL_UNIT << "⨯au\n";
   std::cout << "Unit of velocity: " << Fluid::unit_velocity << "⨯cm⨯s⁻¹\n";
   std::cout << "Unit of mass: " << Fluid::unit_mass << "⨯g\n";
   std::cout << "Unit of time: " << Fluid::unit_time << "⨯s\n";
   std::cout << "Unit of number density: " << Fluid::unit_number_density << "⨯cm^⁻³\n";
   std::cout << "Unit of temperature: " << Fluid::unit_temperature << "⨯K\n";
   std::cout << "Unit of density: " << Fluid::unit_density << "⨯g⨯cm^⁻³\n";
   std::cout << "Unit of pressure: " << Fluid::unit_pressure << "⨯dyn⨯cm⁻²\n";
   std::cout << "Unit of magnetic field: " << Fluid::unit_magnetic << "⨯G\n";
   std::cout << std::endl;

   std::cout << "-----------------------------Testing fluid routines-----------------------------\n";

   std::string message;
   double res1, res2;

   Specie<SpecieId::proton_core> specie;
   double m_cgs = specie.mass * Particle::unit_mass;
   double n_cgs = 0.14;
   double rho_cgs = n_cgs * SPC_CONST_CGSM_MASS_PROTON;
   double T_cgs = 7.0E4;
   double P_cgs = n_cgs * SPC_CONST_CGSM_BOLTZMANN * T_cgs;
   double u_cgs = 1.0E7;
   GeoVector B_cgs = {1.5E-6, -2.5E-6, -4.2E-6};
   double Bmag_cgs = B_cgs.Norm();
   double e_cgs = 10.5 * P_cgs;
   double Cs_cgs = 8.0E6;
   double Va_cgs = 4.6E6;
   double alpha = DegToRad(40.0);
   double vth_cgs = 1.3E7;

   message = "Pressure at " + std::format("n={:.3e}⨯cm⁻³, ",  n_cgs)
                            + std::format("T={:.3e}⨯K", T_cgs);
   res1 = Fluid::Pressure(n_cgs / Fluid::unit_number_density, T_cgs / Fluid::unit_temperature);
   res2 = P_cgs / Fluid::unit_pressure;
   PrintPassFail(message, res1, res2);

   message = "Energy density at " + std::format("ρ={:.3e}⨯g⨯cm⁻³, ", rho_cgs)
                                  + std::format("u={:.3e}⨯cm⨯s⁻¹, ", u_cgs)
                                  + std::format("P={:.3e}⨯dyn⨯cm⁻², ", P_cgs)
                                  + std::format("B={:.3e}⨯G", Bmag_cgs);
   
   res1 = Fluid::EnergyDensity<specie>(rho_cgs / Fluid::unit_density,  Sqr(u_cgs / Fluid::unit_velocity),
                                       P_cgs / Fluid::unit_pressure, Sqr(Bmag_cgs / Fluid::unit_magnetic));
   res2 = (0.5 * rho_cgs * Sqr(u_cgs) + P_cgs / (specie.pt_idx - 1.0) + Sqr(Bmag_cgs) / M_8PI) / Fluid::unit_pressure;
   PrintPassFail(message, res1, res2);

   message = "Pressure at " + std::format("ρ={:.3e}⨯g⨯cm⁻³, ", rho_cgs)
                            + std::format("u={:.3e}⨯cm⨯s⁻¹, ", u_cgs)
                            + std::format("ε={:.3e}⨯dyn⨯cm⁻², ", e_cgs)
                            + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Fluid::Pressure<specie>(rho_cgs / Fluid::unit_density,  Sqr(u_cgs / Fluid::unit_velocity),
                                  e_cgs / Fluid::unit_pressure, Sqr(Bmag_cgs / Fluid::unit_magnetic));
   res2 = (e_cgs - 0.5 * rho_cgs * Sqr(u_cgs) - Sqr(Bmag_cgs) / M_8PI) * (specie.pt_idx - 1.0) / Fluid::unit_pressure;
   PrintPassFail(message, res1, res2);


   message = "Sound speed at " + std::format("ρ={:.3e}⨯g⨯cm⁻³, ", rho_cgs)
                               + std::format("P={:.3e}⨯dyn⨯cm⁻², ", P_cgs);
   res1 = Fluid::SoundSpeed<specie>(rho_cgs / Fluid::unit_density, P_cgs / Fluid::unit_pressure);
   res2 = sqrt(specie.pt_idx * P_cgs / rho_cgs) / Fluid::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Alfven speed at " + std::format("ρ={:.3e}⨯g⨯cm⁻³, ", rho_cgs)
                                + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Fluid::AlfvenSpeed(rho_cgs / Fluid::unit_density, Sqr(Bmag_cgs) / Fluid::unit_pressure);
   res2 = Bmag_cgs / sqrt(M_4PI * rho_cgs) / Fluid::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Fast magnetosonic speed at " + std::format("Cs={:.3e}⨯cm⨯s⁻¹, ", Cs_cgs)
                                           + std::format("Va={:.3e}⨯cm⨯s⁻¹, ", Va_cgs)
                                           + std::format("alpha={:.3e}", alpha);
   res1 = Fluid::FastMagnetosonicSpeed(Sqr(Cs_cgs / Fluid::unit_velocity), Sqr(Va_cgs / Fluid::unit_velocity), Sqr(cos(alpha)));
   res2 = sqrt(Sqr(Cs_cgs) + Sqr(Va_cgs) + sqrt(Sqr(Sqr(Cs_cgs) + Sqr(Va_cgs)) - Sqr(2.0 * Cs_cgs * Va_cgs * cos(alpha))))
          / M_SQRT2 / Fluid::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Slow magnetosonic speed at " + std::format("Cs={:.3e}⨯cm⨯s⁻¹, ", Cs_cgs)
                                           + std::format("Va={:.3e}⨯cm⨯s⁻¹, ", Va_cgs)
                                           + std::format("alpha={:.3e}", alpha);
   res1 = Fluid::SlowMagnetosonicSpeed(Sqr(Cs_cgs / Fluid::unit_velocity), Sqr(Va_cgs / Fluid::unit_velocity), Sqr(cos(alpha)));
   res2 = sqrt(Sqr(Cs_cgs) + Sqr(Va_cgs) - sqrt(Sqr(Sqr(Cs_cgs) + Sqr(Va_cgs)) - Sqr(2.0 * Cs_cgs * Va_cgs * cos(alpha))))
          / M_SQRT2 / Fluid::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Temperature at " + std::format("v_th={:.3e}⨯cm⨯s⁻¹", vth_cgs);
   res1 = Fluid::EffectiveTemperature<specie>(vth_cgs / Fluid::unit_velocity);
   res2 = m_cgs * Sqr(vth_cgs) / (2.0 * SPC_CONST_CGSM_BOLTZMANN) / Fluid::unit_temperature;
   PrintPassFail(message, res1, res2);

   std::cout << std::endl;

   std::cout << "----------------------Testing particle/particle routines------------------------\n";

   const CoordinateSystem vel_sys = CoordinateSystem::Cartesian;
   double q_cgs = specie.charge * Particle::unit_charge;
   GeoVector v_cgs = GeoVector(0.23, 0.05, 0.56) * SPC_CONST_CGSM_SPEED_OF_LIGHT;
   double vmag_cgs = v_cgs.Norm();
   GeoVector p_cgs = GeoVector(0.27, 0.41, -0.13) * m_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT;
   double pmag_cgs = p_cgs.Norm();
   double K_cgs = 0.23 * m_cgs * Sqr(SPC_CONST_CGSM_SPEED_OF_LIGHT);
   GeoVector E_cgs = 1.0E-5 * B_cgs;
   double Emag_cgs = E_cgs.Norm();
   double mag_mom_cgs = q_cgs * 2.7E11;
   double ppara_cgs = (p_cgs * B_cgs) / B_cgs.Norm();

   message = "Relativistic factor at " + std::format("v={:.3e}⨯cm⨯s⁻¹", vmag_cgs);
   res1 = Particle::RelFactor(vmag_cgs / Particle::unit_velocity);
   double rel = 1.0 / sqrt(1.0 - Sqr(vmag_cgs / SPC_CONST_CGSM_SPEED_OF_LIGHT));
   res2 = rel;
   PrintPassFail(message, res1, res2);

   message = "Relativistic factor at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹", pmag_cgs);
   res1 = Particle::RelFactor1<specie>(pmag_cgs / Particle::unit_momentum);
   res2 = sqrt(1.0 + Sqr(pmag_cgs / (m_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT)));
   PrintPassFail(message, res1, res2);

   message = "Momentum at " + std::format("K={:.3e}⨯erg", K_cgs);
   res1 = Particle::Mom<specie>(K_cgs / Particle::unit_energy);
   res2 = sqrt(K_cgs * (K_cgs / Sqr(SPC_CONST_CGSM_SPEED_OF_LIGHT) + 2.0 * m_cgs)) / Particle::unit_momentum;
   PrintPassFail(message, res1, res2);

   message = "Total energy at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹", pmag_cgs);
   res1 = Particle::EnrTot<specie>(pmag_cgs / Particle::unit_momentum);
   double Et_cgs = sqrt(Sqr(pmag_cgs) + Sqr(m_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT)) * SPC_CONST_CGSM_SPEED_OF_LIGHT;
   res2 = Et_cgs / Particle::unit_energy;
   PrintPassFail(message, res1, res2);

   message = "Kinetic energy at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹", pmag_cgs);
   res1 = Particle::EnrKin<specie>(pmag_cgs / Particle::unit_momentum);
   res2 = (Et_cgs - m_cgs * Sqr(SPC_CONST_CGSM_SPEED_OF_LIGHT)) / Particle::unit_energy;
   PrintPassFail(message, res1, res2);

   message = "Velocity at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹", pmag_cgs);
   res1 = Particle::Vel<specie>(pmag_cgs / Particle::unit_momentum);
   res2 = pmag_cgs * Sqr(SPC_CONST_CGSM_SPEED_OF_LIGHT) / Et_cgs / Particle::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Velocity at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹", pmag_cgs);
   res1 = Particle::Vel<specie, vel_sys>(p_cgs / Particle::unit_momentum).Norm();
   res2 = p_cgs.Norm() * Sqr(SPC_CONST_CGSM_SPEED_OF_LIGHT) / Et_cgs / Particle::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Momentum at " + std::format("v={:.3e}⨯cm⨯s⁻¹", vmag_cgs);
   res1 = Particle::Mom<specie, vel_sys>(v_cgs / Particle::unit_velocity).Norm();
   res2 = m_cgs * v_cgs.Norm() / sqrt(1.0 - Sqr(v_cgs.Norm() / SPC_CONST_CGSM_SPEED_OF_LIGHT)) / Particle::unit_momentum;
   PrintPassFail(message, res1, res2);

   message = "Rigidity at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹", pmag_cgs);
   res1 = Particle::Rigidity<specie>(pmag_cgs / Particle::unit_momentum);
   res2 = pmag_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT / q_cgs / Particle::unit_rigidity;
   PrintPassFail(message, res1, res2);

   message = "Temperature at " + std::format("v_th={:.3e}⨯cm⨯s⁻¹", vth_cgs);
   res1 = Particle::EffectiveTemperature<specie>(vth_cgs / Particle::unit_velocity);
   res2 = m_cgs * Sqr(vth_cgs) / (2.0 * SPC_CONST_CGSM_BOLTZMANN) / Particle::unit_temperature;
   PrintPassFail(message, res1, res2);

   std::cout << std::endl;

   std::cout << "------------------------Testing particle/fluid routines-------------------------\n";

   message = "Induced Electric field at " + std::format("v={:.3e}⨯cm⨯s⁻¹, ", vmag_cgs)
                                          + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Particle::InducedEfield(v_cgs / Particle::unit_velocity, B_cgs / Particle::unit_magnetic).Norm();
   res2 = (v_cgs ^ B_cgs).Norm() / SPC_CONST_CGSM_SPEED_OF_LIGHT / Particle::unit_magnetic;
   PrintPassFail(message, res1, res2);

   message = "Lorentz force at " + std::format("v={:.3e}⨯cm⨯s⁻¹, ", vmag_cgs)
                                 + std::format("E={:.3e}⨯G, ", Emag_cgs)
                                 + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Particle::LorentzForce<specie>(v_cgs / Particle::unit_velocity, E_cgs / Particle::unit_magnetic, B_cgs / Particle::unit_magnetic).Norm();
   res2 = q_cgs * (E_cgs + (v_cgs ^ B_cgs) / SPC_CONST_CGSM_SPEED_OF_LIGHT).Norm() / Particle::unit_force;
   PrintPassFail(message, res1, res2);

   message = "Thermal speed at " + std::format("T={:.3e}⨯K", T_cgs);
   res1 = Particle::ThermalSpeed<specie>(T_cgs / Particle::unit_temperature);
   res2 = sqrt(2.0 * SPC_CONST_CGSM_BOLTZMANN * T_cgs / m_cgs) / Particle::unit_velocity;
   PrintPassFail(message, res1, res2);

   message = "Cyclotron frequency at " + std::format("v={:.3e}⨯cm⨯s⁻¹, ", vmag_cgs)
                                       + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Particle::CyclotronFrequency<specie>(vmag_cgs / Particle::unit_velocity, Bmag_cgs / Particle::unit_magnetic);
   res2 = q_cgs * Bmag_cgs / (rel * m_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT) * Particle::unit_time;
   PrintPassFail(message, res1, res2);

   message = "Larmor radius at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹, ", pmag_cgs)
                                 + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Particle::LarmorRadius<specie>(pmag_cgs / Particle::unit_momentum, Bmag_cgs / Particle::unit_magnetic);
   res2 = pmag_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT / q_cgs / Bmag_cgs / Particle::unit_length;
   PrintPassFail(message, res1, res2);

   message = "Magnetic moment at " + std::format("p={:.3e}⨯g⨯cm⨯s⁻¹, ", pmag_cgs)
                                   + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Particle::MagneticMoment<specie>(pmag_cgs / Particle::unit_momentum, Bmag_cgs / Particle::unit_magnetic);
   res2 = Sqr(pmag_cgs) / (2.0 * m_cgs * Bmag_cgs) / Particle::unit_magnetic_moment;
   PrintPassFail(message, res1, res2);

   message = "Perpendicular momentum at " + std::format("M={:.3e}⨯Fr⨯cm, ", mag_mom_cgs)
                                          + std::format("B={:.3e}⨯G", Bmag_cgs);
   res1 = Particle::PerpMomentum<specie>(mag_mom_cgs / Particle::unit_magnetic_moment, Bmag_cgs / Particle::unit_magnetic);
   res2 = sqrt(2.0 * m_cgs * mag_mom_cgs * Bmag_cgs) / Particle::unit_momentum;
   PrintPassFail(message, res1, res2);

   message = "Relativistic factor at " + std::format("p∥={:.3e}⨯g⨯cm⨯s⁻¹, ", ppara_cgs)
                                       + std::format("M={:.3e}⨯Fr⨯cm, ", mag_mom_cgs)
                                       + std::format("B={:.3e}⨯G", Bmag_cgs);
   
   res1 = Particle::RelFactor2<specie>(ppara_cgs / Particle::unit_momentum, mag_mom_cgs / Particle::unit_magnetic_moment,
                                                                            Bmag_cgs / Particle::unit_magnetic);
   res2 = sqrt(1.0 + (2.0 * m_cgs * mag_mom_cgs * Bmag_cgs + Sqr(ppara_cgs)) / Sqr(m_cgs * SPC_CONST_CGSM_SPEED_OF_LIGHT));
   PrintPassFail(message, res1, res2);

   std::cout << std::endl;
};

};
