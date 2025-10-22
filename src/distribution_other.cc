/*!
\file distribution_other.cc
\brief Implements several event distribution classes
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "distribution_other.hh"
#include "common/physics.hh"

namespace Spectrum {


//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/19/2024
*/
template <typename Trajectory, class distroClass>
DistributionUniform<Trajectory, distroClass>::DistributionUniform()
                                : DistributionTemplated("", STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 06/18/2021
*/
template <typename Trajectory, class distroClass>
DistributionUniform<Trajectory, distroClass>::DistributionUniform(const std::string& name_in, uint16_t status_in)
                                : DistributionTemplated(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename Trajectory, class distroClass>
DistributionUniform<Trajectory, distroClass>::DistributionUniform(const DistributionUniform& other)
                                : DistributionTemplated(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename Trajectory, class distroClass>
void DistributionUniform<Trajectory, distroClass>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(val_hot);
   container.Read(val_cold);

// Place the actions into the table
   ActionTable.push_back([this]() {UniformHot();});
   ActionTable.push_back([this]() {UniformCold();});
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
template <typename Trajectory, class distroClass>
void DistributionUniform<Trajectory, distroClass>::UniformHot(void)
{
   _weight = val_hot;
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
template <typename Trajectory, class distroClass>
void DistributionUniform<Trajectory, distroClass>::UniformCold(void)
{
   _weight = val_cold;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionTimeUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2021
*/
template <typename HConfig>
DistributionTimeUniform<HConfig>::DistributionTimeUniform(void)
                       : DistributionUniform(dist_name_time_uniform, DISTRO_TIME)
{
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionTimeUniform<HConfig>::DistributionTimeUniform(const DistributionTimeUniform& other)
                       : DistributionUniform(other)
{
   RAISE_BITS(_status, DISTRO_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/11/2024
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionTimeUniform<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(val_time);

// Check that ONLY the first dimension is active.
   if (dims != 1) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/11/2024
*/
template <typename HConfig>
void DistributionTimeUniform<HConfig>::EvaluateValue(void)
{
   if (val_time == 0) _value[0] = _coords1.Time();
   else _value[0] = _coords2.Time();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
*/
template <typename HConfig>
DistributionPositionUniform<HConfig>::DistributionPositionUniform(void)
                           : DistributionUniform(dist_name_position_uniform, DISTRO_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionPositionUniform<HConfig>::DistributionPositionUniform(const DistributionPositionUniform& other)
                           : DistributionUniform(other)
{
   RAISE_BITS(_status, DISTRO_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionPositionUniform<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(val_time);
   container.Read(val_coord);

// Check that ALL three dimensions are active.
   if (dims != 7) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
*/
template <typename HConfig>
void DistributionPositionUniform<HConfig>::EvaluateValue(void)
{
   if (val_time == 0) _value = _coords1.Pos();
   else _value = _coords2.Pos();

   if (val_coord == 1) _value.XYZ_RTP();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionMomentumUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/27/2024
*/
template <typename HConfig>
DistributionMomentumUniform<HConfig>::DistributionMomentumUniform(void)
                           : DistributionUniform(dist_name_momentum_uniform, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/27/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionMomentumUniform<HConfig>::DistributionMomentumUniform(const DistributionMomentumUniform& other)
                           : DistributionUniform(other)
{
   RAISE_BITS(_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 02/27/2024
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionMomentumUniform<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(val_time);
   container.Read(val_coord);

// Check that ALL three dimensions are active.
   if (dims != 7) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
 \author Lucius Schoenbaum
\date 08/09/2025

If val_coord == 0, then "native coordinates" are used, meaning that the momentum vector is left in whatever coordinates are used for trajectory integration.
If val_coord == 1, then the momentum vector is converted to locally spherical coordinates with B || z. In this case, only the momentum magnitude and pitch angle (if available) are recorded, so the distribution effectively becomes 2D.
*/
template <typename HConfig>
void DistributionMomentumUniform<HConfig>::EvaluateValue(void)
{
   GeoVector momentum, bhat;
   if (val_time == 0) {
      momentum = _coords1.Mom();
      bhat = _fields1.HatMag();
   }
   else {
      momentum = _coords2.Mom();
      bhat = _fields2.HatMag();
   };

   if (val_coord == 0) _value = momentum;
   else {
      if constexpr (std::same_as<Trajectory, TrajectoryFocused<HConfig>> || std::same_as<Trajectory, TrajectoryParker<HConfig>>) {
// Focused and Parker trajectories are already in locally spherical coordinates
         _value = momentum;
      }
      else if constexpr (std::derived_from<Trajectory, TrajectoryFieldlineBase<HConfig>>) {
         _value[0] = momentum[2];
         _value[1] = 0.0;
         _value[2] = 0.0;
      }
      else if constexpr (std::derived_from<Trajectory, TrajectoryGuiding<HConfig>>) {
         _value[0] = momentum.Norm();
         _value[1] = momentum[2] / _value[0];
         _value[2] = 0.0;
      }
      else if constexpr (std::same_as<Trajectory, TrajectoryLorentz<HConfig>>) {
         _value[0] = momentum.Norm();
         _value[1] = momentum * bhat / _value[0];
         _value[2] = 0.0;
      }
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionMomentumUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
template <typename HConfig>
DistributionPositionMomentumUniform<HConfig>::DistributionPositionMomentumUniform(void)
                                   : DistributionUniform(dist_name_position_momentum_uniform, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionPositionMomentumUniform<HConfig>::DistributionPositionMomentumUniform(const DistributionPositionMomentumUniform& other)
                                   : DistributionUniform(other)
{
   RAISE_BITS(_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionPositionMomentumUniform<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(val_time);
   container.Read(pos_idx);
   container.Read(mom_idx);

// Check that ONLY the first two dimensions are active.
   if (dims != 5) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
template <typename HConfig>
void DistributionPositionMomentumUniform<HConfig>::EvaluateValue(void)
{
   if (val_time == 0) {
      _value[0] = _coords1.Pos()[pos_idx];
      _value[1] = _coords1.Mom()[mom_idx];
   }
   else {
      _value[0] = _coords2.Pos()[pos_idx];
      _value[1] = _coords2.Mom()[mom_idx];
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionAnisotropyLISM methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
DistributionAnisotropyLISM<HConfig>::DistributionAnisotropyLISM(void)
                          : DistributionTemplated(dist_name_anisotropy_LISM, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionAnisotropyLISM<HConfig>::DistributionAnisotropyLISM(const DistributionAnisotropyLISM& other)
                          : DistributionTemplated(other)
{
   RAISE_BITS(_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(rot_matrix[0]);
   container.Read(rot_matrix[1]);
   container.Read(rot_matrix[2]);
   container.Read(U_LISM);
   container.Read(mom_pow_law);
   container.Read(grad_perp_dens);

// Place the actions into the table
   ActionTable.push_back([this]() {ComptonGettingFactor();});
   ActionTable.push_back([this]() {MomPowerLawAnisotropy();});
   ActionTable.push_back([this]() {FirstLegendreAnisotropy();});
   ActionTable.push_back([this]() {SecondLegendreAnisotropy();});
   ActionTable.push_back([this]() {bCrossGradientAnisotropy();});

// Check that ONLY the first two dimensions are active.
   if (dims != 5) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::EvaluateValue(void)
{
// Find incoming direction in specified coordinate frame
   mom_rel = _coords1.Mom();
   mom_rel.ChangeToBasis(rot_matrix);
   mom_rel.XYZ_RTP();
   _value[0] = mom_rel[1];
   _value[1] = mom_rel[2];

// Find relative momentum in LISM. Reuse "mom_rel" in EvaluateWeight().
   mom_rel = _coords2.Mom() - RelFactor1<specie>(_coords2.Mom().Norm()) * specie.mass * U_LISM;
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::ComptonGettingFactor(void)
{
   double vel;
   GeoVector mom_hat;

   vel = Vel<specie>(_coords1.Mom().Norm());
   mom_hat = UnitVec(_coords1.Mom());

// The Comptom-Getting factor is an approximation of the momentum power law anisotropy for "U_LISM" << "c_code"
   _weight = 1.0 - mom_pow_law * U_LISM * mom_hat / vel;
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::MomPowerLawAnisotropy(void)
{
   _weight = pow(mom_rel.Norm() / _coords1.Mom().Norm(), mom_pow_law);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::FirstLegendreAnisotropy(void)
{
   _weight = UnitVec(mom_rel) * _fields2.HatMag();
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::SecondLegendreAnisotropy(void)
{
   _weight = 0.5 * (3.0 * Sqr(UnitVec(mom_rel) * _fields2.HagMat()) - 1.0);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
template <typename HConfig>
void DistributionAnisotropyLISM<HConfig>::bCrossGradientAnisotropy(void)
{
// FIXME: This is according to Zhang et al. 2020, but (perhaps) differs from Zhang et al. 2014. We should investigate this further.
   _weight = grad_perp_dens * (_coords2.Pos() + LarmorRadius<specie>(mom_rel.Norm(), _fields2.AbsMag()) * (UnitVec(mom_rel) ^ _fields2.HatMag()));
};

//#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumKineticEnergyPowerLaw
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2021
*/
template <typename HConfig>
DistributionSpectrumKineticEnergyPowerLaw<HConfig>::DistributionSpectrumKineticEnergyPowerLaw(void)
                                         : DistributionTemplated(dist_name_spectrum_kinetic_energy_power_law, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
DistributionSpectrumKineticEnergyPowerLaw<HConfig>::DistributionSpectrumKineticEnergyPowerLaw(const std::string& name_in, uint16_t status_in)
                                         : DistributionTemplated(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/13/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionSpectrumKineticEnergyPowerLaw<HConfig>::DistributionSpectrumKineticEnergyPowerLaw(const DistributionSpectrumKineticEnergyPowerLaw& other)
                                         : DistributionTemplated(other)
{
   RAISE_BITS(_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionSpectrumKineticEnergyPowerLaw<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(J0);
   container.Read(T0);
   container.Read(pow_law);
   container.Read(val_cold);

// Place the actions into the table
   ActionTable.push_back([this]() {SpectrumKineticEnergyPowerLawHot();});
   ActionTable.push_back([this]() {SpectrumKineticEnergyPowerLawCold();});

// Check that ONLY the first dimension is active.
   if (dims != 1) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/09/2025
*/
template <typename HConfig>
void DistributionSpectrumKineticEnergyPowerLaw<HConfig>::EvaluateValue(void)
{
   if constexpr (std::same_as<Trajectory, TrajectoryFocused<HConfig>> || std::same_as<Trajectory, TrajectoryParker<HConfig>>) {
      _value[0] = EnrKin<specie>(_coords1.Mom()[0]);
   }
   else if constexpr (std::derived_from<Trajectory, TrajectoryFieldlineBase<Trajectory, HConfig>>) {
      _value[0] = EnrKin<specie>(_coords1.Mom()[2]);
   }
   else if constexpr (std::same_as<Trajectory, TrajectoryLorentz<HConfig>> || std::derived_from<Trajectory, TrajectoryGuiding<HConfig>>) {
      _value[0] = EnrKin<specie>(_coords1.Mom().Norm());
   }
   else {
// stub
      ;
   }
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/09/2025
*/
template <typename HConfig>
void DistributionSpectrumKineticEnergyPowerLaw<HConfig>::SpectrumKineticEnergyPowerLawHot(void)
{
   double mom2mag;
   if constexpr (std::same_as<Trajectory, TrajectoryFocused<HConfig>> || std::same_as<Trajectory, TrajectoryParker<HConfig>>) {
      mom2mag = _coords2.Mom()[0];
   }
   else if constexpr (std::derived_from<Trajectory, TrajectoryFieldlineBase<HConfig>>) {
      mom2mag = _coords2.Mom()[2];
   }
   else if constexpr (std::same_as<Trajectory, TrajectoryLorentz<HConfig>> || std::derived_from<Trajectory, TrajectoryGuiding<HConfig>>) {
      mom2mag = _coords2.Mom().Norm();
   }
   else {
// stub
      ;
   }
   kin_energy = EnrKin<specie>(mom2mag);

#if DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE == 0
   double velocity = Vel<specie>(mom2mag);
// The power law is the differential density U=f(p)*p^2/v, but the weighting function is f(p) itself, so a division by p^2 and multiplication by v is required here.
   _weight = J0 * velocity * pow(kin_energy / T0, pow_law) / Sqr(mom2mag);
#elif DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE == 1
// The power law is the differential intensity J=f(p)*p^2, but the weighting function is f(p) itself, so a division by p^2 is required here.
   _weight = J0 * pow(kin_energy / T0, pow_law) / Sqr(mom2mag);
#elif DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE == 2
// The power law is the distribution function f(p), no weighting necessary
   _weight = J0 * pow(kin_energy / T0, pow_law);
#else
   std::cerr << "DistributionSpectrumKineticEnergyPowerLaw Error: TYPE not recognized." << std::endl;
#endif
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
template <typename HConfig>
void DistributionSpectrumKineticEnergyPowerLaw<HConfig>::SpectrumKineticEnergyPowerLawCold(void)
{
   _weight = val_cold;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumKineticEnergyBentPowerLaw
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
*/

template <typename HConfig>
DistributionSpectrumKineticEnergyBentPowerLaw<HConfig>::DistributionSpectrumKineticEnergyBentPowerLaw(void)
                                             : DistributionSpectrumKineticEnergyPowerLaw(dist_name_spectrum_kinetic_energy_bent_power_law, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionSpectrumKineticEnergyBentPowerLaw<HConfig>::DistributionSpectrumKineticEnergyBentPowerLaw(const DistributionSpectrumKineticEnergyBentPowerLaw& other)
                                             : DistributionSpectrumKineticEnergyPowerLaw(other)
{
   RAISE_BITS(_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionSpectrumKineticEnergyBentPowerLaw<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionSpectrumKineticEnergyPowerLaw::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   double pow_law_b;
   container.Read(T_b);
   container.Read(pow_law_b);
   container.Read(bend_smoothness);
   pow_law_comb = (pow_law - pow_law_b) / bend_smoothness;

// Check that ONLY the first dimension is active.
   if (dims != 1) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
*/
template <typename HConfig>
void DistributionSpectrumKineticEnergyBentPowerLaw<HConfig>::SpectrumKineticEnergyPowerLawHot(void)
{
   DistributionSpectrumKineticEnergyPowerLaw::SpectrumKineticEnergyPowerLawHot();
// The power law is the differential intensity J=f(p)*p^2, but the weighting function is f(p) itself, so a division by p^2 is required here.
   _weight /= pow(1.0 + pow(kin_energy / T_b, pow_law_comb), bend_smoothness);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroPositionCumulativeO1
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <typename HConfig>
DistributionPositionCumulativeOrder1<HConfig>::DistributionPositionCumulativeOrder1(void)
                                    : DistributionTemplated(dist_name_pos_cumulative_O1, DISTRO_TIME)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionPositionCumulativeOrder1<HConfig>::DistributionPositionCumulativeOrder1(const DistributionPositionCumulativeOrder1& other)
                                    : DistributionTemplated(other)
{
   RAISE_BITS(_status, DISTRO_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionPositionCumulativeOrder1<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

// Place the actions into the table
   ActionTable.push_back([this]() {RecordPosition();});

// Check that ONLY the first dimension is active.
   if (dims != 1) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <typename HConfig>
void DistributionPositionCumulativeOrder1<HConfig>::EvaluateValue(void)
{
   _value[0] = _coords2.Time();
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <typename HConfig>
void DistributionPositionCumulativeOrder1<HConfig>::RecordPosition(void)
{
   _weight = _coords2.Pos() - _coords1.Pos();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroPositionCumulativeO2
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <typename HConfig>
DistributionPositionCumulativeOrder2<HConfig>::DistributionPositionCumulativeOrder2(void)
                                    : DistributionTemplated(dist_name_pos_cumulative_O2, DISTRO_TIME)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionPositionCumulativeOrder2<HConfig>::DistributionPositionCumulativeOrder2(const DistributionPositionCumulativeOrder2& other)
                                    : DistributionTemplated(other)
{
   RAISE_BITS(_status, DISTRO_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionPositionCumulativeOrder2<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

// Place the actions into the table
   ActionTable.push_back([this]() {RecordPosition();});

// Check that ONLY the first dimension is active.
   if (dims != 1) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <typename HConfig>
void DistributionPositionCumulativeOrder2<HConfig>::EvaluateValue(void)
{
   _value[0] = _coords2.Time();
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <typename HConfig>
void DistributionPositionCumulativeOrder2<HConfig>::RecordPosition(void)
{
   _weight.Dyadic(_coords2.Pos() - _coords1.Pos());
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroLossCone
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
*/
template <typename HConfig>
DistributionLossCone<HConfig>::DistributionLossCone(void)
                    : DistributionTemplated(dist_name_loss_cone, DISTRO_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionLossCone<HConfig>::DistributionLossCone(const DistributionLossCone& other)
                    : DistributionTemplated(other)
{
   RAISE_BITS(_status, DISTRO_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
\param[in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionLossCone<HConfig>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated::SetupDistribution(false);
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   container.Read(val_time);
   container.Read(val_coord);

// Place the actions into the table
   ActionTable.push_back([this]() {RecordLossCone();});

// Check that ALL three dimensions are active.
   if (dims != 7) LOWER_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
*/
template <typename HConfig>
void DistributionLossCone<HConfig>::EvaluateValue(void)
{
   if (val_time == 0) _value = _coords1.Pos();
   else _value = _coords2.Pos();

   if (val_coord == 1) _value.XYZ_RTP();
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/09/2025
*/
template <typename HConfig>
void DistributionLossCone<HConfig>::RecordLossCone(void)
{
   if (val_time == 0) _weight = GeoVector(_edata2.Bmag_min, _edata2.Bmag_max, asin(sqrt(_fields1.Mag() / _edata2.Bmag_max)));
   else _weight = GeoVector(_edata2.Bmag_min, _edata2.Bmag_max, asin(sqrt(_fields2.Mag() / _edata2.Bmag_max)));
};

};
