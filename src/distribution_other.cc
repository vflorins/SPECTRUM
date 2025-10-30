/*!
\file distribution_other.cc
\brief Implements several event distribution classes
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "distribution_other.hh"
#include "geometry/coordinates.hh"
#include "common/physics.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/19/2024
*/
template <class distroClass>
DistributionUniform<distroClass>::DistributionUniform()
                                : DistributionTemplated<distroClass>("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 06/18/2021
*/
template <class distroClass>
DistributionUniform<distroClass>::DistributionUniform(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                                : DistributionTemplated<distroClass>(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <class distroClass>
DistributionUniform<distroClass>::DistributionUniform(const DistributionUniform& other)
                                : DistributionTemplated<distroClass>(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] construct Whether called from a copy constructor or separately
*/
template <class distroClass>
void DistributionUniform<distroClass>::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<distroClass>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(val_hot);
   this->container.Read(val_cold);

// Place the actions into the table
   this->ActionTable.push_back([this]() {UniformHot();});
   this->ActionTable.push_back([this]() {UniformCold();});
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
template <class distroClass>
void DistributionUniform<distroClass>::UniformHot(void)
{
   this->_weight = val_hot;
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
template <class distroClass>
void DistributionUniform<distroClass>::UniformCold(void)
{
   this->_weight = val_cold;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionTimeUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2021
*/
DistributionTimeUniform::DistributionTimeUniform(void)
                       : DistributionUniform<double>(dist_name_time_uniform, 0, DISTRO_TIME)
{
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionTimeUniform::DistributionTimeUniform(const DistributionTimeUniform& other)
                       : DistributionUniform<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/11/2024
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionTimeUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(val_time);

// Check that ONLY the first dimension is active.
   if (this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/11/2024
*/
void DistributionTimeUniform::EvaluateValue(void)
{
   if (val_time == 0) this->_value[0] = this->_t;
   else this->_value[0] = this->_t2;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
*/
DistributionPositionUniform::DistributionPositionUniform(void)
                           : DistributionUniform<double>(dist_name_position_uniform, 0, DISTRO_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionPositionUniform::DistributionPositionUniform(const DistributionPositionUniform& other)
                           : DistributionUniform<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionPositionUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(val_time);
   this->container.Read(val_coord);

// Check that ALL three dimensions are active.
   if (this->dims != 7) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
*/
void DistributionPositionUniform::EvaluateValue(void)
{
   if (val_time == 0) this->_value = this->_pos;
   else this->_value = this->_pos2;

   if (val_coord == 1) Metric<CoordinateSystem::Spherical>::PosToCurv(this->_value);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionMomentumUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/27/2024
*/
DistributionMomentumUniform::DistributionMomentumUniform(void)
                           : DistributionUniform<double>(dist_name_momentum_uniform, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/27/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionMomentumUniform::DistributionMomentumUniform(const DistributionMomentumUniform& other)
                           : DistributionUniform<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 02/27/2024
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionMomentumUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(val_time);
   this->container.Read(val_coord);

// Check that ALL three dimensions are active.
   if (this->dims != 7) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 02/27/2024

If val_coord == 0, then "native coordinates" are used, meaning that the momentum vector is left in whatever coordinates are used for trajectory integration.
If val_coord == 1, then the momentum vector is converted to locally spherical coordinates with B || z. In this case, only the momentum magnitude and pitch angle (if available) are recorded, so the distribution effectively becomes 2D.
*/
void DistributionMomentumUniform::EvaluateValue(void)
{
   GeoVector momentum, bhat;
   if (val_time == 0) {
      momentum = this->_mom;
      bhat = this->_spdata.bhat;
   }
   else {
      momentum = this->_mom2;
      bhat = this->_spdata2.bhat;
   };

   if (val_coord == 0) this->_value = momentum;
   else {
#if (TRAJ_TYPE == TRAJ_FOCUSED) || (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE)
// Focused and Parker trajectories are already in locally spherical coordinates
   this->_value = momentum;
#elif TRAJ_TYPE == TRAJ_FIELDLINE
   this->_value[0] = momentum[2];
   this->_value[1] = 0.0;
   this->_value[2] = 0.0;
#elif (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   this->_value[0] = momentum.Norm();
   this->_value[1] = momentum[2] / this->_value[0];
   this->_value[2] = 0.0;
#elif TRAJ_TYPE == TRAJ_LORENTZ
   this->_value[0] = momentum.Norm();
   this->_value[1] = momentum * bhat / this->_value[0];
   this->_value[2] = 0.0;
#endif
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionMomentumUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
DistributionPositionMomentumUniform::DistributionPositionMomentumUniform(void)
                                   : DistributionUniform<double>(dist_name_position_momentum_uniform, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionPositionMomentumUniform::DistributionPositionMomentumUniform(const DistributionPositionMomentumUniform& other)
                                   : DistributionUniform<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionPositionMomentumUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionUniform<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(val_time);
   this->container.Read(pos_idx);
   this->container.Read(mom_idx);

// Check that ONLY the first two dimensions are active.
   if (this->dims != 5) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
void DistributionPositionMomentumUniform::EvaluateValue(void)
{
   if (val_time == 0) {
      this->_value[0] = (pos_idx == 3 ? this->_pos.Norm() : this->_pos[pos_idx]);
      this->_value[1] = this->_mom[mom_idx];
   }
   else {
      this->_value[0] = (pos_idx == 3 ? this->_pos2.Norm() : this->_pos2[pos_idx]);
      this->_value[1] = this->_mom2[mom_idx];
   };
};

#if TRAJ_TYPE == TRAJ_LORENTZ

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionAnisotropyLISM methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
DistributionAnisotropyLISM::DistributionAnisotropyLISM(void)
                          : DistributionTemplated<double>(dist_name_anisotropy_LISM, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionAnisotropyLISM::DistributionAnisotropyLISM(const DistributionAnisotropyLISM& other)
                          : DistributionTemplated<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionAnisotropyLISM::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(rot_matrix[0]);
   this->container.Read(rot_matrix[1]);
   this->container.Read(rot_matrix[2]);
   this->container.Read(U_LISM);
   this->container.Read(mom_pow_law);
   this->container.Read(grad_perp_dens);

// Place the actions into the table
   this->ActionTable.push_back([this]() {ComptonGettingFactor();});
   this->ActionTable.push_back([this]() {MomPowerLawAnisotropy();});
   this->ActionTable.push_back([this]() {FirstLegendreAnisotropy();});
   this->ActionTable.push_back([this]() {SecondLegendreAnisotropy();});
   this->ActionTable.push_back([this]() {bCrossGradientAnisotropy();});

// Check that ONLY the first two dimensions are active.
   if (this->dims != 5) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
void DistributionAnisotropyLISM::EvaluateValue(void)
{
// Find incoming direction in specified coordinate frame
   mom_rel = this->_mom;
   mom_rel.ChangeToBasis(rot_matrix);
   Metric<CoordinateSystem::Spherical>::PosToCurv(mom_rel);
   this->_value[0] = mom_rel[1];
   this->_value[1] = mom_rel[2];

// Find relative momentum in LISM. Reuse "mom_rel" in EvaluateWeight().
   mom_rel = this->_mom2 - RelFactor1(this->_mom2.Norm(), this->specie) * mass[this->specie] * U_LISM;
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
void DistributionAnisotropyLISM::ComptonGettingFactor(void)
{
   double vel;
   GeoVector mom_hat;

   vel = Vel(this->_mom.Norm(), this->specie);
   mom_hat = UnitVec(this->_mom);

// The Comptom-Getting factor is an approximation of the momentum power law anisotropy for "U_LISM" << "c_code"
   this->_weight = 1.0 - mom_pow_law * U_LISM * mom_hat / vel;
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
void DistributionAnisotropyLISM::MomPowerLawAnisotropy(void)
{
   this->_weight = pow(mom_rel.Norm() / this->_mom.Norm(), mom_pow_law);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
void DistributionAnisotropyLISM::FirstLegendreAnisotropy(void)
{
   this->_weight = UnitVec(mom_rel) * this->_spdata2.bhat;
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
void DistributionAnisotropyLISM::SecondLegendreAnisotropy(void)
{
   this->_weight = 0.5 * (3.0 * Sqr(UnitVec(mom_rel) * this->_spdata2.bhat) - 1.0);
};

/*!
\author Juan G Alonso Guzman
\date 12/26/2023
*/
void DistributionAnisotropyLISM::bCrossGradientAnisotropy(void)
{
// FIXME: This is according to Zhang et al. 2020, but (perhaps) differs from Zhang et al. 2014. We should investigate this further.
   this->_weight = grad_perp_dens * (this->_pos2 + LarmorRadius(mom_rel.Norm(), this->_spdata2.Bmag, this->specie) * (UnitVec(mom_rel) ^ this->_spdata2.bhat));
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumKineticEnergyPowerLaw
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2021
*/

DistributionSpectrumKineticEnergyPowerLaw::DistributionSpectrumKineticEnergyPowerLaw(void)
                                         : DistributionTemplated<double>(dist_name_spectrum_kinetic_energy_power_law, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/

DistributionSpectrumKineticEnergyPowerLaw::DistributionSpectrumKineticEnergyPowerLaw(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                                         : DistributionTemplated<double>(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/13/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionSpectrumKineticEnergyPowerLaw::DistributionSpectrumKineticEnergyPowerLaw(const DistributionSpectrumKineticEnergyPowerLaw& other)
                                         : DistributionTemplated<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionSpectrumKineticEnergyPowerLaw::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(J0);
   this->container.Read(T0);
   this->container.Read(pow_law);
   this->container.Read(val_cold);

// Place the actions into the table
   this->ActionTable.push_back([this]() {SpectrumKineticEnergyPowerLawHot();});
   this->ActionTable.push_back([this]() {SpectrumKineticEnergyPowerLawCold();});

// Check that ONLY the first dimension is active.
   if (this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionSpectrumKineticEnergyPowerLaw::EvaluateValue(void)
{
#if (TRAJ_TYPE == TRAJ_FOCUSED) || (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE)
   this->_value[0] = EnrKin(this->_mom[0], this->specie);
#elif TRAJ_TYPE == TRAJ_FIELDLINE
   this->_value[0] = EnrKin(this->_mom[2], this->specie);
#elif (TRAJ_TYPE == TRAJ_LORENTZ) || (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   this->_value[0] = EnrKin(this->_mom.Norm(), this->specie);
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
void DistributionSpectrumKineticEnergyPowerLaw::SpectrumKineticEnergyPowerLawHot(void)
{
   double mom2mag;
#if (TRAJ_TYPE == TRAJ_FOCUSED) || (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE)
   mom2mag = this->_mom2[0];
#elif TRAJ_TYPE == TRAJ_FIELDLINE
   mom2mag = this->_mom2[2];
#elif (TRAJ_TYPE == TRAJ_LORENTZ) || (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   mom2mag = this->_mom2.Norm();
#endif
   kin_energy = EnrKin(mom2mag, this->specie);

#if DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE == 0
   double velocity = Vel(mom2mag, this->specie);
// The power law is the differential density U=f(p)*p^2/v, but the weighting function is f(p) itself, so a division by p^2 and multiplication by v is required here.
   this->_weight = J0 * velocity * pow(kin_energy / T0, pow_law) / Sqr(mom2mag);
#elif DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE == 1
// The power law is the differential intensity J=f(p)*p^2, but the weighting function is f(p) itself, so a division by p^2 is required here.
   this->_weight = J0 * pow(kin_energy / T0, pow_law) / Sqr(mom2mag);
#elif DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE == 2
// The power law is the distribution function f(p), no weighting necessary
   this->_weight = J0 * pow(kin_energy / T0, pow_law);
#else
   std::cerr << "DistributionSpectrumKineticEnergyPowerLaw Error: TYPE not recognized." << std::endl;
#endif
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
void DistributionSpectrumKineticEnergyPowerLaw::SpectrumKineticEnergyPowerLawCold(void)
{
   this->_weight = val_cold;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumKineticEnergyBentPowerLaw
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
*/

DistributionSpectrumKineticEnergyBentPowerLaw::DistributionSpectrumKineticEnergyBentPowerLaw(void)
                                             : DistributionSpectrumKineticEnergyPowerLaw(dist_name_spectrum_kinetic_energy_bent_power_law, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionSpectrumKineticEnergyBentPowerLaw::DistributionSpectrumKineticEnergyBentPowerLaw(const DistributionSpectrumKineticEnergyBentPowerLaw& other)
                                             : DistributionSpectrumKineticEnergyPowerLaw(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionSpectrumKineticEnergyBentPowerLaw::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionSpectrumKineticEnergyPowerLaw::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   double pow_law_b;
   this->container.Read(T_b);
   this->container.Read(pow_law_b);
   this->container.Read(bend_smoothness);
   pow_law_comb = (pow_law - pow_law_b) / bend_smoothness;

// Check that ONLY the first dimension is active.
   if (this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
*/
void DistributionSpectrumKineticEnergyBentPowerLaw::SpectrumKineticEnergyPowerLawHot(void)
{
   DistributionSpectrumKineticEnergyPowerLaw::SpectrumKineticEnergyPowerLawHot();
// The power law is the differential intensity J=f(p)*p^2, but the weighting function is f(p) itself, so a division by p^2 is required here.
   this->_weight /= pow(1.0 + pow(kin_energy / T_b, pow_law_comb), bend_smoothness);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroPositionCumulativeO1
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/

DistributionPositionCumulativeOrder1::DistributionPositionCumulativeOrder1(void)
                                    : DistributionTemplated<GeoVector>(dist_name_pos_cumulative_O1, 0, DISTRO_TIME)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionPositionCumulativeOrder1::DistributionPositionCumulativeOrder1(const DistributionPositionCumulativeOrder1& other)
                                    : DistributionTemplated<GeoVector>(other)
{
   RAISE_BITS(this->_status, DISTRO_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionPositionCumulativeOrder1::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<GeoVector>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

// Place the actions into the table
   this->ActionTable.push_back([this]() {RecordPosition();});

// Check that ONLY the first dimension is active.
   if (this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void DistributionPositionCumulativeOrder1::EvaluateValue(void)
{
   this->_value[0] = this->_t2;
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void DistributionPositionCumulativeOrder1::RecordPosition(void)
{
   this->_weight = this->_pos2 - this->_pos;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroPositionCumulativeO2
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/

DistributionPositionCumulativeOrder2::DistributionPositionCumulativeOrder2(void)
                                    : DistributionTemplated<GeoMatrix>(dist_name_pos_cumulative_O2, 0, DISTRO_TIME)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionPositionCumulativeOrder2::DistributionPositionCumulativeOrder2(const DistributionPositionCumulativeOrder2& other)
                                    : DistributionTemplated<GeoMatrix>(other)
{
   RAISE_BITS(this->_status, DISTRO_TIME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionPositionCumulativeOrder2::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<GeoMatrix>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

// Place the actions into the table
   this->ActionTable.push_back([this]() {RecordPosition();});

// Check that ONLY the first dimension is active.
   if (this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void DistributionPositionCumulativeOrder2::EvaluateValue(void)
{
   this->_value[0] = this->_t2;
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void DistributionPositionCumulativeOrder2::RecordPosition(void)
{
   this->_weight.Dyadic(this->_pos2 - this->_pos);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroLossCone
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
*/

DistributionLossCone::DistributionLossCone(void)
                    : DistributionTemplated<GeoVector>(dist_name_loss_cone, 0, DISTRO_SPACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionLossCone::DistributionLossCone(const DistributionLossCone& other)
                    : DistributionTemplated<GeoVector>(other)
{
   RAISE_BITS(this->_status, DISTRO_SPACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionLossCone::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<GeoVector>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(val_time);
   this->container.Read(val_coord);

// Place the actions into the table
   this->ActionTable.push_back([this]() {RecordLossCone();});

// Check that ALL three dimensions are active.
   if (this->dims != 7) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
*/
void DistributionLossCone::EvaluateValue(void)
{
   if (val_time == 0) this->_value = this->_pos;
   else this->_value = this->_pos2;

   if (val_coord == 1) Metric<CoordinateSystem::Spherical>::PosToCurv(this->_value);
};

/*!
\author Juan G Alonso Guzman
\date 06/23/2023
*/
void DistributionLossCone::RecordLossCone(void)
{
   if (val_time == 0) this->_weight = GeoVector(_spdata2.Bmag_min, _spdata2.Bmag_max, asin(sqrt(_spdata.Bmag / _spdata2.Bmag_max)));
   else this->_weight = GeoVector(_spdata2.Bmag_min, _spdata2.Bmag_max, asin(sqrt(_spdata2.Bmag / _spdata2.Bmag_max)));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionDivergenceFlowPowerLaw
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 06/18/2021
*/

DistributionDivergenceFlowPowerLaw::DistributionDivergenceFlowPowerLaw(void)
                                 : DistributionTemplated<double>(dist_name_divergence_flow_power_law, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionDivergenceFlowPowerLaw::DistributionDivergenceFlowPowerLaw(const DistributionDivergenceFlowPowerLaw& other)
                                  : DistributionTemplated<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2025
\param[in] construct Whether called from a copy constructor or separately
*/
void DistributionDivergenceFlowPowerLaw::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistributionTemplated<double>::SetupDistribution(false);
   if (BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(A0);
   this->container.Read(pow_law);
   this->container.Read(val_cold);

// Place the actions into the table
   this->ActionTable.push_back([this]() {DivergenceFlowPowerLawHot();});
   this->ActionTable.push_back([this]() {DivergenceFlowPowerLawCold();});

// Check that ONLY the first dimension is active.
   if (this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2025
*/
void DistributionDivergenceFlowPowerLaw::EvaluateValue(void)
{
   this->_value = this->_mom;
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2025
*/
void DistributionDivergenceFlowPowerLaw::DivergenceFlowPowerLawHot(void)
{
   this->_weight = A0 * pow(fabs(_spdata2.divU()), pow_law);
};

/*!
\author Juan G Alonso Guzman
\date 07/01/2025
*/
void DistributionDivergenceFlowPowerLaw::DivergenceFlowPowerLawCold(void)
{
   this->_weight = val_cold;
};

};
