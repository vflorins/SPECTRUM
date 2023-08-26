/*!
\file distribution_other.cc
\brief Implements several event distribution classes
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "distribution_other.hh"
#include "common/physics.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
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
   if(!construct) DistributionTemplated<distroClass>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(&val_hot);
   this->container.Read(&val_cold);

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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionTimeUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionUniform<double>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   if(this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionTimeUniform::EvaluateValue(void)
{
   this->_value[0] = this->_t2;
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
*/
void DistributionPositionUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionUniform<double>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(&val_time);
   this->container.Read(&val_coord);
};

/*!
\author Juan G Alonso Guzman
\date 06/13/2022
*/
void DistributionPositionUniform::EvaluateValue(void)
{
   if(val_time == 0) this->_value = this->_pos;
   else this->_value = this->_pos2;

   if(val_coord == 1) this->_value.XYZ_RTP();
};

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPitchUniform methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2021
*/
DistributionPitchUniform::DistributionPitchUniform(void)
                        : DistributionUniform<double>(dist_name_pitch_uniform, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionPitchUniform::DistributionPitchUniform(const DistributionPitchUniform& other)
                        : DistributionUniform<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionPitchUniform::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionUniform<double>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(&val_time);

   this->unit_val[0] = 1.0;
   if(this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionPitchUniform::EvaluateValue(void)
{
#if TRAJ_TYPE == TRAJ_FOCUSED
   if(val_time == 0) this->_value[0] = this->_mom[1];
   else this->_value[0] = this->_mom2[1];
#elif (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   if(val_time == 0) this->_value[0] = this->_mom[2] / this->_mom.Norm();
   else this->_value[0] = this->_mom2[2] / this->_mom2.Norm();
#elif TRAJ_TYPE == TRAJ_LORENTZ
   if(val_time == 0) this->_value[0] = (this->_mom * this->_spdata.bhat) / this->_mom.Norm();
   else this->_value[0] = (this->_mom2 * this->_spdata2.bhat) / this->_mom2.Norm();
#endif
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistroSpectrumLISM
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/18/2021
*/

DistributionSpectrumLISM::DistributionSpectrumLISM(void)
                        : DistributionTemplated<double>(dist_name_spectrum_lism, 0, DISTRO_MOMENTUM)
{
};

/*!
\author Vladimir Florinski
\date 05/13/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionSpectrumLISM::DistributionSpectrumLISM(const DistributionSpectrumLISM& other)
                        : DistributionTemplated<double>(other)
{
   RAISE_BITS(this->_status, DISTRO_MOMENTUM);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionSpectrumLISM::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionTemplated<double>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(&J0);
   this->container.Read(&T0);
   this->container.Read(&pow_law);
   this->container.Read(&val_cold);

// Place the actions into the table
   this->ActionTable.push_back([this]() {SpectrumLISMHot();});
   this->ActionTable.push_back([this]() {SpectrumLISMCold();});

   if(this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\date 05/05/2022
*/
void DistributionSpectrumLISM::EvaluateValue(void)
{
   this->_value[0] = EnrKin(this->_mom.Norm(), this->specie);
};

/*!
\author Vladimir Florinski
\date 05/17/2022
*/
void DistributionSpectrumLISM::SpectrumLISMHot(void)
{
   double kin_energy = EnrKin(this->_mom2.Norm(), this->specie);

// The distribution is the intensity J=f(p)*p^2, but the weighting function is f(p) itself, so a division by p^2 is required here.
   this->_weight = J0 * pow(kin_energy / T0, pow_law) / this->_mom2.Norm2();
};


/*!
\author Vladimir Florinski
\date 05/17/2022
*/
void DistributionSpectrumLISM::SpectrumLISMCold(void)
{
   this->_weight = val_cold;
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void DistributionPositionCumulativeOrder1::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionTemplated<GeoVector>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

// Place the actions into the table
   this->ActionTable.push_back([this]() {RecordPosition();});

   if(this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
void DistributionPositionCumulativeOrder2::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionTemplated<GeoMatrix>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

// Place the actions into the table
   this->ActionTable.push_back([this]() {RecordPosition();});

   if(this->dims != 1) LOWER_BITS(this->_status, STATE_SETUP_COMPLETE);
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
\author Juan G Alonso
\date 06/23/2023
*/

DistributionLossCone::DistributionLossCone(void)
                    : DistributionTemplated<GeoVector>(dist_name_loss_cone, 0, DISTRO_SPACE)
{
};

/*!
\author Juan G Alonso
\date 06/23/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionLossCone::DistributionLossCone(const DistributionLossCone& other)
                    : DistributionTemplated<GeoVector>(other)
{
   RAISE_BITS(this->_status, DISTRO_SPACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso
\date 06/23/2023
*/
void DistributionLossCone::SetupDistribution(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistributionTemplated<GeoVector>::SetupDistribution(false);
   if(BITS_LOWERED(this->_status, STATE_SETUP_COMPLETE)) return;

   this->container.Read(&val_time);
   this->container.Read(&val_coord);

// Place the actions into the table
   this->ActionTable.push_back([this]() {RecordLossCone();});
};

/*!
\author Juan G Alonso
\date 06/23/2023
*/
void DistributionLossCone::EvaluateValue(void)
{
   if(val_time == 0) this->_value = this->_pos;
   else this->_value = this->_pos2;

   if(val_coord == 1) this->_value.XYZ_RTP();
};

/*!
\author Juan G Alonso
\date 06/23/2023
*/
void DistributionLossCone::RecordLossCone(void)
{
   if(val_time == 0) this->_weight = GeoVector(_spdata2.Bmag_min, _spdata2.Bmag_max, asin(sqrt(_spdata.Bmag / _spdata2.Bmag_max)));
   else this->_weight = GeoVector(_spdata2.Bmag_min, _spdata2.Bmag_max, asin(sqrt(_spdata2.Bmag / _spdata2.Bmag_max)));
};

};
