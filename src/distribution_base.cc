/*!
\file distribution_base.cc
\brief Implements a base class to record distributions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "distribution_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
DistributionBase<HConfig>::DistributionBase(void)
                : Params("", STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
DistributionBase<HConfig>::DistributionBase(const std::string& name_in, uint16_t status_in)
                : Params(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <typename HConfig>
DistributionBase<HConfig>::DistributionBase(const DistributionBase<HConfig>& other)
                : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param [in] construct Whether called from a copy constructor or separately
*/
template <typename HConfig>
void DistributionBase<HConfig>::SetupDistribution(bool construct)
{
};

/*!
\author Vladimir Florinski
\date 06/18/2021
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupDistribution()" method.
*/
template <typename HConfig>
void DistributionBase<HConfig>::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupDistribution(false);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 09/13/2022
*/
template <typename HConfig>
void DistributionBase<HConfig>::AddEvent(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
*/
template <typename HConfig>
void DistributionBase<HConfig>::AddRecord(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] n_records_in Total number of events
*/
template <typename HConfig>
void DistributionBase<HConfig>::SetNRecords(int n_records_in)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
*/
template <typename HConfig>
void DistributionBase<HConfig>::ResetDistribution(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
*/
template <typename HConfig>
void DistributionBase<HConfig>::ResetRecords(void)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] other Second distribution
\return Reference to this object
*/
template <typename HConfig>
DistributionBase<HConfig>& DistributionBase<HConfig>::operator +=(const DistributionBase<HConfig>& other)
{
   return *this;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] other Second distribution
*/
template <typename HConfig>
void DistributionBase<HConfig>::CopyRecords(const DistributionBase<HConfig>& other)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[out] size Size of an element of the "distro" array
\return Pointer to the distribution storage
*/
template <typename HConfig>
void* DistributionBase<HConfig>::GetDistroAddress(size_t& size)
{
   size = 0;
   return nullptr;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[out] size Size of an element of the "weights_record" array
\return Pointer to the weights storage
*/
template <typename HConfig>
void* DistributionBase<HConfig>::GetWeightsRecordAddress(size_t& size)
{
   size = 0;
   return nullptr;
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/05/2025
\param[in] coords1        First coordinates data class
\param[in] fields1   First fields data class
\param[in] edata1     First extrema data class
\param[in] coords2        Second coordinates data class
\param[in] fields2   Second fields data class
\param[in] edata2     Second extrema data class
\param[in] action_in Action to calculate the weight
*/
template <typename HConfig>
void DistributionBase<HConfig>::ProcessTrajectory(const Coordinates& coords1, const Fields& fields1, const MagExtrema& edata1,
                                         const Coordinates& coords2, const Fields& fields2, const MagExtrema& edata2,
                                         int action_in)
{
   _coords1 = coords1;
   _fields1 = fields1;
   _edata1 = edata1;
   _coords2 = coords2;
   _fields2 = fields2;
   _edata2 = edata2;
   EvaluateValue();
   EvaluateWeight(action_in);
   AddEvent();
   if (keep_records) AddRecord();
};

/*!
\author Vladimir Florinski
\date 09/13/2022
*/
template <typename HConfig>
void DistributionBase<HConfig>::EvaluateValue(void)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] action_in Action index from "ActionTable"
*/
template <typename HConfig>
void DistributionBase<HConfig>::EvaluateWeight(int action_in)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] file_name Distribution file name
*/
template <typename HConfig>
void DistributionBase<HConfig>::Dump(const std::string& file_name) const
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] file_name Distribution file name
*/
template <typename HConfig>
void DistributionBase<HConfig>::Restore(const std::string& file_name)
{
};
   
/*!
\author Vladimir Florinski
\date 09/13/2019
\param[in] ijk        Which component to print (the remaining dimensions are collapsed)
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
template <typename HConfig>
void DistributionBase<HConfig>::Print1D(int ijk, const std::string& dist_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/29/2022
\param[in] ijk1       Component 1 to print
\param[in] ijk2       Component 2 to print
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
template <typename HConfig>
void DistributionBase<HConfig>::Print2D(int ijk1, int ijk2, const std::string& dist_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
template <typename HConfig>
void DistributionBase<HConfig>::PrintRecords(const std::string& dist_name, bool phys_units) const
{
};

};
