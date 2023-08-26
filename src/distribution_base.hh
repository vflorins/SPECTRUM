/*!
\file distribution_base.hh
\brief Declares a base class to record distributions
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DISTRIBUTION_BASE_HH
#define SPECTRUM_DISTRIBUTION_BASE_HH

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "config.h"
#include "common/params.hh"
#include "common/spatial_data.hh"
#include <functional>

#ifndef TRAJ_TYPE
#error Trajectory type is undefined!
#endif

namespace Spectrum {

//! Distribution is in time, allowed to combine with space and momentum
const uint16_t DISTRO_TIME = 0x0010;

//! Distribution is in space, allowed to combine with time and momentum
const uint16_t DISTRO_SPACE = 0x0020;

//! Distribution is in momentum, allowed to combine with time and space
const uint16_t DISTRO_MOMENTUM = 0x0040;

//! Function type to respond to different actions
using WeightAction = std::function<void(void)>;

/*!
\brief A base class describing a generic binned distribution
\author Vladimir Florinski

The class is optimized for convenience, not speed. All distributions are treated as 3D and addressed through a multi-index. Access to 1D and 2D distributions is therefore slower than using one/two integers. This should not be an issue for most applications because a single event requires integrating a trajectory, which is a lengthy procedure. A distribution could mix space and momentum coordinates.
*/
class DistributionBase : public Params {

protected:

//! Number of dimensions (1-3). To provide a uniform interface all distributions are treated as 3D (persistent)
   int dims;

//! Number of bins in the distribution (persistent)
   MultiIndex n_bins = MultiIndex(1, 1, 1);

//! Smallest value of the variable (persistent)
   GeoVector minval;

//! Largest value of the variable (persistent)
   GeoVector maxval;

//! Use linear or logarithmic bins (persistent)
   MultiIndex log_bins;

//! Whether to add outlying events to the end bins (persistent)
   MultiIndex bin_outside;

//! Physical units of the bin variable (persistent)
   GeoVector unit_val;

//! Lower and upper bin limits; could be logarithmic, etc. (persistent)
   GeoVector limits[2];

//! Range of the parameter (persistent)
   GeoVector range;

//! Table of functions to respond to different boundary events (persistent)
   std::vector <WeightAction> ActionTable;

//! Whether or not to keep records of values and weights (persistent)
   bool keep_records;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Counts array - all dimensions rolled into one (transient)
   std::vector <int> counts;

//! Size of the current bin or its logarithm (transient)
   GeoVector bin_size;

//! Total event count (transient)
   int n_events = 0;

//! Spatial data (transient)
   SpatialData _spdata;

//! Second time value (transient)
   double _t2;

//! Second spatial position (transient)
   GeoVector _pos2;

//! Second momentum vector (transient)
   GeoVector _mom2;

//! Second spatial data (transient)
   SpatialData _spdata2;

//! Internal value to be binned (transient)
   GeoVector _value;

//! Record of values (transient)
   std::vector <GeoVector> values_record;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Default constructor (protected, class not designed to be instantiated)
   DistributionBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DistributionBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   DistributionBase(const DistributionBase& other);

//! Set up the distribution accumulator based on "params" (stub)
   virtual void SetupDistribution(bool construct);

//! Determine the value to be binned from a phase space position and other arguments (stub)
   virtual void EvaluateValue(void);

//! Dispatch routine to call the appropriate weight function (stub)
   virtual void EvaluateWeight(int action_in);

//! Add a single processed event (stub)
   virtual void AddEvent(void);

//! Add a record (stub)
   virtual void AddRecord(void);

public:

//! Destructor
   virtual ~DistributionBase() = default;

//! Clone function (stub)
   virtual std::shared_ptr <DistributionBase> Clone(void) const = 0;

//! Set up the object's persistent class data members (generic)
   void SetupObject(const DataContainer& cont_in);

//! Return the number of bins in all dimensions
   MultiIndex NBins(void) const;

//! Return the minimum value for this bin
   double BinLeft(int ijk, int bin) const;

//! Return the "center" value for this bin
   double BinCent(int ijk, int bin) const;

//! Return the maximum value for this bin
   double BinRght(int ijk, int bin) const;

//! Return total number of events
   int NEvents(void) const;

//! Return toral number of records
   int NRecords(void) const;

//! Set the total number of events
   void SetNEvents(int n_events_in);

//! Set the total number of records
   virtual void SetNRecords(int n_records_in);

//! Return the event count in a bin
   int GetEvents(const MultiIndex& bin) const;

//! Return if records are being kept
   bool GetKeepRecords(void);

//! Add another distribution to this (stub)
   virtual DistributionBase& operator +=(const DistributionBase& other);

//! Copy records from other distribution (stub)
   virtual void CopyRecords(const DistributionBase& other);

//! Clear the distribution counts and weights
   virtual void ResetDistribution(void);

//! Clear values and weights records
   virtual void ResetRecords(void);

//! Return the address of "distro" (stub)
   virtual void* GetDistroAddress(size_t& size);

//! Return the address of "weights_record" (stub)
   virtual void* GetWeightsRecordAddress(size_t& size);

//! Return the address of "counts"
   int* GetCountsAddress(void);

//! Return the address of "values_record"
   GeoVector* GetValuesRecordAddress(void);

//! Analyze the trajectory outcome and record it in the distribution (generic)
   void ProcessTrajectory(double t1, const GeoVector& pos1, const GeoVector& mom1, const SpatialData& spdata1,
                          double t2, const GeoVector& pos2, const GeoVector& mom2, const SpatialData& spdata2,
                          int action_in);

//! Dump the complete distribution to a file (stub)
   virtual void Dump(const std::string& file_name) const;

//! Restore the distribution from a dump file (stub)
   virtual void Restore(const std::string& file_name);

//! Print the reduced distribution in 1D (stub)
   virtual void Print1D(int ijk, const std::string& file_name, bool phys_units) const;

//! Print the reduced distribution in 2D (stub)
   virtual void Print2D(int ijk1, int ijk2, const std::string& file_name, bool phys_units) const;

//! Print records
   virtual void PrintRecords(const std::string& file_name, bool phys_units) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionBase inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
\return Number of bins in each dimension
*/
inline MultiIndex DistributionBase::NBins(void) const
{
   return n_bins;
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] ijk Which dimension
\param[in] bin Bin number
\return Smallest value for this bin
*/
inline double DistributionBase::BinLeft(int ijk, int bin) const
{
   if((bin < 0) || (bin >= n_bins[ijk])) return 0.0;
   if(log_bins[ijk]) return pow(10.0, limits[0][ijk] + bin * bin_size[ijk]);
   else return limits[0][ijk] + bin * bin_size[ijk];
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] ijk Which dimension
\param[in] bin Bin number
\return Center value for this bin (log-center for logarithmic bins)
*/
inline double DistributionBase::BinCent(int ijk, int bin) const
{
   if((bin < 0) || (bin >= n_bins[ijk])) return 0.0;
   if(log_bins[ijk]) return pow(10.0, limits[0][ijk] + (bin + 0.5) * bin_size[ijk]);
   else return limits[0][ijk] + (bin + 0.5) * bin_size[ijk];
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] ijk Which dimension
\param[in] bin Bin number
\return Largest value for this bin
*/
inline double DistributionBase::BinRght(int ijk, int bin) const
{
   if((bin < 0) || (bin >= n_bins[ijk])) return 0.0;
   if(log_bins[ijk]) return pow(10.0, limits[0][ijk] + (bin + 1) * bin_size[ijk]);
   else return limits[0][ijk] + (bin + 1) * bin_size[ijk];
};

/*!
\author Vladimir Florinski
\date 10/30/2019
\return Total number of events
*/
inline int DistributionBase::NEvents(void) const
{
   return n_events;
};

/*!
\author Juan G Alonso Guzman
\date 12/02/2022
\return Total number of records
*/
inline int DistributionBase::NRecords(void) const
{
   return values_record.size();
};

/*!
\author Vladimir Florinski
\date 06/10/2022
\param[in] n_events_in Total number of events
*/
inline void DistributionBase::SetNEvents(int n_events_in)
{
   n_events = n_events_in;
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] bin Bin multi-index
\return Number of events in the bin
*/
inline int DistributionBase::GetEvents(const MultiIndex& bin) const
{
   return counts[n_bins.LinIdx(bin)];
};

/*!
\author Juan G Alonso Guzman
\date 12/02/2022
\return keep_records variable
*/
inline bool DistributionBase::GetKeepRecords(void)
{
   return keep_records;
};

/*!
\author Vladimir Florinski
\date 09/18/2020
\return Pointer to the counts storage
*/
inline int* DistributionBase::GetCountsAddress(void)
{
   return counts.data();
};

/*!
\author Juan G Alonso Guzman
\date 12/02/2022
\return Pointer to the counts storage
*/
inline GeoVector* DistributionBase::GetValuesRecordAddress(void)
{
   return values_record.data();
};

};

#endif
