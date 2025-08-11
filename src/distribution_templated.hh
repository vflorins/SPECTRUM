/*!
\file distribution_templated.hh
\brief An extension of the base distribution class for arbitrary type of variables
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DISTRIBUTION_TEMPLATED_HH
#define SPECTRUM_DISTRIBUTION_TEMPLATED_HH

#include "distribution_base.hh"
#include "common/matrix.hh"

namespace Spectrum {

/*!
\brief A base class describing a generic binned distribution
\author Vladimir Florinski

The class is optimized for convenience, not speed. All distributions are treated as 3D and addressed through a multi-index. Access to 1D and 2D distributions is therefore slower than using one/two integers. This should not be an issue for most applications because a single event requires integrating a trajectory, which is a lengthy procedure. A distribution could mix space and momentum coordinates.

Parameters: MultiIndex n_bins, GeoVector minval, GeoVector maxval, MultiIndex log_bins, MultiIndex bin_outside, distroClass unit_distro,
            GeoVector unit_val
*/
template <typename Trajectory_, class distroClass>
class DistributionTemplated : public DistributionBase<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;

   using DistributionBase::container;
   using DistributionBase::unit_val;
   using DistributionBase::_status;
   using DistributionBase::n_bins;
   using DistributionBase::n_events;
   using DistributionBase::minval;
   using DistributionBase::maxval;
   using DistributionBase::log_bins;
   using DistributionBase::bin_outside;
   using DistributionBase::keep_records;
   using DistributionBase::dims;
   using DistributionBase::limits;
   using DistributionBase::bin_size;
   using DistributionBase::range;
   using DistributionBase::counts;
   using DistributionBase::_value;
   using DistributionBase::values_record;
   using DistributionBase::ActionTable;
   using DistributionBase::class_name;
   using DistributionBase::specie;
   using DistributionBase::BinCent;

protected:

//! Physical units of the distro variable (persistent)
   distroClass unit_distro;

//! Distribution array - all dimensions rolled into one (transient)
   std::vector <distroClass> distro;

//! Internal weight to be added to the bin (transient)
   distroClass _weight;

//! Record of weights (transient)
   std::vector <distroClass> weights_record;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Default constructor (protected, class not designed to be instantiated)
   DistributionTemplated(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DistributionTemplated(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   DistributionTemplated(const DistributionTemplated& other);

// Class dependent portion of distribution setup routine
   void ReadUnitDistro(void);

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Dispatch routine to call the appropriate weight function
   void EvaluateWeight(int action_in) override;

//! Add a single processed event
   void AddEvent(void) override;

//! Add a record
   void AddRecord(void) override;

//! Obtain a normalized value in a bin
   distroClass operator [](const MultiIndex& bin) const;

//! Class dependent portion of printing routines
   void PrintSumDistro(std::ofstream& distfile, long sum_counts, distroClass sum_distro, bool phys_units) const;

//! Class dependent portion of record printing routines
   void PrintWeight(std::ofstream& distfile, int record, bool phys_units) const;

public:

//! Add another distribution to this
   DistributionBase& operator +=(const DistributionBase& other) override;

//! Copy records from other distribution
   void CopyRecords(const DistributionBase& other) override;

//! Clear the distribution counts and weights
   void ResetDistribution(void) override;

//! Clear the values and weights records
   void ResetRecords(void) override;

//! Set the total number of records
   void SetNRecords(int n_records_in) override;

//! Return the address of "distro"
   void* GetDistroAddress(size_t& size) override;

//! Return the address of "weights_record"
   void* GetWeightsRecordAddress(size_t& size) override;

//! Dump the complete distribution to a file
   void Dump(const std::string& file_name) const override;

//! Restore the distribution from a dump file
   void Restore(const std::string& file_name) override;

//! Print the reduced distribution in 1D
   void Print1D(int ijk, const std::string& file_name, bool phys_units) const override;

//! Print the reduced distribution in 2D
   void Print2D(int ijk1, int ijk2, const std::string& file_name, bool phys_units) const override;

//! Print records
   void PrintRecords(const std::string& file_name, bool phys_units) const override;
};

};

// Something like this is needed for templated classes
#include "distribution_templated.cc"

#endif
