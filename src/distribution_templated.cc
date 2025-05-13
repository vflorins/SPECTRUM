/*!
\file distribution_base.cc
\brief Implements a base class to record distributions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "distribution_templated.hh"
#include "common/print_warn.hh"
#include <algorithm>
#include <fstream>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <class distroClass>
DistributionTemplated<distroClass>::DistributionTemplated(void)
                                  : DistributionBase()
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
template <class distroClass>
DistributionTemplated<distroClass>::DistributionTemplated(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                                  : DistributionBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
template <class distroClass>
DistributionTemplated<distroClass>::DistributionTemplated(const DistributionTemplated& other)
                                  : DistributionBase(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Juan G Alonso Guzman
\date 09/15/2022
*/
template <class distroClass>
void DistributionTemplated<distroClass>::ReadUnitDistro(void)
{
// This template method should only work for <double> because of &
   container.Read(unit_distro);
};

// Method specialization for GeoVector
template <>
void DistributionTemplated<GeoVector>::ReadUnitDistro(void)
{
   container.Read(unit_distro);
};

// Method specialization for GeoMatrix
template <>
void DistributionTemplated<GeoMatrix>::ReadUnitDistro(void)
{
   container.Read(unit_distro);
};

/*!
\author Vladimir Florinski
\date 08/27/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <class distroClass>
void DistributionTemplated<distroClass>::SetupDistribution(bool construct)
{
   int ijk;

   LOWER_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);

// Unpack the parameters. For 1D and 2D distros the arguments are allowed to have junk for the ignorable dimensions.
   container.Reset();
   container.Read(n_bins);
   container.Read(minval);
   container.Read(maxval);
   container.Read(log_bins);
   container.Read(bin_outside);
   ReadUnitDistro();
   container.Read(unit_val);
   container.Read(keep_records);

// The number of bins in active dimensions cannot be less than one. Set the value to 0 or a negative number for each ignorable dimension. Internally, the code changes that number to 1 to make linear addresing possible. Active dimensions are flagged as bits in "dims".
   dims = 0;
   for (ijk = 0; ijk < 3; ijk++) {
      if (n_bins[ijk] > 0) RAISE_BITS(dims, 1 << ijk);
      n_bins[ijk] = std::max(n_bins[ijk], 1);
   };

// "limits" and "range" keep the transformed min/max values (linear, log, or possibly other functional forms).
   for (ijk = 0; ijk < 3; ijk++) {
      if (log_bins[ijk]) {
         limits[0][ijk] = log10(minval[ijk]);
         limits[1][ijk] = log10(maxval[ijk]);
      }
      else {
         limits[0][ijk] = minval[ijk];
         limits[1][ijk] = maxval[ijk];
      };
   };

   bin_size = (limits[1] - limits[0]) / n_bins;
   range = limits[1] - limits[0];

// Set up the vectors
   distro.resize(n_bins.Prod());
   counts.resize(n_bins.Prod());

   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 08/27/2024
*/
template <class distroClass>
void DistributionTemplated<distroClass>::AddEvent(void)
{
   int ijk, lin_bin;
   MultiIndex bin = mi_zeros;
   double val;

// The event has been already processed and "_values" and "_weight" computed
   for (ijk = 0; ijk < 3; ijk++) {

// Skip if dimension is ignorable
      if (BITS_LOWERED(dims, 1 << ijk)) continue;

// Shift value relative to lower limit
      val = (log_bins[ijk] ? log10(_value[ijk]) : _value[ijk]) - limits[0][ijk];
     
// Check for "left" outlier event
      if (val < 0.0) {
         if (bin_outside[ijk]) bin[ijk] = 0;
         else return;
      };

// Compute bin for shifted value. Bin will not be negative but could be >= "n_bins".
      bin[ijk] = val / bin_size[ijk];

// Check for "right" outlier events
      if (bin[ijk] >= n_bins[ijk]) {
         if (bin_outside[ijk]) bin[ijk] = n_bins[ijk] - 1;
         else return;
      };
   };

// Update the counts
   lin_bin = n_bins.LinIdx(bin);
   distro[lin_bin] += _weight;
   counts[lin_bin]++;
   n_events++;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
*/
template <class distroClass>
void DistributionTemplated<distroClass>::AddRecord(void)
{
   values_record.push_back(_value);
   weights_record.push_back(_weight);
};

/*!
\author Vladimir Florinski
\date 05/17/2022
\param[in] bin Bin multi-index
\return Processed distribution value for this bin
*/
template <class distroClass>
distroClass DistributionTemplated<distroClass>::operator [](const MultiIndex& bin) const
{
// This template method should only work for <double> because of 0.0
   int lin_bin = n_bins.LinIdx(bin);

// This is the default routine that should fit most needs. The distribution is assumed to be zero in empty bins. Derived classes could override it if the weighting is not the same as the distribution, or to optimize for reduced dimensionality.
   return (counts[lin_bin] ? distro[lin_bin] / counts[lin_bin] : 0.0);
};

// Method specialization for GeoVector
template <>
GeoVector DistributionTemplated<GeoVector>::operator [](const MultiIndex& bin) const
{
   int lin_bin = n_bins.LinIdx(bin);

// This is the default routine that should fit most needs. The distribution is assumed to be zero in empty bins. Derived classes could override it if the weighting is not the same as the distribution, or to optimize for reduced dimensionality.
   return (counts[lin_bin] ? distro[lin_bin] / (double)counts[lin_bin] : gv_zeros);
};

// Method specialization for GeoMatrix
template <>
GeoMatrix DistributionTemplated<GeoMatrix>::operator [](const MultiIndex& bin) const
{
   int lin_bin = n_bins.LinIdx(bin);

// This is the default routine that should fit most needs. The distribution is assumed to be zero in empty bins. Derived classes could override it if the weighting is not the same as the distribution, or to optimize for reduced dimensionality.
   return (counts[lin_bin] ? distro[lin_bin] / (double)counts[lin_bin] : gm_zeros);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/04/2022
*/
template <class distroClass>
void DistributionTemplated<distroClass>::ResetDistribution(void)
{
// This template method should only work for <double> because of 0.0
   std::fill(distro.begin(), distro.end(), 0.0);
   std::fill(counts.begin(), counts.end(), 0  );
   n_events = 0;
};

// Method specialization for GeoVector
template <>
void DistributionTemplated<GeoVector>::ResetDistribution(void)
{
   std::fill(distro.begin(), distro.end(), gv_zeros);
   std::fill(counts.begin(), counts.end(), 0  );
   n_events = 0;
};

// Method specialization for GeoMatrix
template <>
void DistributionTemplated<GeoMatrix>::ResetDistribution(void)
{
   std::fill(distro.begin(), distro.end(), gm_zeros);
   std::fill(counts.begin(), counts.end(), 0  );
   n_events = 0;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
*/
template <class distroClass>
void DistributionTemplated<distroClass>::ResetRecords(void)
{
   values_record.clear();
   weights_record.clear();
};

/*!
\author Vladimir Florinski
\date 05/04/2022
\param[in] other Second distribution
\return Reference to this object
*/
template <class distroClass>
DistributionBase& DistributionTemplated<distroClass>::operator +=(const DistributionBase& other)
{
// We need to cast from the base class to the derived class because the non-templated DistributionBase does not have "distro"
   const DistributionTemplated<distroClass>& other_cast = dynamic_cast<const DistributionTemplated<distroClass>&>(other);
   std::transform(distro.begin(), distro.end(), other_cast.distro.begin(), distro.begin(), std::plus<distroClass>{});
   std::transform(counts.begin(), counts.end(), other_cast.counts.begin(), counts.begin(), std::plus<int>{});
   n_events += other_cast.n_events;
   return *this;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] other Second distribution
*/
template <class distroClass>
void DistributionTemplated<distroClass>::CopyRecords(const DistributionBase& other)
{
   // We need to cast from the base class to the derived class because the non-templated DistributionBase does not have "distro"
   const DistributionTemplated<distroClass>& other_cast = dynamic_cast<const DistributionTemplated<distroClass>&>(other);
   values_record.insert(values_record.end(), other_cast.values_record.begin(), other_cast.values_record.end());
   weights_record.insert(weights_record.end(), other_cast.weights_record.begin(), other_cast.weights_record.end());
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] n_records_in Total number of events
*/
template <class distroClass>
void DistributionTemplated<distroClass>::SetNRecords(int n_records_in)
{
   values_record.resize(n_records_in);
   weights_record.resize(n_records_in);
};

/*!
\author Vladimir Florinski
\date 09/13/2020
\param[out] size Size of an element of the "distro" array
\return Pointer to the distribution storage
*/
template <class distroClass>
void* DistributionTemplated<distroClass>::GetDistroAddress(size_t& size)
{
   size = sizeof(distroClass);
   return distro.data();
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[out] size Size of an element of the "weights_record" array
\return Pointer to the weights storage
*/
template <class distroClass>
void* DistributionTemplated<distroClass>::GetWeightsRecordAddress(size_t& size)
{
   size = sizeof(distroClass);
   return weights_record.data();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/17/2022
\param[in] action_in Action index from "ActionTable"
*/
template <class distroClass>
void DistributionTemplated<distroClass>::EvaluateWeight(int action_in)
{
// This template method should only work for <double> because of 0.0
// TODO Implement exception throwing for bad action index
   if ((action_in < 0) || (action_in >= ActionTable.size())) {
      PrintError(__FILE__, __LINE__, "Invalid action index", true);
      _weight = 0.0;
      return;
   };
   ActionTable[action_in]();
};

// Method specialization for GeoVector
template <>
void DistributionTemplated<GeoVector>::EvaluateWeight(int action_in)
{
// TODO Implement exception throwing for bad action index
   if ((action_in < 0) || (action_in >= ActionTable.size())) {
      PrintError(__FILE__, __LINE__, "Invalid action index", true);
      _weight = gv_zeros;
      return;
   };
   ActionTable[action_in]();
};

// Method specialization for GeoMatrix
template <>
void DistributionTemplated<GeoMatrix>::EvaluateWeight(int action_in)
{
// TODO Implement exception throwing for bad action index
   if ((action_in < 0) || (action_in >= ActionTable.size())) {
      PrintError(__FILE__, __LINE__, "Invalid action index", true);
      _weight = gm_zeros;
      return;
   };
   ActionTable[action_in]();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2022
\param[in] file_name Distribution file name
*/
template <class distroClass>
void DistributionTemplated<distroClass>::Dump(const std::string& file_name) const
{
   unsigned int datalen, datasize;

   std::ofstream distfile(file_name.c_str(), std::ofstream::binary);

   datalen = class_name.length();
   distfile.write((char*)&datalen, sizeof(datalen));
   distfile.write(class_name.data(), datalen);

   distfile.write((char*)&specie, sizeof(specie));

   datasize = container.size();
   distfile.write((char*)&datasize, sizeof(datasize));
   datalen = container.length();
   distfile.write((char*)&datalen, sizeof(datalen));
   distfile.write((char*)container.data(), datalen);

// Dump the distribution's data
   distfile.write((char*)&n_events, sizeof(n_events));
   distfile.write((char*)counts.data(), counts.size() * sizeof(int));
   distfile.write((char*)distro.data(), distro.size() * sizeof(distroClass));

// Dump records
   datasize = values_record.size();
   distfile.write((char*)&datasize, sizeof(datasize));
   distfile.write((char*)values_record.data(), values_record.size() * sizeof(GeoVector));
   datasize = weights_record.size();
   distfile.write((char*)&datasize, sizeof(datasize));
   distfile.write((char*)weights_record.data(), weights_record.size() * sizeof(distroClass));
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2022
\param[in] file_name Distribution file name
*/
template <class distroClass>
void DistributionTemplated<distroClass>::Restore(const std::string& file_name)
{
   unsigned int datalen, datasize;

   std::ifstream distfile(file_name.c_str(), std::ifstream::binary);
   if (!distfile.is_open()) return;

   distfile.read((char*)&datalen, sizeof(datalen));
   std::string class_name_in;
   class_name_in.resize(datalen);
   distfile.read((char*)class_name_in.data(), datalen);
   if (class_name_in != class_name) return;

   distfile.read((char*)&specie, sizeof(specie));
   
   distfile.read((char*)&datasize, sizeof(datasize));
   distfile.read((char*)&datalen, sizeof(datalen));
   container.resize(datasize, datalen);
   distfile.read((char*)container.data(), datalen);

// Restore the distribution's data
   SetupDistribution(false);
   distfile.read((char*)&n_events, sizeof(n_events));
   distfile.read((char*)counts.data(), counts.size() * sizeof(int));
   distfile.read((char*)distro.data(), distro.size() * sizeof(distroClass));

// Restore records
   distfile.read((char*)&datasize, sizeof(datasize));
   values_record.resize(datasize);
   distfile.read((char*)values_record.data(), values_record.size() * sizeof(GeoVector));
   distfile.read((char*)&datasize, sizeof(datasize));
   weights_record.resize(datasize);
   distfile.read((char*)weights_record.data(), weights_record.size() * sizeof(distroClass));
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] distfile   Opened file to which to write the distribution
\param[in] sum_counts Total counts after collapsing dimensions
\param[in] sum_distro Total distro after collapsing dimensions
\param[in] phys_units Use physical units for output
*/
template <class distroClass>
void DistributionTemplated<distroClass>::PrintSumDistro(std::ofstream& distfile, long sum_counts, distroClass sum_distro, bool phys_units) const
{
// This template method should only work for <double> because of 0.0
   distfile << std::setw(20) << (sum_counts ? sum_distro / sum_counts : 0.0) * (phys_units ? unit_distro : 1.0)
            << std::setw(20) << sum_distro * (phys_units ? unit_distro : 1.0);
};

// Method specialization for GeoVector
template <>
void DistributionTemplated<GeoVector>::PrintSumDistro(std::ofstream& distfile, long sum_counts, GeoVector sum_distro, bool phys_units) const
{
   int i;

   if (sum_counts) {
      if (phys_units) {
         for (i = 0; i < 3; i++) distfile << std::setw(20) << sum_distro[i] / sum_counts * unit_distro[i];
         for (i = 0; i < 3; i++) distfile << std::setw(20) << sum_distro[i] * unit_distro[i];
      }
      else {
         for (i = 0; i < 3; i++) distfile << std::setw(20) << sum_distro[i] / sum_counts;
         for (i = 0; i < 3; i++) distfile << std::setw(20) << sum_distro[i];
      };
   }
   else {
      for (i = 0; i < 3; i++) distfile << std::setw(20) << 0.0;
      for (i = 0; i < 3; i++) distfile << std::setw(20) << 0.0;
   };
};

// Method specialization for GeoMatrix
template <>
void DistributionTemplated<GeoMatrix>::PrintSumDistro(std::ofstream& distfile, long sum_counts, GeoMatrix sum_distro, bool phys_units) const
{
   int i, j;

   if (sum_counts) {
      if (phys_units) {
         for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) distfile << std::setw(20) << sum_distro[i][j] / sum_counts * unit_distro[i][j];
         };
         for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) distfile << std::setw(20) << sum_distro[i][j] * unit_distro[i][j];
         };
      }
      else {
         for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) distfile << std::setw(20) << sum_distro[i][j] / sum_counts;
         };
         for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) distfile << std::setw(20) << sum_distro[i][j];
         };
      };
   }
   else {
      for (i = 0; i < 9; i++) distfile << std::setw(20) << 0.0;
      for (i = 0; i < 9; i++) distfile << std::setw(20) << 0.0;
   };
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] ijk        Which component to print (the remaining dimensions are collapsed)
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
template <class distroClass>
void DistributionTemplated<distroClass>::Print1D(int ijk, const std::string& dist_name, bool phys_units) const
{
   int ijk1, ijk2;
   long lin_bin, sum_counts;
   distroClass sum_distro;
   MultiIndex bin;

   if ((ijk < 0) || (ijk >= 3)) return;

   ijk1 = (ijk + 1) % 3;
   ijk2 = (ijk + 2) % 3;

   std::ofstream distfile(dist_name.c_str());
   distfile << "# Total number of events: " << n_events << std::endl << std::endl;
   
   distfile << std::setprecision(10);
   for (bin[ijk] = 0; bin[ijk] < n_bins[ijk]; bin[ijk]++) {
      sum_counts = 0;
      sum_distro = 0.0;

// Collapse the remaining two dimensions
      for (bin[ijk1] = 0; bin[ijk1] < n_bins[ijk1]; bin[ijk1]++) {
         for (bin[ijk2] = 0; bin[ijk2] < n_bins[ijk2]; bin[ijk2]++) {
            lin_bin = n_bins.LinIdx(bin);
            sum_counts += counts[lin_bin];
            sum_distro += distro[lin_bin];
         };
      };

      distfile << std::setw(20) << BinCent(ijk, bin[ijk]) * (phys_units ? unit_val[ijk] : 1.0);
      PrintSumDistro(distfile, sum_counts, sum_distro, phys_units);
      distfile << std::setw(20) << sum_counts << std::endl;
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 09/15/2022
\param[in] ijk1       Component 1 to print
\param[in] ijk2       Component 2 to print
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
template <class distroClass>
void DistributionTemplated<distroClass>::Print2D(int ijk1, int ijk2, const std::string& dist_name, bool phys_units) const
{
   int ijk;
   long lin_bin, sum_counts;
   distroClass sum_distro;
   MultiIndex bin;

   if ((ijk1 < 0) || (ijk1 >= 3) || (ijk2 < 0) || (ijk2 >= 3)) return;

   if (ijk1 == ijk2) {
      Print1D(ijk1, dist_name, phys_units);
      return;
   };

   if ((ijk1 != 0) && (ijk2 != 0)) ijk = 0;
   else if ((ijk1 != 1) && (ijk2 != 1)) ijk = 1;
   else if ((ijk1 != 2) && (ijk2 != 2)) ijk = 2;
   else return;

   std::ofstream distfile(dist_name.c_str());
   distfile << "# Total number of events: " << n_events << std::endl << std::endl;
   
   distfile << std::setprecision(10);
   for (bin[ijk1] = 0; bin[ijk1] < n_bins[ijk1]; bin[ijk1]++) {
      for (bin[ijk2] = 0; bin[ijk2] < n_bins[ijk2]; bin[ijk2]++) {
         sum_counts = 0;
         sum_distro = 0.0;

// Collapse the remaining dimension
         for (bin[ijk] = 0; bin[ijk] < n_bins[ijk]; bin[ijk]++) {
            lin_bin = n_bins.LinIdx(bin);
            sum_counts += counts[lin_bin];
            sum_distro += distro[lin_bin];
         };

         distfile << std::setw(20) << BinCent(ijk1, bin[ijk1]) * (phys_units ? unit_val[ijk1] : 1.0)
                  << std::setw(20) << BinCent(ijk2, bin[ijk2]) * (phys_units ? unit_val[ijk2] : 1.0);
         PrintSumDistro(distfile, sum_counts, sum_distro, phys_units);
         distfile << std::setw(20) << sum_counts << std::endl;
      };
   };
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 09/15/2022
\param[in] distfile   Opened file to which to write the distribution
\param[in] record     Which record to print
\param[in] phys_units Use physical units for output
*/
template <class distroClass>
void DistributionTemplated<distroClass>::PrintWeight(std::ofstream& distfile, int record, bool phys_units) const
{
// This template method should only work for <double>
   distfile << std::setw(20) << weights_record[record] * (phys_units ? unit_distro : 1.0);
};

// Method specialization for GeoVector
template <>
void DistributionTemplated<GeoVector>::PrintWeight(std::ofstream& distfile, int record, bool phys_units) const
{
   int i;
   for (i = 0; i < 3; i++) distfile << std::setw(20) << weights_record[record][i] * (phys_units ? unit_distro[i] : 1.0);
};

// Method specialization for GeoMatrix
template <>
void DistributionTemplated<GeoMatrix>::PrintWeight(std::ofstream& distfile, int record, bool phys_units) const
{
   int i, j;
   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) distfile << std::setw(20) << weights_record[record][i][j] * (phys_units ? unit_distro[i][j] : 1.0);
   };
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 08/11/2024
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
template <class distroClass>
void DistributionTemplated<distroClass>::PrintRecords(const std::string& dist_name, bool phys_units) const
{
   int i, j;
   std::ofstream distfile(dist_name.c_str());
   
   distfile << "# Total number of records: " << values_record.size() << std::endl << std::endl;
   for (i = 0; i < values_record.size(); i++) {
      for (j = 0; j < 3; j++) distfile << std::setw(20) << values_record[i][j] * (phys_units ? unit_val[j] : 1.0);
      PrintWeight(distfile, i, phys_units);
      distfile << std::endl;
   };
};

template class DistributionTemplated<double>;
template class DistributionTemplated<GeoVector>;
template class DistributionTemplated<GeoMatrix>;

};
