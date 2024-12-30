/*!
\file data_container.cc
\brief Implements a container class for heterogeneous parameters.
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <iomanip>
#include <iostream>

#include "common/data_container.hh"
#include "common/spatial_data.hh"
#include "common/vectors.hh"
#include "common/matrix.hh"
#include "common/turb_prop.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DataContainer methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 12/03/2020
\param[in] other Container to be copied into this one
*/
DataContainer::DataContainer(const DataContainer& other)
{
   operator =(other);
};

/*!
\author Vladimir Florinski
\date 10/16/2024
\param[in] other Container to be moved into this one
*/
DataContainer::DataContainer(DataContainer&& other)
{
   operator =(other);
};

/*!
\author Vladimir Florinski
\date 12/03/2020
*/
DataContainer::~DataContainer(void)
{
   delete[] storage;
};

/*!
\author Vladimir Florinski
\date 09/27/2021
\param[in] other Container to be copied into this one
*/
DataContainer& DataContainer::operator =(const DataContainer& other)
{
// Don't change anything if "other" is empty or self
   if (!other.total_size || (this == &other)) return *this;

   delete[] storage;
   n_records = other.n_records;
   total_size = other.total_size;
   storage = new uint8_t[total_size];
   std::memcpy(storage, other.storage, total_size);
   return *this;
};

/*!
\author Vladimir Florinski
\date 10/16/2024
\param[in] other Container to be moved into this one
*/
DataContainer& DataContainer::operator =(DataContainer&& other)
{
// Don't change anything if "other" is empty or self
   if (!other.total_size || (this == &other)) return *this;

   delete[] storage;
   n_records = other.n_records;
   total_size = other.total_size;
   storage = other.storage;
   other.n_records = 0;
   other.total_size = 0;
   other.storage = nullptr;
   return *this;
};

/*!
\author Vladimir Florinski
\date 07/27/2022
\param[in] n_records_in  Number of records
\param[in] total_size_in Container size in bytes
*/
void DataContainer::resize(size_t n_records_in, size_t total_size_in)
{
   n_records = n_records_in;
   total_size = total_size_in;
   delete[] storage;
   storage = new uint8_t[total_size];
};

/*!
\author Vladimir Florinski
\date 12/03/2020
*/
void DataContainer::Clear(void)
{
   delete[] storage;
   storage = nullptr;
   n_records = 0;
   total_size = 0;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
*/
void DataContainer::Reset(void)
{
   _param = storage;
};

/*!
\author Vladimir Florinski
\date 07/27/2022
\param[in] arg Parameter to insert
*/
template <typename T>
void DataContainer::Insert(const T& arg)
{
   uint8_t arg_size = sizeof(arg);

// Copy existing arguments into a new memory region
   uint8_t* storage_new = new uint8_t[total_size + arg_size + 1];
   std::memcpy(storage_new, storage, total_size);

// Insert the size of the new argument and the argument itself
   storage_new[total_size] = arg_size;
   std::memcpy(storage_new + total_size + 1, &arg, arg_size);

// Replace "storage" with "storage_new" and increment the size
   delete[] storage;
   storage = storage_new;
   n_records++;
   total_size += arg_size + 1;
};

/*!
\author Vladimir Florinski
\date 03/22/2023
\param[in] arg Parameter to insert
*/
template <>
void DataContainer::Insert(const std::string& arg)
{
   uint8_t arg_size = arg.size();

// Copy existing arguments into a new memory region
   uint8_t* storage_new = new uint8_t[total_size + arg_size + 1];
   std::memcpy(storage_new, storage, total_size);

// Insert the size of the new argument and the argument itself
   storage_new[total_size] = arg_size;
   std::memcpy(storage_new + total_size + 1, arg.data(), arg_size);

// Replace "storage" with "storage_new" and increment the size
   delete[] storage;
   storage = storage_new;
   n_records++;
   total_size += arg_size + 1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/27/2022
\param[in] arg Parameter to insert
*/
template <typename T>
void DataContainer::Insert(const std::vector<T>& arg)
{
   uint8_t arg_size = arg.size() * sizeof(T);

// Copy existing arguments into a new memory region
   uint8_t* storage_new = new uint8_t[total_size + arg_size + 1];
   std::memcpy(storage_new, storage, total_size);

// Insert the size of the new argument and the argument itself
   storage_new[total_size] = arg_size;
   std::memcpy(storage_new + total_size + 1, arg.data(), arg_size);

// Point "storage" to the new memory region
   delete[] storage;
   storage = storage_new;
   n_records++;
   total_size += arg_size + 1;
};

/*!
\author Vladimir Florinski
\date 11/28/2024
\param[out] arg Parameter to read
*/
template <typename T>
void DataContainer::Read(T& arg)
{
   std::memcpy(&arg, _param + 1, *_param);
   _param += *_param + 1;
};

/*!
\author Vladimir Florinski
\date 10/16/2024
\param[out] arg Parameter to read
*/
template <>
void DataContainer::Read(std::string& arg)
{
   arg.assign((char*)(_param + 1), *_param);
   _param += *_param + 1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/16/2024
\param[out] arg Parameter to read
*/
template <typename T>
void DataContainer::Read(std::vector<T>& arg)
{
   arg.resize(*_param / sizeof(T));
   std::memcpy(arg.data(), _param + 1, *_param);
   _param += *_param + 1;
};

/*!
\author Keyvan Ghanbari
\author Juan G Alonso Guzman
\date 04/23/2021
*/ 
void DataContainer::Print(void)
{
   std::cout << "Number of records: " << n_records << std::endl;
   std::cout << "Total length in bytes: " << total_size << std::endl;
   std::cout << "Data: " << std::endl;
   _param = storage;

// Print the byte decimal representation of the parameters
   std::cout << std::hex;
   std::cout << std::setfill('0');
   for (size_t i = 0; i < n_records; i++) {
      for (uint8_t j = 1; j <= *_param; j++) std::cout << " " << std::setw(2) << (int)_param[j];
      std::cout << std::endl;
      _param += *_param + 1;
   };
};

template void DataContainer::Insert<bool>(const bool& arg);
template void DataContainer::Insert<char>(const char& arg);
template void DataContainer::Insert<int>(const int& arg);
template void DataContainer::Insert<double>(const double& arg);
template void DataContainer::Insert<MultiIndex>(const MultiIndex& arg);
template void DataContainer::Insert<GeoVector>(const GeoVector& arg);
template void DataContainer::Insert<GeoMatrix>(const GeoMatrix& arg);
template void DataContainer::Insert<TurbProp>(const TurbProp& arg);

template void DataContainer::Insert<int>(const std::vector<int>& arg);
template void DataContainer::Insert<double>(const std::vector<double>& arg);

template void DataContainer::Read<char>(char& arg);
template void DataContainer::Read<bool>(bool& arg);
template void DataContainer::Read<int>(int& arg);
template void DataContainer::Read<double>(double& arg);
template void DataContainer::Read<MultiIndex>(MultiIndex& arg);
template void DataContainer::Read<GeoVector>(GeoVector& arg);
template void DataContainer::Read<GeoMatrix>(GeoMatrix& arg);
template void DataContainer::Read<TurbProp>(TurbProp& arg);

template void DataContainer::Read<int>(std::vector<int>& arg);
template void DataContainer::Read<double>(std::vector<double>& arg);

};
