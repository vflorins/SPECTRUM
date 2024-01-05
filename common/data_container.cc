/*!
\file data_container.cc
\brief Implements a container class for heterogeneous parameters.
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/data_container.hh"
#include "common/spatial_data.hh"
#include "common/vectors.hh"
#include "common/matrix.hh"
#include "common/turb_prop.hh"
#include "common/mpi_config.hh"
#include <iomanip>
#include <iostream>
#include <cstring>
#include <memory>

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
// Don't change anything if "other" is empty
   if(!other.total_size) return *this;

   delete[] storage;
   n_records = other.n_records;
   total_size = other.total_size;
   storage = new uint8_t[total_size];
   std::memcpy(storage, other.storage, total_size);
   return *this;
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
\return Number of records in the container
*/
size_t DataContainer::size(void) const
{
   return n_records;
};

/*!
\author Vladimir Florinski
\date 07/27/2022
\return Number of bytes in the container
*/
size_t DataContainer::length(void) const
{
   return total_size;
};

/*!
\author Vladimir Florinski
\date 07/27/2022
\return Pointer to the binary chunk
*/
const uint8_t* DataContainer::data(void) const
{
   return storage;
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
void DataContainer::Reset(void)
{
   _param = storage;
};

/*!
\author Vladimir Florinski
\param[in] data Parameter to insert
\date 07/27/2022
*/
template <typename T> void DataContainer::Insert(T data)
{
   uint8_t size = sizeof(data);

// Copy existing arguments into a new memory region
   uint8_t* params_new = new uint8_t[total_size + size + 1];
   std::memcpy(params_new, storage, total_size);
   params_new[total_size] = size;
   std::memcpy(params_new + total_size + 1, &data, size);

// Point "storage" to the new memory region
   delete[] storage;
   storage = params_new;
   n_records++;
   total_size += size + 1;
};

/*!
\author Vladimir Florinski
\param[in] data Parameter to insert
\date 03/22/2023
*/
template <> void DataContainer::Insert(std::string data)
{
   uint8_t size = data.size();

// Copy existing arguments into a new memory region
   uint8_t* params_new = new uint8_t[total_size + size + 1];
   std::memcpy(params_new, storage, total_size);
   params_new[total_size] = size;
   std::memcpy(params_new + total_size + 1, data.data(), size);

// Point "storage" to the new memory region
   delete[] storage;
   storage = params_new;
   n_records++;
   total_size += size + 1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\param[in] data Parameter to insert
\date 07/27/2022
*/
template <typename T> void DataContainer::Insert(std::vector<T> data)
{
   uint8_t size = data.size() * sizeof(T);

// Copy existing arguments into a new memory region
   uint8_t* params_new = new uint8_t[total_size + size + 1];
   std::memcpy(params_new, storage, total_size);
   params_new[total_size] = size;
   std::memcpy(params_new + total_size + 1, data.data(), size);

// Point "storage" to the new memory region
   delete[] storage;
   storage = params_new;
   n_records++;
   total_size += size + 1;
};

/*!
\author Vladimir Florinski
\param[out] data_ptr Pointer to the parameter to be read
\date 07/27/2022
*/
template <typename T> void DataContainer::Read(T* data_ptr)
{
   std::memcpy(data_ptr, _param + 1, *_param);
   _param += *_param + 1;
};

/*!
\author Vladimir Florinski
\param[out] data_ptr Pointer to the parameter to be read
\date 03/22/2023
*/
template <> void DataContainer::Read(std::string* data_ptr)
{
   data_ptr->assign((char*)(_param + 1), *_param);
   _param += *_param + 1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\param[out] data_ptr Pointer to the parameter to be read
\date 07/27/2022
*/
template <typename T> void DataContainer::Read(std::vector<T>* data_ptr)
{
   data_ptr->resize(*_param / sizeof(T));
   std::memcpy(data_ptr->data(), _param + 1, *_param);
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
   for(size_t i = 0; i < n_records; i++) {
      for(uint8_t j = 1; j <= *_param; j++) std::cout << std::setw(4) << _param[j];
      std::cout << std::endl;
      _param += *_param + 1;
   };
};

template void DataContainer::Insert <bool>(bool data);
template void DataContainer::Insert <char>(char data);
template void DataContainer::Insert <int>(int data);
template void DataContainer::Insert <int>(std::vector<int> data);
template void DataContainer::Insert <double>(double data);
template void DataContainer::Insert <double>(std::vector<double> data);
template void DataContainer::Insert <MultiIndex>(MultiIndex data);
template void DataContainer::Insert <GeoVector>(GeoVector data);
template void DataContainer::Insert <GeoMatrix>(GeoMatrix data);
template void DataContainer::Insert <TurbProp>(TurbProp data);
template void DataContainer::Insert <std::shared_ptr<MPI_Config>*>(std::shared_ptr<MPI_Config>* data);
template void DataContainer::Read <char>(char* data_ptr);
template void DataContainer::Read <bool>(bool* data_ptr);
template void DataContainer::Read <int>(int* data_ptr);
template void DataContainer::Read <int>(std::vector<int>* data_ptr);
template void DataContainer::Read <double>(double* data_ptr);
template void DataContainer::Read <double>(std::vector<double>* data_ptr);
template void DataContainer::Read <MultiIndex>(MultiIndex* data_ptr);
template void DataContainer::Read <GeoVector>(GeoVector* data_ptr);
template void DataContainer::Read <GeoMatrix>(GeoMatrix* data_ptr);
template void DataContainer::Read <TurbProp>(TurbProp* data_ptr);
template void DataContainer::Read <std::shared_ptr<MPI_Config>*>(std::shared_ptr<MPI_Config>** data_ptr);

};
