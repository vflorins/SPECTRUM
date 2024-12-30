/*!
\file data_container.hh
\brief Declares a container class for heterogeneous parameters.
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DATA_CONTAINER_HH
#define SPECTRUM_DATA_CONTAINER_HH

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DataContainer class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Container class for storing arbitrary parameters as a binary chunk
\author Vladimir Florinski

This is a very simple class to store an arbitrary set of parameters to pass to a function. Parameters can be inserted and read sequentially, which is all that is needed to convey the parameters. The largest byte size of a parameter is 256 bytes.
*/
class DataContainer
{
private:

//! A pointer into "storage"
   uint8_t* _param = nullptr;

protected:

//! Number of records
   size_t n_records = 0;

//! Total length of the storage in bytes
   size_t total_size = 0;

//! Binary parameter chunk. Each field is stored as a byte size length field followed by the data.
   uint8_t* storage = nullptr;

public:

//! Default constructor
   DataContainer(void) {};

//! Copy constructor
   DataContainer(const DataContainer& other);

//! Move constructor
   DataContainer(DataContainer&& other);

//! Destructor
   ~DataContainer(void);

//! Copy assignment operator
   DataContainer& operator =(const DataContainer& other);

//! Move assignment operator
   DataContainer& operator =(DataContainer&& other);

//! Return pointer to the memory
   uint8_t* data(void) const {return storage;};

//! Resize the container to import data from a binary chunk
   void resize(size_t n_records_in, size_t total_size_in);

//! Return number of records
   size_t size(void) const {return n_records;};

//! Get the total size in bytes
   size_t length(void) const {return total_size;};

//! Clear the content
   void Clear(void);

//! Set the pointer to the start of the data chunk
   void Reset(void);

//! Resize and add a new parameter at the end
   template <typename T> void Insert(const T& arg);

//! Resize and add a new std::vector parameter at the end
   template <typename T> void Insert(const std::vector<T>& arg);

//! Read a parameter from the current position given by "_param"
   template <typename T> void Read(T& arg);

//! Read a std::vector parameter from the current position given by "_param"
   template <typename T> void Read(std::vector<T>& arg);

//! Print data container info (mainly for debugging purposes)
   void Print(void);
};

};

#endif
