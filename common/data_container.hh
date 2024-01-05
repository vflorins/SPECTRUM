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
class DataContainer {

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
   DataContainer(void) = default;

//! Copy constructor
   DataContainer(const DataContainer& other);

//! Destructor
   ~DataContainer(void);

//! Assignment operator
   DataContainer& operator =(const DataContainer& other);

//! Clear the content
   void Clear(void);

//! Return number of records
   size_t size(void) const;

//! Get the total size in bytes
   size_t length(void) const;

//! Return the pointer to the data
   const uint8_t* data(void) const;

//! Resize the container to import data from a binary chunk
   void resize(size_t n_records_in, size_t total_size_in);

//! Set the pointer to the start of the data chunk
   void Reset(void);

//! Copy insertion operation with resizing (generic)
   template <typename T> void Insert(T data);

//! Copy insertion operation with resizing for vectors
   template <typename T> void Insert(std::vector<T> data);

//! Read a parameter from the current pointer position (generic)
   template <typename T> void Read(T* data_ptr);

//! Read a parameter from the current pointer position for vectors
   template <typename T> void Read(std::vector<T>* data_ptr);

//! Print data container info (mainly for debugging purposes)
   void Print(void);
};

};

#endif
