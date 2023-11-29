/*!
\file multi_index.hh
\brief Declares a three-component index class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MULTI_INDEX_HH
#define SPECTRUM_MULTI_INDEX_HH

#include "common/gpu_config.hh"
#include <cstring>
#include <ostream>

namespace Spectrum {

//! Size of this class (should be 12)
#define SZMI sizeof(MultiIndex)

//! A multi-index with zero components
#define mi_zeros MultiIndex(0, 0, 0)

#define mi_ones MultiIndex(1, 1, 1)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MultiIndex class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three component index class
\author Vladimir Florinski
*/
struct MultiIndex {

//! Storage
   union {
      struct {
         int i, j, k;
      };
      int ijk[3];
   };

//! Default constructor
   SPECTRUM_DEVICE_FUNC MultiIndex(void);

//! Constructor from indices
   SPECTRUM_DEVICE_FUNC constexpr MultiIndex(int i_in, int j_in, int k_in);

//! Copy constructor
   SPECTRUM_DEVICE_FUNC MultiIndex(const MultiIndex& other);

//! Sum of the indices
   SPECTRUM_DEVICE_FUNC int Sum(void) const;

//! Product of the indices
   SPECTRUM_DEVICE_FUNC long Prod(void) const;

//! Largest index
   SPECTRUM_DEVICE_FUNC int Largest(void) const;

//! Smallest index
   SPECTRUM_DEVICE_FUNC int Smallest(void) const;

//! Compute a linear index in a 3D array using a second multi-index
   SPECTRUM_DEVICE_FUNC long LinIdx(const MultiIndex& other) const;

//! Switch the first and third index
   SPECTRUM_DEVICE_FUNC void Flip(void);

//! Access to components for reading
   SPECTRUM_DEVICE_FUNC const int& operator [](int xyz) const;

//! Access to components for writing
   SPECTRUM_DEVICE_FUNC int& operator [](int xyz);

//! Assignment operator
   SPECTRUM_DEVICE_FUNC MultiIndex& operator =(const MultiIndex& other);

//! Set all three indices to the given value
   SPECTRUM_DEVICE_FUNC MultiIndex& operator =(int val);

//! Increment all three indices by one
   SPECTRUM_DEVICE_FUNC MultiIndex& operator ++(void);

//! Decrement all three indices by one
   SPECTRUM_DEVICE_FUNC MultiIndex& operator --(void);

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Add two multi-indices together
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator +(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Add a constant (on the right) to all three indices
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator +(const MultiIndex& midx_l, int sclr_r);

//! Add a constant (on the left) to all three indices
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator +(int sclr_l, const MultiIndex& midx_r);

//! Subtract one multi-index from another
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator -(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Subtract a constant from all three indices
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator -(const MultiIndex& midx_l, int sclr_r);

//! Multiply a multi-index by a constant from the right
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator *(const MultiIndex& midx_l, int sclr_r);

//! Multiply a multi-index by a constant from the left
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator *(int sclr_l, const MultiIndex& midx_r);

//! Divide a multi-index by a constant
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator /(const MultiIndex& midx_l, int sclr_r);

//! Modulo division of one multi-index by another
   SPECTRUM_DEVICE_FUNC friend MultiIndex operator %(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Determine whether the first multi-index is larger than the second
   SPECTRUM_DEVICE_FUNC friend bool operator >(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Determine whether the first multi-index is smaller than the second
   SPECTRUM_DEVICE_FUNC friend bool operator <(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Determine whether the first multi-index is larger of equial to the second
   SPECTRUM_DEVICE_FUNC friend bool operator >=(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Determine whether the first multi-index is smaller of equial to the second
   SPECTRUM_DEVICE_FUNC friend bool operator <=(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Determine whether two multi-indices are equal
   SPECTRUM_DEVICE_FUNC friend bool operator ==(const MultiIndex& midx_l, const MultiIndex& midx_r);

//! Determine whether two multi-indices are not equal
   SPECTRUM_DEVICE_FUNC friend bool operator !=(const MultiIndex& midx_l, const MultiIndex& midx_r);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MultiIndex inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/05/2020
*/
SPECTRUM_DEVICE_FUNC inline MultiIndex::MultiIndex(void)
{
};

/*!
\author Vladimir Florinski
\date 05/05/2020
\param[in] i_in First index
\param[in] j_in Second index
\param[in] k_in Third index
*/
SPECTRUM_DEVICE_FUNC inline constexpr MultiIndex::MultiIndex(int i_in, int j_in, int k_in)
                                                : i(i_in),
                                                  j(j_in),
                                                  k(k_in)
{
};

/*!
\author Vladimir Florinski
\date 05/10/2022
\param[in] other Multi-index to create a copy of
*/
SPECTRUM_DEVICE_FUNC inline MultiIndex::MultiIndex(const MultiIndex& other)
{
   std::memcpy(ijk, other.ijk, 3 * sizeof(int));
};

/*!
\author Vladimir Florinski
\date 05/05/2020
\return Sum of components
*/
SPECTRUM_DEVICE_FUNC inline int MultiIndex::Sum(void) const
{
   return i + j + k;
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\return Product of components
*/
SPECTRUM_DEVICE_FUNC inline long MultiIndex::Prod(void) const
{
   return i * j * k;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\return Largest index of the three
*/
SPECTRUM_DEVICE_FUNC inline int MultiIndex::Largest(void) const
{
   return (i > j ? (i > k ? i : k) : (j > k ? j : k));
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\return Smallest index of the three
*/
SPECTRUM_DEVICE_FUNC inline int MultiIndex::Smallest(void) const
{
   return (i < j ? (i < k ? i : k) : (j < k ? j : k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] other A multi-index corresponding to a triplet of coordinate indices
\return Linear index in a 3D contiguous array
*/
SPECTRUM_DEVICE_FUNC inline long MultiIndex::LinIdx(const MultiIndex& other) const
{
   return k * (j * other.i + other.j) + other.k;
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
SPECTRUM_DEVICE_FUNC inline void MultiIndex::Flip(void)
{
   int l = k;
   k = i;
   i = l;
};

/*!
\author Vladimir Florinski
\date 06/10/2020
\param[in] i The desired component
\return Component of the multi-index
*/
SPECTRUM_DEVICE_FUNC inline const int& MultiIndex::operator [](int xyz) const
{
   return ijk[xyz];
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\param[in] i The desired component
\return Component of the multi-index
*/
SPECTRUM_DEVICE_FUNC inline int& MultiIndex::operator [](int xyz)
{
   return ijk[xyz];
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\param[in] other A multi-index that will be copied
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline MultiIndex& MultiIndex::operator =(const MultiIndex& other)
{
   if(this != &other) std::memcpy(ijk, other.ijk, 3 * sizeof(int));
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] val A number to be assigned to all three components
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline MultiIndex& MultiIndex::operator =(int val)
{
   i = j = k = val;
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline MultiIndex& MultiIndex::operator ++(void)
{
   i++;
   j++;
   k++;
   return *this;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\return Reference to this object
*/
SPECTRUM_DEVICE_FUNC inline MultiIndex& MultiIndex::operator --(void)
{
   i--;
   j--;
   k--;
   return *this;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Stream insertion operator (host only)
std::ostream& operator <<(std::ostream& os, const MultiIndex& midx_r);

};

#endif
