/*!
\file multi_index.hh
\brief Declares a three-component index class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MULTI_INDEX_HH
#define SPECTRUM_MULTI_INDEX_HH

#include "common/simple_array.hh"

namespace Spectrum {

//! Size of this class (should be 12)
#define SZMI sizeof(MultiIndex)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MultiIndex class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three component index class
\author Vladimir Florinski
\author Juan G Alonso Guzman
*/
struct MultiIndex : public SimpleArray<int, 3>
{
   using SimpleArray::operator=;
   using SimpleArray::operator+=;
   using SimpleArray::operator-=;
   using SimpleArray::operator*=;
   using SimpleArray::operator/=;

//! Default constructor
   SPECTRUM_DEVICE_FUNC constexpr MultiIndex(void) {};

//! Constructor from a single value
   SPECTRUM_DEVICE_FUNC explicit constexpr MultiIndex(int a);

//! Constructor from an array
   SPECTRUM_DEVICE_FUNC explicit constexpr MultiIndex(const int* other);

//! Constructor from indices
   SPECTRUM_DEVICE_FUNC constexpr MultiIndex(int i_in, int j_in, int k_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC constexpr MultiIndex(const SimpleArray<int, 3>& other);

//! Compute a linear index in a 3D array using a second multi-index
   SPECTRUM_DEVICE_FUNC size_t LinIdx(const MultiIndex& other) const;

//! Switch the first and third index
   SPECTRUM_DEVICE_FUNC void Flip(void);

//! Increment all three indices by one
   SPECTRUM_DEVICE_FUNC MultiIndex& operator ++(void);

//! Decrement all three indices by one
   SPECTRUM_DEVICE_FUNC MultiIndex& operator --(void);

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Determine whether the first multi-index is larger than the second
   SPECTRUM_DEVICE_FUNC friend bool operator >(const MultiIndex& left, const MultiIndex& right);

//! Determine whether the first multi-index is smaller than the second
   SPECTRUM_DEVICE_FUNC friend bool operator <(const MultiIndex& left, const MultiIndex& right);

//! Determine whether the first multi-index is larger of equal to the second
   SPECTRUM_DEVICE_FUNC friend bool operator >=(const MultiIndex& left, const MultiIndex& right);

//! Determine whether the first multi-index is smaller of equal to the second
   SPECTRUM_DEVICE_FUNC friend bool operator <=(const MultiIndex& left, const MultiIndex& right);

//! Determine whether two multi-indices are equal
   SPECTRUM_DEVICE_FUNC friend bool operator ==(const MultiIndex& left, const MultiIndex& right);

//! Determine whether two multi-indices are not equal
   SPECTRUM_DEVICE_FUNC friend bool operator !=(const MultiIndex& left, const MultiIndex& right);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MultiIndex inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/09/2024
\param[in] a Number to be asigned to each index
*/
SPECTRUM_DEVICE_FUNC inline constexpr MultiIndex::MultiIndex(int a)
{
   i = j = k = a;
};

/*!
\author Vladimir Florinski
\date 03/10/2024
\param[in] other Array to initialize from
*/
SPECTRUM_DEVICE_FUNC inline constexpr MultiIndex::MultiIndex(const int* other)
{
   memcpy(data, other, 3 * sizeof(int));
};

/*!
\author Vladimir Florinski
\date 03/09/2024
\param[in] i_in First index
\param[in] j_in Second index
\param[in] k_in Third index
*/
SPECTRUM_DEVICE_FUNC inline constexpr MultiIndex::MultiIndex(int i_in, int j_in, int k_in)
{
   i = i_in;
   j = j_in;
   k = k_in;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline constexpr MultiIndex::MultiIndex(const SimpleArray<int, 3>& other)
{
   memcpy(data, other.data, 3 * sizeof(int));
};

/*!
\author Vladimir Florinski
\date 03/09/2024
\param[in] other A multi-index corresponding to a triplet of coordinate indices
\return Linear index in a 3D contiguous array
*/
SPECTRUM_DEVICE_FUNC inline size_t MultiIndex::LinIdx(const MultiIndex& other) const
{
   return k * (j * other.i + other.j) + other.k;
};

/*!
\author Vladimir Florinski
\date 03/09/2024
*/
SPECTRUM_DEVICE_FUNC inline void MultiIndex::Flip(void)
{
   std::swap(i, k);
};

/*!
\author Vladimir Florinski
\date 03/09/2024
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
\date 03/09/2024
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
// Friend methods of MultiIndex
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] left  Left operand
\param[in] right Right operand
\return True if all three components of "left" are larger than those of "right"
*/
SPECTRUM_DEVICE_FUNC inline bool operator >(const MultiIndex& left, const MultiIndex& right)
{
   return ((left.i > right.i) && (left.j > right.j) && (left.k > right.k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] left  Left operand
\param[in] right Right operand
\return True if all three components of "left" are smaller than those of "right"
*/
SPECTRUM_DEVICE_FUNC inline bool operator <(const MultiIndex& left, const MultiIndex& right)
{
   return ((left.i < right.i) && (left.j < right.j) && (left.k < right.k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] left  Left operand
\param[in] right Right operand
\return True if all three components of "left" are greater or equal to those of "right"
*/
SPECTRUM_DEVICE_FUNC inline bool operator >=(const MultiIndex& left, const MultiIndex& right)
{
   return ((left.i >= right.i) && (left.j >= right.j) && (left.k >= right.k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] left  Left operand
\param[in] right Right operand
\return True if all three components of "left" are smaller or equal to those of "right"
*/
SPECTRUM_DEVICE_FUNC inline bool operator <=(const MultiIndex& left, const MultiIndex& right)
{
   return ((left.i <= right.i) && (left.j <= right.j) && (left.k <= right.k));
};

/*!
\author Vladimir Florinski
\date 06/15/2021
\param[in] left  Left operand
\param[in] right Right operand
\return True if all three components of "left" are equal to those of "right"
*/
SPECTRUM_DEVICE_FUNC inline bool operator ==(const MultiIndex& left, const MultiIndex& right)
{
   return ((left.i == right.i) && (left.j == right.j) && (left.k == right.k));
};

/*!
\author Vladimir Florinski
\date 06/15/2021
\param[in] left  Left operand
\param[in] right Right operand
\return True if any of the three components of "left" is not equal to those of "right"
*/
SPECTRUM_DEVICE_FUNC inline bool operator !=(const MultiIndex& left, const MultiIndex& right)
{
   return ((left.i != right.i) || (left.j != right.j) || (left.k != right.k));
};

//! A multi-index with zero components
SPECTRUM_CONSTEXPR MultiIndex mi_zeros = {0, 0, 0};

//! A multi-index with unit components
SPECTRUM_CONSTEXPR MultiIndex mi_ones = {1, 1, 1};

};

#endif
