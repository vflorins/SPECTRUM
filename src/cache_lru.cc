/*!
\file cache_lru.cc
\brief Implements an LRU cache class to store multiple blocks
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "cache_lru.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockCache methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/26/2023
*/
template <typename HConfig>
void BlockCache<HConfig>::DeleteOldest(void)
{
// Nothing to do if the cache is empty
   if (blocks.empty()) return;

   int bidx = queue.back();
   blocks.erase(bidx);
   helper.erase(bidx);
   queue.pop_back();
};

/*!
\author Vladimir Florinski
\date 01/26/2023
\param[in] block Shared pointer to a block
\return Block index
*/
template <typename HConfig>
int BlockCache<HConfig>::AddBlock(const BlockPtrType& block)
{
// Check if the cache is full
   if (blocks.size() >= max_cache_size) DeleteOldest();

   int bidx = block->GetNode();

   if (blocks.emplace(bidx, block).second) {
      queue.push_front(bidx);
      helper.emplace(bidx, queue.cbegin());
   }
   else bidx = -1;

   return bidx;
};

/*!
\author Vladimir Florinski
\date 01/26/2023
*/
template <typename HConfig>
void BlockCache<HConfig>::Empty(void)
{
   helper.clear();
   queue.clear();
   blocks.clear();
};

/*!
\author Vladimir Florinski
\date 11/21/2023
\param[in] pos Position to test
\return Block index or -1 if no cached block owns the position
*/
template <typename HConfig>
int BlockCache<HConfig>::PosOwner(const GeoVector& pos)
{
   int bidx;

// Probe all blocks starting from the most recent
   QueueIterType iter = queue.cbegin();
   while (iter != queue.cend()) {
      bidx = *iter;
      if (blocks[bidx]->PositionInside(pos)) break;
      iter++;
   };

// If the interator is past the end, then the position is not in the cache. If the iterator is at the beginning, the newest block owns the position, so no renewal is needed.
   if (iter == queue.cend()) bidx = -1;
   else if (iter != queue.cbegin()) Renew(bidx);

   return bidx;
};

/*!
\author Vladimir Florinski
\date 01/26/2023
*/
template <typename HConfig>
void BlockCache<HConfig>::PrintAllIndices(void) const
{
   int count = 1;
   QueueIterType iter;

   for (iter = queue.cbegin(); iter != queue.cend(); iter++) {
      std::cerr << std::setw(6) << count << std::setw(12) << *iter << std::endl;
      count++;
   };
};

};
