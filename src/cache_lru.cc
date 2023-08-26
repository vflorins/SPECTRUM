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
void BlockCache::DeleteOldest(void)
{
// Nothing to do if the cache is empty
   if(!blocks.size()) return;

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
int BlockCache::AddBlock(const BlockPtrType& block)
{
// Check if the cache is full
   if(blocks.size() >= max_cache_size) DeleteOldest();

   int bidx = block->GetNode();

   if(blocks.emplace(bidx, block).second) {
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
void BlockCache::Empty(void)
{
   helper.clear();
   queue.clear();
   blocks.clear();
};

/*!
\author Vladimir Florinski
\date 01/26/2023
\param[in] pos Position to test
\return Block index or -1 if no cached block owns the position
*/
int BlockCache::PosOwner(const GeoVector& pos)
{
   int bidx;
   BlockMapIterType iter = blocks.cbegin();

// A slow method (requires testing of all blocks in the worst case)
   while(iter != blocks.cend()) {
      if(iter->second->PositionInside(pos)) {
         bidx = iter->first;
         break;
      };
      iter++;
   };
   if(iter == blocks.cend()) bidx = -1;
   else Renew(bidx);
   
   return bidx;
};

/*!
\author Vladimir Florinski
\date 01/26/2023
*/
void BlockCache::PrintAllIndices(void) const
{
   int count = 1;
   QueueIterType iter;

   for(iter = queue.cbegin(); iter != queue.cend(); iter++) {
      std::cerr << std::setw(6) << count << std::setw(12) << *iter << std::endl;
      count++;
   };
};

};
