/*!
\file cache_lru.hh
\brief Defines an LRU cache class to store multiple blocks
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CACHE_LRU_HH
#define SPECTRUM_CACHE_LRU_HH

#include <unordered_map>
#include <list>
#include <string>
#include <memory>
#include "block_base.hh"

namespace Spectrum {

//! The size of the cache.
const int max_cache_size = 100;

typedef std::shared_ptr<BlockBase> BlockPtrType;
typedef std::unordered_map<int, BlockPtrType> BlockMapType;
typedef BlockMapType::const_iterator BlockMapIterType;
typedef std::list<int> QueueType;
typedef QueueType::const_iterator QueueIterType;
typedef std::unordered_map<int, QueueIterType> HelperMapType;

/*!
\brief Cache with Least Recently Used (LRU) deletion policy
\author Vladimir Florinski
*/
class BlockCache {

protected:

//! Shared pointers to blocks
   BlockMapType blocks;
   
//! List of block IDs sorted by access time
   QueueType queue;

//! Helper map to locate the element in the queue by ID
   HelperMapType helper;

//! Renew a block
   void Renew(int bidx);

//! DeleteBlock
   void DeleteOldest(void);

public:

//! Default constructor
   BlockCache(void) = default;

//! Return the number of blocks in cache
   int size(void) const;

//! Check if the block is cached
   int Present(int bidx);

//! Add a new block
   int AddBlock(const BlockPtrType& block);

//! Empty the cache
   void Empty(void);

//! Determine which cached block owns a position
   int PosOwner(const GeoVector& pos);

//! Access element
   BlockPtrType& operator[](int bidx);

//! Print all block indices, newest first
   void PrintAllIndices(void) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockCache inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/26/2023
\return The number of blocks stored in cache
*/
inline int BlockCache::size(void) const
{
   return blocks.size();
};

/*!
\author Vladimir Florinski
\date 01/26/2023
\param[in] bidx Block index
*/
inline void BlockCache::Renew(int bidx)
{
// Remove and reinsert the block index in one operation
   queue.splice(queue.begin(), queue, helper[bidx]);
};

/*!
\author Vladimir Florinski
\date 01/26/2023
\param[in] bidx Block index
\return Same block index or -1 if block was not found
*/
inline int BlockCache::Present(int bidx)
{
   if(blocks.find(bidx) != blocks.cend()) {
      Renew(bidx);
      return bidx;
   }
   else return -1;
};

/*!
\author Vladimir Florinski
\date 01/26/2023
\param[in] bidx Block index
\return Shared pointer to the block
*/
inline BlockPtrType& BlockCache::operator[](int bidx)
{
// It is assumed the block is present; otherwise the result is indeterminate
   return blocks[bidx];
};

};

#endif
