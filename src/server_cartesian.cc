/*!
\file server_server.cc
\brief Implements a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman
*/

#include "server_cartesian.hh"
#include "block_cartesian.hh"
#include "reader_cartesian.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void StencilCartesian::Print(void)
{
   for(auto iz = 0; iz < n_elements; iz++) {
      std::cerr << std::setw(10) << blocks[iz] << "  " << zones[iz] << "  " << weights[iz] << std::endl;
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
\param[out] block_ptr pointer to block type
*/
void ServerCartesian::InitializeBlockPtr(BlockBase* &block_ptr)
{
   block_ptr = new BlockCartesian();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::CreateBlockDatatype(void)
{
   int i;

// Temporary object for displacement calculations
   BlockBase* block_tmp;
   InitializeBlockPtr(block_tmp);
   n_variables = block_tmp->GetVariableCount();

// Set up MPI data type for "BlockXXXX" (without the dynamic storage)
   MPI_Datatype block_types[] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
   int block_lengths[] = {1, 3, 3};
   MPI_Aint block_start, block_displ[3];

// Figure out field displacements using "block_tmp" as template
   MPI_Get_address(block_tmp, &block_start);
   MPI_Get_address(block_tmp->GetNodeAddress(), &block_displ[0]);
   MPI_Get_address(block_tmp->GetFaceMinAddress(), &block_displ[1]);
   MPI_Get_address(block_tmp->GetFaceMaxAddress(), &block_displ[2]);
   for(i = 2; i >= 0; i--) block_displ[i] -= block_start;

// Commit the type
   MPI_Type_create_struct(3, block_lengths, block_displ, block_types, &MPIBlockType);
   MPI_Type_commit(&MPIBlockType);
   delete block_tmp;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
*/
void ServerCartesian::CreateStencilDatatype(void)
{
   int i;

// Set up MPI data type for "InterpolationStencil"
   MPI_Datatype stencil_types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
   int stencil_lengths[] = {1, 8, 24, 8, 24};
   MPI_Aint stencil_start, stencil_displ[5];

// Figure out field displacements using "stencil" for template
   MPI_Get_address(&stencil, &stencil_start);
   MPI_Get_address(&stencil.n_elements , &stencil_displ[0]);
   MPI_Get_address(&stencil.blocks     , &stencil_displ[1]);
   MPI_Get_address(&stencil.zones      , &stencil_displ[2]);
   MPI_Get_address(&stencil.weights    , &stencil_displ[3]);
   MPI_Get_address(&stencil.derivatives, &stencil_displ[4]);
   for(i = 4; i >= 0; i--) stencil_displ[i] -= stencil_start;

// Commit the type
   MPI_Type_create_struct(5, stencil_lengths, stencil_displ, stencil_types, &MPIStencilType);
   MPI_Type_commit(&MPIStencilType);
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::CreateMPIDatatypes(void)
{
// Create block datatype
   CreateBlockDatatype();
// Create stencil datatype
   CreateStencilDatatype();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::ServerStart(void)
{
// Call the parent version
   ServerBase::ServerStart();
   CreateMPIDatatypes();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesian::ServerFinish(void)
{
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);

// Call the parent version
   ServerBase::ServerFinish();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianFront methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianFront::ServerStart(void) {

// No need to call the ServerBaseFront version because it merely calls the ServerBase version
   ServerCartesian::ServerStart();

   cache_line.Empty();
   stencil_outcomes[0] = stencil_outcomes[1] = 0;
   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, mpi_config->node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, mpi_config->node_comm);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianFront::ServerFinish(void)
{
   MPI_Send(nullptr, 0, MPI_BYTE, 0, tag_stopserve, mpi_config->node_comm);
   ServerCartesian::ServerFinish();
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[out] block_new pointer to block type
*/
void ServerCartesianFront::MakeSharedBlock(BlockPtrType &block_new)
{
   block_new = std::make_shared<BlockCartesian>();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\return Index of the block in "cache_line"
*/
int ServerCartesianFront::RequestBlock(void)
{
   int bidx;
   BlockPtrType block_new;

// Test whether the block is cached. Either call will renew the block if it is present.
   if(_inquiry.type) bidx = cache_line.PosOwner(_inquiry.pos);
   else bidx = cache_line.Present(_inquiry.node);

// Block is not in the cache, request it from the server.
   if(bidx == -1) {
      MPI_Send(&_inquiry, 1, MPIInquiryType, 0, tag_needblock, mpi_config->node_comm);

// Allocate memory for block.
      MakeSharedBlock(block_new);
      
// Adjust for ghost cells if necessary
#if NUM_GHOST_CELLS > 0
      block_new->SetGhostCells(NUM_GHOST_CELLS);
#endif

// Receive the block in 4 parts (own data plus 3 dynamic arrays)
      MPI_Recv(block_new.get(), 1, MPIBlockType, 0, tag_sendblock, mpi_config->node_comm, MPI_STATUS_IGNORE);
      MPI_Recv(block_new->GetNeighborNodesAddress(), block_new->GetNeighborCount(), MPI_INT, 0,
               tag_sendblock, mpi_config->node_comm, MPI_STATUS_IGNORE);
      MPI_Recv(block_new->GetVariablesAddress(), block_new->GetVariableCount() * block_new->GetZoneCount(), MPI_DOUBLE, 0,
               tag_sendblock, mpi_config->node_comm, MPI_STATUS_IGNORE);
      MPI_Recv(block_new->GetNeighborLevelsAddress(), block_new->GetNeighborLevelCount(), MPI_INT, 0,
               tag_sendblock, mpi_config->node_comm, MPI_STATUS_IGNORE);

      block_new->ConfigureProperties();
      cache_line.AddBlock(block_new);
      bidx = block_new->GetNode();
   };

   return bidx;
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
*/
void ServerCartesianFront::RequestStencil(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2023
*/
void ServerCartesianFront::RequestVariables(double* vars)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\param[in] zone_lo   lower interpolation zone
\param[in] zone_hi   higher interpolation zone
\param[in] offset_lo lower interpolation offset
\param[in] offset_hi higher interpolation offset
*/
void ServerCartesianFront::InteriorInterpolationStencil(const MultiIndex zone_lo, const MultiIndex zone_hi,
                                                        const GeoVector offset_lo, const GeoVector offset_hi,
                                                        const GeoVector delta)
{
   stencil.zones[0] = zone_lo;
   stencil.weights[0] = offset_hi[0] * offset_hi[1] * offset_hi[2];
   stencil.derivatives[ 0] = -offset_hi[1] * offset_hi[2] / delta[0];
   stencil.derivatives[ 1] = -offset_hi[0] * offset_hi[2] / delta[1];
   stencil.derivatives[ 2] = -offset_hi[0] * offset_hi[1] / delta[2];

   stencil.zones[1] = stencil.zones[0];
   stencil.zones[1].i++;
   stencil.weights[1] = offset_lo[0] * offset_hi[1] * offset_hi[2];
   stencil.derivatives[ 3] =  offset_hi[1] * offset_hi[2] / delta[0];
   stencil.derivatives[ 4] = -offset_lo[0] * offset_hi[2] / delta[1];
   stencil.derivatives[ 5] = -offset_lo[0] * offset_hi[1] / delta[2];
   
   stencil.zones[2] = stencil.zones[0];
   stencil.zones[2].j++;
   stencil.weights[2] = offset_hi[0] * offset_lo[1] * offset_hi[2];
   stencil.derivatives[ 6] = -offset_lo[1] * offset_hi[2] / delta[0];
   stencil.derivatives[ 7] =  offset_hi[0] * offset_hi[2] / delta[1];
   stencil.derivatives[ 8] = -offset_hi[0] * offset_lo[1] / delta[2];

   stencil.zones[3] = stencil.zones[1];
   stencil.zones[3].j++;
   stencil.weights[3] = offset_lo[0] * offset_lo[1] * offset_hi[2];
   stencil.derivatives[ 9] =  offset_lo[1] * offset_hi[2] / delta[0];
   stencil.derivatives[10] =  offset_lo[0] * offset_hi[2] / delta[1];
   stencil.derivatives[11] = -offset_lo[0] * offset_lo[1] / delta[2];

   stencil.zones[4] = stencil.zones[0];
   stencil.zones[4].k++;
   stencil.weights[4] = offset_hi[0] * offset_hi[1] * offset_lo[2];
   stencil.derivatives[12] = -offset_hi[1] * offset_lo[2] / delta[0];
   stencil.derivatives[13] = -offset_hi[0] * offset_lo[2] / delta[1];
   stencil.derivatives[14] =  offset_hi[0] * offset_hi[1] / delta[2];

   stencil.zones[5] = stencil.zones[4];
   stencil.zones[5].i++;
   stencil.weights[5] = offset_lo[0] * offset_hi[1] * offset_lo[2];
   stencil.derivatives[15] =  offset_hi[1] * offset_lo[2] / delta[0];
   stencil.derivatives[16] = -offset_lo[0] * offset_lo[2] / delta[1];
   stencil.derivatives[17] =  offset_lo[0] * offset_hi[1] / delta[2];
   
   stencil.zones[6] = stencil.zones[4];
   stencil.zones[6].j++;
   stencil.weights[6] = offset_hi[0] * offset_lo[1] * offset_lo[2];
   stencil.derivatives[18] = -offset_lo[1] * offset_lo[2] / delta[0];
   stencil.derivatives[19] =  offset_hi[0] * offset_lo[2] / delta[1];
   stencil.derivatives[20] =  offset_hi[0] * offset_lo[1] / delta[2];

   stencil.zones[7] = zone_hi;
   stencil.weights[7] = offset_lo[0] * offset_lo[1] * offset_lo[2];
   stencil.derivatives[21] =  offset_lo[1] * offset_lo[2] / delta[0];
   stencil.derivatives[22] =  offset_lo[0] * offset_lo[2] / delta[1];
   stencil.derivatives[23] =  offset_lo[0] * offset_lo[1] / delta[2];

   stencil.n_elements = 8;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\param[in] pos Interpolation point position
\return Status: 0 if local interpolation was used, 1, 2, and 3 if plane interpolation was used, 4 if external interpolation was used
*/
int ServerCartesianFront::BuildInterpolationStencil(const GeoVector& pos)
{
   int pri_idx, sec_idx, xyz, iz;
   MultiIndex block_size, quadrant, node_idx, zone_lo, zone_hi;
   GeoVector offset_lo, offset_hi, delta;

// Find the block that owns the interpolation point
   _inquiry.type = 1;
   _inquiry.pos = pos;
   pri_idx = RequestBlock();
   BlockPtrType block_pri = cache_line[pri_idx];

   block_pri->GetZoneOffset(pos, zone_lo, offset_lo);
   delta = block_pri->GetZoneLength();
   block_size = block_pri->GetBlockSize();
   quadrant = block_pri->GetQuadrant(pos);

   offset_hi = 1.0 - offset_lo;
   zone_hi = zone_lo + 1;

// Default to local interpolation in the primary block
   InteriorInterpolationStencil(zone_lo, zone_hi, offset_lo, offset_hi, delta);

// Correct stencil zones and blocks for multi-block interpolation
// The weights do not change because of the uniformity of the grid
   for(iz = 0; iz < 8; iz++) {
      node_idx = mi_ones;

// Find which indices fall outside of the primary block's boundaries
      for(xyz = 0; xyz < 3; xyz++) {
// Use "previous" block
         if(stencil.zones[iz][xyz] < 0) {
            stencil.zones[iz][xyz] = block_size[xyz] - 1;
            node_idx[xyz]--;
         }
// Use "next" block
         else if(stencil.zones[iz][xyz] == block_size[xyz]) {
            stencil.zones[iz][xyz] = 0;
            node_idx[xyz]++;
         };
      };

// Assign the appropriate block for each zone
      sec_idx = block_pri->GetNeighborNode(node_idx);
      if(sec_idx == pri_idx) stencil.blocks[iz] = pri_idx;
      else {
         _inquiry.type = 0;
         _inquiry.node = sec_idx;
         stencil.blocks[iz] = RequestBlock();
      };
   };

   return 0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in]  t      Time
\param[in]  pos    Position
\param[out] spdata Fields, dmax, etc.
*/
void ServerCartesianFront::GetVariables(double t, const GeoVector& pos, SpatialData& spdata)
{
   double rho;
   int bidx, vidx, xyz, uvw;

// Request the block and the zone size
   _inquiry.type = 1;
   _inquiry.pos = pos;
   bidx = RequestBlock();

   BlockPtrType block = cache_line[bidx];

   double vars[n_variables] = {0.0};
   double grads[n_variables][3] = {0.0};
   spdata.dmax = block->GetZoneLength().Smallest();

//----------------------------------------------------------------------------------------------------------------------------------------------------

#if INTERP_ORDER == -1

   RequestVariables(vars);
   stencil_outcomes[1]++;

#ifdef VAR_INDEX_RHO
   rho = vars[VAR_INDEX_RHO];
#endif

   for(xyz = 0; xyz < 3; xyz++) {

#if defined(VAR_INDEX_MOM) && defined(VAR_INDEX_RHO)
      spdata.Uvec[xyz] = vars[VAR_INDEX_MOM + xyz] / rho;
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = vars[VAR_INDEX_MAG + xyz];
   };
   spdata.region = 1.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------

#elif INTERP_ORDER == 0

   MultiIndex zone;

// Take the nearest cell value
   zone = block->GetZone(pos);

#ifdef VAR_INDEX_RHO
   rho = block->GetValue(zone, VAR_INDEX_RHO);
#endif

// Convert the variables to SPECTRUM format
   for(xyz = 0; xyz < 3; xyz++) {

#if defined(VAR_INDEX_MOM) && defined(VAR_INDEX_RHO)
      spdata.Uvec[xyz] = block->GetValue(zone, VAR_INDEX_MOM + xyz) / rho;
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = block->GetValue(zone, VAR_INDEX_MAG + xyz);
   };
   spdata.region = 1.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------

#elif INTERP_ORDER == 1

   int iz;
   double var;

// Build the stencil. This is a time consuming operation.
   stencil_status = BuildInterpolationStencil(pos);
   if(stencil_status == 4) stencil_outcomes[1]++;
   else stencil_outcomes[0]++;

// Compute the weighted average
   for(iz = 0; iz < stencil.n_elements; iz++) {
      for(vidx = 0; vidx < n_variables; vidx++) {
         var = cache_line[stencil.blocks[iz]]->GetValue(stencil.zones[iz], vidx);
         vars[vidx] += stencil.weights[iz] * var;
      };
   };

#ifdef VAR_INDEX_RHO
   rho = vars[VAR_INDEX_RHO];
#endif

// Convert the variables to SPECTRUM format
   for(xyz = 0; xyz < 3; xyz++) {

#if defined(VAR_INDEX_MOM) && defined(VAR_INDEX_RHO)
      spdata.Uvec[xyz] = vars[VAR_INDEX_MOM + xyz] / rho;
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = vars[VAR_INDEX_MAG + xyz];
   };
   spdata.region = 1.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------

#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for fields
   spdata.Uvec *= unit_velocity_server / unit_velocity_fluid;
   spdata.Bvec *= unit_magnetic_server / unit_magnetic_fluid;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[out] spdata Field greadients
*/
void ServerCartesianFront::GetGradients(SpatialData& spdata)
{
#if INTERP_ORDER == -1

   RAISE_BITS(spdata._mask, BACKGROUND_grad_FAIL);
   return;

//----------------------------------------------------------------------------------------------------------------------------------------------------

// All gradients are explicitly set to zero, and the background must not attempt to compute them using "NumericalDerivatives()"
#elif INTERP_ORDER == 0

   spdata.gradUvec = 0.0;
   spdata.gradBvec = 0.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------

#elif INTERP_ORDER == 1

   double var, rho = 0.0;
   int vidx, xyz, uvw, iz;

   double grads[n_variables][3] = {0.0};

   LOWER_BITS(spdata._mask, BACKGROUND_grad_FAIL);

// Compute the weighted average
   for(iz = 0; iz < stencil.n_elements; iz++) {
      for(vidx = 0; vidx < n_variables; vidx++) {
         var = cache_line[stencil.blocks[iz]]->GetValue(stencil.zones[iz], vidx);

#ifdef VAR_INDEX_RHO
         if(vidx == VAR_INDEX_RHO) rho += stencil.weights[iz] * var;
#endif

         for(xyz = 0; xyz < 3; xyz++) {
            grads[vidx][xyz] += stencil.derivatives[3 * iz + xyz] * var;
         };
      };
   };

// Convert the gradients to SPECTRUM format
   for(xyz = 0; xyz < 3; xyz++) {
      for(uvw = 0; uvw < 3; uvw++) {

#if defined(VAR_INDEX_MOM) && defined(VAR_INDEX_RHO)
         spdata.gradUvec[xyz][uvw] = (grads[VAR_INDEX_MOM + xyz][uvw] - spdata.Uvec[xyz] * grads[VAR_INDEX_RHO][uvw]) / rho;
#else
         spdata.gradUvec[xyz][uvw] = 0.0;
#endif

// The magnetic field must be always provided
         spdata.gradBvec[xyz][uvw] =  grads[VAR_INDEX_MAG + xyz][uvw];
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for gradients
   spdata.gradUvec *= unit_velocity_server / unit_velocity_fluid;
   spdata.gradBvec *= unit_magnetic_server / unit_magnetic_fluid;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianFront::PrintStencilOutcomes(void)
{
   std::cerr << "Stencil outcomes: " << std::setw(10) << stencil_outcomes[0] << std::setw(10) << stencil_outcomes[1] << std::endl;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianBack public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
ServerCartesianBack::ServerCartesianBack(const std::string& file_name_pattern_in)
                   : ServerBaseBack(file_name_pattern_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/01/2023
*/
void ServerCartesianBack::ReadData(const std::string data_file)
{
   ReadCartesianHeader(data_file.c_str(), line_width, 1);
   ReadCartesianData(data_file.c_str(), line_width, 0, 0);
   ReadCartesianGetDomain(domain_min.Data(), domain_max.Data());
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianBack::ServerStart(void)
{
   ServerBaseBack::ServerStart();
   CreateMPIDatatypes();
   InitializeBlockPtr(block_served);

// Adjust for ghost cells if necessary
#if NUM_GHOST_CELLS > 0
   block_served->SetGhostCells(NUM_GHOST_CELLS);
#endif

// Initialize the Cartesian library and read the data into memory
   std::string data_file = file_name_pattern + ".out";

   ReadData(data_file);
   domain_min *= unit_length_server / unit_length_fluid;
   domain_max *= unit_length_server / unit_length_fluid;

   MPI_Bcast(domain_min.Data(), 3, MPI_DOUBLE, 0, mpi_config->node_comm);
   MPI_Bcast(domain_max.Data(), 3, MPI_DOUBLE, 0, mpi_config->node_comm);
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
void ServerCartesianBack::CleanReader(void)
{
   ReadCartesianClean();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void ServerCartesianBack::ServerFinish(void)
{
   CleanReader();
   delete block_served;

// This is a repeat of the function from ServerCartesian
   MPI_Type_free(&MPIBlockType);
   MPI_Type_free(&MPIStencilType);

   ServerBaseBack::ServerFinish();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\return Number of clients that completed their tasks during this cycle
*/
int ServerCartesianBack::ServerFunctions(void)
{
// Handle "needblock" requests
   HandleNeedBlockRequests();
// Handle "stopserve" requests
   return HandleStopServeRequests();
};

/*!
\author Juan G Alonso Guzman
\date 07/28/2023
*/
void ServerCartesianBack::GetBlock(const double* pos, int* node)
{
   ReadCartesianGetNode(pos, node);
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
void ServerCartesianBack::HandleNeedBlockRequests(void)
{
   GeoVector pos_cartesian;
   int cpu, cpu_idx, count_needblock = 0;

   // Service the "needblock" requests
   MPI_Testsome(mpi_config->node_comm_size, req_needblock, &count_needblock, index_needblock, MPI_STATUSES_IGNORE);

// Load the block requested. If the requestor does not know the node, figure it out.
   for(cpu_idx = 0; cpu_idx < count_needblock; cpu_idx++) {
      cpu = index_needblock[cpu_idx];

      if(buf_needblock[cpu].type) {
         pos_cartesian = buf_needblock[cpu].pos / unit_length_server * unit_length_fluid;
         GetBlock(pos_cartesian.Data(), &buf_needblock[cpu].node);
      };

      block_served->SetNode(buf_needblock[cpu].node);
      block_served->LoadDimensions(unit_length_server);
      block_served->ConfigureProperties();
      block_served->LoadNeighbors();
      block_served->LoadVariables();

// Send the block to a worker. We use a blocking Send to ensure that the buffer can be reused.
      MPI_Send(block_served, 1, MPIBlockType, cpu, tag_sendblock, mpi_config->node_comm);
      MPI_Send(block_served->GetNeighborNodesAddress(), block_served->GetNeighborCount(),
               MPI_INT, cpu, tag_sendblock, mpi_config->node_comm);
      MPI_Send(block_served->GetVariablesAddress(), block_served->GetVariableCount() * block_served->GetZoneCount(),
               MPI_DOUBLE, cpu, tag_sendblock, mpi_config->node_comm);
      MPI_Send(block_served->GetNeighborLevelsAddress(), block_served->GetNeighborLevelCount(),
               MPI_INT, cpu, tag_sendblock, mpi_config->node_comm);

// Post the receive for the next block request from this worker
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, tag_needblock, mpi_config->node_comm, &req_needblock[cpu]);
   };
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
\return Number of clients that completed their tasks during this cycle
*/
int ServerCartesianBack::HandleStopServeRequests(void)
{
   int cpu, cpu_idx, count_stopserve = 0;

   // Service the "stopserve" requests. We assume that each worker sends a single request at the end of the simulation. 
   MPI_Testsome(mpi_config->node_comm_size, req_stopserve, &count_stopserve, index_stopserve, MPI_STATUSES_IGNORE);
   
// Cancel all "needblock" and "needstencil" receive requests from the cpus that have finished.
   for(cpu_idx = 0; cpu_idx < count_stopserve; cpu_idx++) {
      cpu = index_stopserve[cpu_idx];
      MPI_Cancel(&req_needblock[cpu]);
      MPI_Cancel(&req_needstencil[cpu]);
      MPI_Cancel(&req_needvars[cpu]);
      MPI_Request_free(&req_needblock[cpu]);
      MPI_Request_free(&req_needstencil[cpu]);
      MPI_Request_free(&req_needvars[cpu]);
   };

   return count_stopserve;
};

};
