/*!
\file server_base.cc
\brief Implements a base class of a data server from an external source
\author Vladimir Florinski
*/

#include "server_base.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/25/2022
\param[in] mpi_config_in A pointer to an MPI_Config object
*/
void ServerBase::ConnectMPIConfig(const std::shared_ptr<MPI_Config> mpi_config_in)
{
   mpi_config = mpi_config_in;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
*/
void ServerBase::ServerStart(void)
{
// Set up MPI data type for "Inquiry"
   MPI_Datatype inquiry_types[] = {MPI_INT, MPI_INT, MPI_DOUBLE};
   int inquiry_lengths[] = {1, 1, 3};
   MPI_Aint inquiry_displ[3];

// Figure out field displacements using "_inquiry" as template
   MPI_Get_address(&_inquiry.type, &inquiry_displ[0]);
   MPI_Get_address(&_inquiry.node, &inquiry_displ[1]);
   MPI_Get_address(&_inquiry.pos , &inquiry_displ[2]);
   for(auto i = 2; i >= 0; i--) inquiry_displ[i] -= inquiry_displ[0];

// Commit the type
   MPI_Type_create_struct(3, inquiry_lengths, inquiry_displ, inquiry_types, &MPIInquiryType);
   MPI_Type_commit(&MPIInquiryType);
};

/*!
\author Vladimir Florinski
\date 10/27/2022
*/
void ServerBase::ServerFinish(void)
{
   MPI_Type_free(&MPIInquiryType);
};

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
int ServerBase::ServerFunctions(void)
{
   return 0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBaseFront methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
void ServerBaseFront::ServerStart(void)
{
   ServerBase::ServerStart();
   cache_line.Empty();
};

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
void ServerBaseFront::ServerFinish(void)
{
   MPI_Send(nullptr, 0, MPI_BYTE, 0, tag_stopserve, mpi_config->node_comm);
   ServerBase::ServerFinish();
};

/*!
\author Vladimir Florinski
\date 10/28/2022
*/
int ServerBaseFront::ServerFunctions(void)
{
   return 0;
};

/*!
\author Vladimir Florinski
\date 06/19/2020
\return Number of blocks in the cache
*/
int ServerBaseFront::GetNCachedBlocks(void) const
{
   return cache_line.size();
};

/*!
\author Vladimir Florinski
\date 01/27/2023
*/
void ServerBaseFront::InvalidateCache(void)
{
   cache_line.Empty();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerBaseBack methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/26/2022
\param[in] file_name_pattern_in A string describing the file naming pattern
*/
ServerBaseBack::ServerBaseBack(const std::string& file_name_pattern_in) 
              : ServerBase(),
                file_name_pattern(file_name_pattern_in)
{
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
void ServerBaseBack::ServerStart(void)
{
   ServerBase::ServerStart();

// Indices of processes returned by "Testsome"
   index_needblock   = new int[mpi_config->node_comm_size];
   index_needstencil = new int[mpi_config->node_comm_size];
   index_needvars    = new int[mpi_config->node_comm_size];
   index_stopserve   = new int[mpi_config->node_comm_size];

// Request arrays
   req_needblock   = new MPI_Request[mpi_config->node_comm_size];
   req_needstencil = new MPI_Request[mpi_config->node_comm_size];
   req_needvars    = new MPI_Request[mpi_config->node_comm_size];
   req_stopserve   = new MPI_Request[mpi_config->node_comm_size];

// Message buffers
   buf_needblock   = new Inquiry[mpi_config->node_comm_size];
   buf_needstencil = new Inquiry[mpi_config->node_comm_size];
   buf_needvars    = new Inquiry[mpi_config->node_comm_size];

// Post initial receives for all request types
   req_needblock[0] = MPI_REQUEST_NULL;
   req_needstencil[0] = MPI_REQUEST_NULL;
   req_needvars[0] = MPI_REQUEST_NULL;
   req_stopserve[0] = MPI_REQUEST_NULL;
   for(int cpu = 1; cpu < mpi_config->node_comm_size; cpu++) {
      MPI_Irecv(&buf_needblock[cpu], 1, MPIInquiryType, cpu, tag_needblock, mpi_config->node_comm, &req_needblock[cpu]);
      MPI_Irecv(&buf_needstencil[cpu], 1, MPIInquiryType, cpu, tag_needstencil, mpi_config->node_comm, &req_needstencil[cpu]);
      MPI_Irecv(&buf_needvars[cpu], 1, MPIInquiryType, cpu, tag_needvars, mpi_config->node_comm, &req_needvars[cpu]);
      MPI_Irecv(nullptr, 0, MPI_BYTE, cpu, tag_stopserve, mpi_config->node_comm, &req_stopserve[cpu]);
   };
};

/*!
\author Vladimir Florinski
\date 10/26/2022
*/
void ServerBaseBack::ServerFinish(void)
{
// Deallocate index arrays
   delete[] index_needblock;
   delete[] index_needstencil;
   delete[] index_needvars;
   delete[] index_stopserve;

// Deallocate request arrays
   delete[] req_needblock;
   delete[] req_needstencil;
   delete[] req_needvars;
   delete[] req_stopserve;
   
// Deallocate buffers
   delete[] buf_needblock;
   delete[] buf_needstencil;
   delete[] buf_needvars;

   ServerBase::ServerFinish();
};

/*!
\author Vladimir Florinski
\date 10/27/2022
\return Unused
*/
int ServerBaseBack::ServerFunctions()
{
   return 0;
};

/*!
\author Vladimir Florinski
\date 06/19/2020
\param[in] pos Position
*/
/*
void AMR_Server::TestInterpolationStencil(const GeoVector& pos)
{
   int iz, xyz, outcome;
   double tot_weight = 0.0;
   GeoVector coord = gv_zeros;
   InterpolationStencil stencil;

   std::cerr << std::endl;
   outcome = BuildInterpolationStencil(pos, stencil);
   if(outcome == -1) {
      std::cerr << "Plane interpolation failed due to block level changes at edges\n";
      return;
   };

   std::cerr << "Plane interpolation stencil:" << std::endl;
   for(iz = 0; iz < 8; iz++) {
      std::cerr << "Node: " << std::setw(8) << cache_line[stencil.blocks[iz]].GetNode();
      std::cerr << "     Zone: " << std::setw(3) << stencil.zones[iz].i
                                 << std::setw(3) << stencil.zones[iz].j
                                 << std::setw(3) << stencil.zones[iz].k;
      std::cerr << "     Weight: " << std::fixed << std::setw(9) << std::setprecision(5) << stencil.weights[iz];
      std::cerr << std::endl;
      tot_weight += stencil.weights[iz];
      for(xyz = 0; xyz < 3; xyz++) coord[xyz] += stencil.weights[iz] * cache_line[stencil.blocks[iz]].GetValue(stencil.zones[iz], xyz);
   };
   std::cerr << "Total weight is: " << std::fixed << std::setw(9) << std::setprecision(5) << tot_weight << std::endl;


   std::cerr << "Exact coordinates are: " << std::setw(10) << std::setprecision(5) << pos[0]
                                          << std::setw(10) << std::setprecision(5) << pos[1]
                                          << std::setw(10) << std::setprecision(5) << pos[2]
                                          << std::endl;
   std::cerr << "Inter coordinates are: " << std::setw(10) << std::setprecision(5) << coord[0]
                                          << std::setw(10) << std::setprecision(5) << coord[1]
                                          << std::setw(10) << std::setprecision(5) << coord[2]
                                          << std::endl;
   std::cerr << std::endl;

};
*/
/*!
\author Vladimir Florinski
\date 06/19/2020
\param[in] pos Position
*/
/*
void AMR_Server::TestInterpolationStencil2(const GeoVector& pos)
{
   int iz;
   double tot_weight = 0.0;
   InterpolationStencil stencil;

   std::cerr << std::endl;
   LoadInterpolationStencil(pos, stencil);

   std::cerr << "External interpolation stencil:" << std::endl;
   for(iz = 0; iz < 8; iz++) {
      std::cerr << "Node: " << std::setw(8) << cache_line[stencil.blocks[iz]].GetNode();
      std::cerr << "     Zone: " << std::setw(3) << stencil.zones[iz].i
                                 << std::setw(3) << stencil.zones[iz].j
                                 << std::setw(3) << stencil.zones[iz].k;
      std::cerr << "     Weight: " << std::fixed << std::setw(9) << std::setprecision(5) << stencil.weights[iz];
      std::cerr << std::endl;
      tot_weight += stencil.weights[iz];
   };
   std::cerr << "Total weigh is: " << std::fixed << std::setw(9) << std::setprecision(5) << tot_weight << std::endl;
};
*/
/*!
\author Vladimir Florinski
\date 06/22/2020
\param[in] n_pts   Number of random points
\param[in] pos_min Smallest position in the cube
\param[in] pos_max Largest position in the cube
*/
/*
void AMR_Server::TestInterpolationRandom(int n_pts, const GeoVector& pos_min, const GeoVector& pos_max)
{
   int xyz, iz, ipt, outcome;
   double tot_weight;
   InterpolationStencil stencil;
   GeoVector pos, coord, box_size;

   std::cerr << "  Point    Error   Outcome   Weight    Blocks\n";
   std::cerr << std::fixed;

   box_size = pos_max - pos_min;
   for(ipt = 0; ipt < n_pts; ipt++) {
      for(xyz = 0; xyz < 3; xyz++) pos[xyz] = pos_min[xyz] + drand48() * box_size[xyz];
      outcome = BuildInterpolationStencil(pos, stencil);
      coord = gv_zeros;
      tot_weight = 0.0;
      if(outcome == -1) tot_weight = 0.0;
      else {
         for(iz = 0; iz < 8; iz++) {
            tot_weight += stencil.weights[iz];
            for(xyz = 0; xyz < 3; xyz++) coord[xyz] += stencil.weights[iz] * cache_line[stencil.blocks[iz]].GetValue(stencil.zones[iz], xyz);
         };
      };
      AgeBlocks();
      if(!(ipt % 10)) PurgeBlocks();

      std::cerr << std::setw(5) << ipt
                << std::setw(12) << std::setprecision(6) << (pos - coord).Norm()
                << std::setw(6) << outcome
                << std::setw(12) << std::setprecision(4) << tot_weight
                << std::setw(8) << cache_line.size()
                << std::endl;
   };
};
*/
/*!
\author Vladimir Florinski
\date 06/22/2020
\param[in] n_pts   Number of random points
\param[in] pos_min Smallest position in the cube
\param[in] pos_max Largest position in the cube
*/
/*
void AMR_Server::TestInterpolationRandom2(int n_pts, const GeoVector& pos_min, const GeoVector& pos_max)
{
   int xyz, iz, ipt, outcome, block_pri;
   double tot_weight, interp1;
   InterpolationStencil stencil;
   GeoVector pos, box_size;

//   std::cerr << "  Point    Value   Outcome   Weight    Blocks\n";
   std::cerr << std::fixed;

   box_size = pos_max - pos_min;
   for(ipt = 0; ipt < n_pts; ipt++) {
      for(xyz = 0; xyz < 3; xyz++) pos[xyz] = pos_min[xyz] + drand48() * box_size[xyz];


      block_pri = CreateOnDemand(pos);
      interp1 = cache_line[block_pri].GetValue(pos, 0);



      outcome = BuildInterpolationStencil(pos, stencil);
      tot_weight = 0.0;
      interp1 = 0.0;
      for(iz = 0; iz < 8; iz++) {
         tot_weight += stencil.weights[iz];
         interp1 += stencil.weights[iz] * cache_line[stencil.blocks[iz]].GetValue(stencil.zones[iz], 0);
      };


      LoadInterpolationStencil(pos, stencil.blocks, stencil.zones, stencil.weights);
      interp1 = 0.0;
      for(iz = 0; iz < 8; iz++) {
         interp1 += stencil.weights[iz] * cache_line[stencil.blocks[iz]].GetValue(stencil.zones[iz], 0);
      };

      AgeBlocks();
      if(!(ipt % 10)) PurgeBlocks();


      std::cerr << std::setw(5) << ipt
                << std::setw(12) << std::setprecision(6) << interp1
                << std::setw(6) << outcome
                << std::setw(12) << std::setprecision(4) << tot_weight
                << std::setw(8) << cache_line.size()
                << std::endl;

   };
};
*/

};
