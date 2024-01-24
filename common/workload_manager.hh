/*!
\file workload_manager.hh
\brief Declares and implements a structure utilized to interface with the workload manager
\author Juan G Alonso Guzman
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_WORKLOAD_MANAGER_HH
#define SPECTRUM_WORKLOAD_MANAGER_HH

#ifdef USE_SLURM
#include <slurm/slurm.h>
#endif
#include <cstdlib>

namespace Spectrum {

//#define TIME_STRICT

/*!
\brief Enumeration of possible workload managers, including none (NO_WM)
\author Juan G Alonso Guzman
*/
enum Workload_Manager {NO_WM, SLURM_WM};

/*!
\brief Class to store the information about the workload manager and handle events
\author Juan G Alonso Guzman
*/
struct Workload_Manager_Handler {

//! Workload manager
   Workload_Manager workload_manager;

//! Job ID
   int jobid;

//! Default constructor
   Workload_Manager_Handler(void);

//! Detect workload manager
   void DetectManager(void);

//! Get time left in simulation (in seconds)
   long int GetRemAllocTime(void);
};

/*!
\brief Default constructor for Workload_Manager_Handler
\author Juan G Alonso Guzman
*/
inline Workload_Manager_Handler::Workload_Manager_Handler(void)
{
   workload_manager = NO_WM;
   jobid = -1;
};

/*!
\brief Detect workload manager
\author Juan G Alonso Guzman
*/
inline void Workload_Manager_Handler::DetectManager(void)
{
   char* jobid_str = NULL;

// Try SLURM
#ifdef USE_SLURM
   jobid_str = getenv("SLURM_JOB_ID");
   if(jobid_str) {
      workload_manager = SLURM_WM;
      jobid = atoi(jobid_str);
      return;
   };
#endif
//TODO: Try other workload managers, like OpenPBS 

// If none of the above job ID checks are successful, no workload manager is detected
   workload_manager = NO_WM;
   jobid = -1;
};

/*!
\brief Get time left allocated for simulation by workload manager
\author Juan G Alonso Guzman
*/
inline long int Workload_Manager_Handler::GetRemAllocTime(void)
{
   if(workload_manager == SLURM_WM) {
#ifdef USE_SLURM
      return slurm_get_rem_time(jobid);
#else
      return -1;
#endif
   }
   else return -1;
};

};

#endif
