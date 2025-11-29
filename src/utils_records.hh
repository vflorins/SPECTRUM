/*!
\file utils_records.hh
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_UTILS_RECORDS_HH
#define SPECTRUM_UTILS_RECORDS_HH


/*
 * Classes to provide the service of recording data during a trajectory.
 * When the trajectories are set without recording enabled, these
 * classes instantiate trivially.
 *
 * Two kinds of records are possible: magnetic field extrema,
 * and trajectory paths.
 *
 */

#include <vector>
#include <fstream>
#include <iostream>
#include "common/vectors.hh"
#include "common/definitions.hh" // LocateInArray
#include "common/physics.hh" // Mom(), Vel()

namespace Spectrum {


struct MagExtrema {

//! Minimum magnetic field magnitude along a trajectory
   double Mag_min;

//! Maximum magnetic field magnitude along a trajectory
   double Mag_max;

};


/*!
\brief A class to handle magnetic field records
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
 */
template <bool record_mag_extrema>
struct MagExtremaRecords;



/*!
\brief A class to handle trajectory records
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
 */
template<typename Coordinates, Specie specie, bool record_trajectory>
struct TrajectoryRecords;




// No records case
template <>
struct MagExtremaRecords<false> {

   MagExtrema GetMagExtrema() {
      return {0.0, 0.0};
   }

   MagExtrema GetMagInitial() {
      return {0.0, 0.0};
   }

};



// No records case
template<typename Coordinates, Specie specie>
class TrajectoryRecords<Coordinates, specie, false> {

public:

   int n_segs = 0;

   explicit TrajectoryRecords(int presize):
         n_segs(0)
   {}

   void SetStart() {
      n_segs = 0;
   }

   void Store(Coordinates& coords) {
      ++n_segs;
   }

   int Segments(void) const {
      return n_segs;
   };

   void PrintTrajectory(const std::string traj_name, bool phys_units, unsigned int output,
                        unsigned int stride, double dt_out) const
   {
      std::cerr << "Cannot print trajectory because it is not being recorded." << std::endl;
   }

   void PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride = 1) const
   {
      std::cerr << "Cannot print trajectory because it is not being recorded." << std::endl;
   }

};










template <>
struct MagExtremaRecords<true> {

//! max and min at initial (start) time
   double Mag_initial = 1e30;

//! Maximum and minimum magnetic field magnitude along a trajectory
   double Mag_min = 1e30;
   double Mag_max = -1e30;

   template <typename Fields>
   void SetStart(Fields& fields);

   template <typename Fields>
   void Store(Fields& fields);

   MagExtrema GetMagExtrema() {
      return {Mag_min, Mag_max};
   }

   MagExtrema GetMagInitial() {
      return {Mag_initial, Mag_initial};
   }

};






template<typename Coordinates, Specie specie>
struct TrajectoryRecords<Coordinates, specie, true> {
public:

   std::vector<Coordinates> traj_coords;

   explicit TrajectoryRecords(int presize);

   void SetStart();

   void Store(Coordinates& coords);

//! Return the number of segments in the trajectory
   [[nodiscard]] int Segments(void) const;

protected:

//! locate a trajectory in time from stored segments, returning index data
   void GetIdx(double t_in, int& pt, double& weight) const;

//! locate a trajectory in time from stored segments, returning position
   GeoVector GetPosition(double t_in) const;

//! locate a trajectory in time from stored segments, returning velocity
   GeoVector GetVelocity(double t_in) const;

//! locate a trajectory in time from stored segments, returning energy
   double GetEnergy(double t_in) const;

//! locate a trajectory in time from stored segments, returning traversed distance (path integral)
   double GetDistance(double t_in) const;

public:

//! generate path record in storage (print)
   void PrintTrajectory(const std::string traj_name, bool phys_units, unsigned int output,
                                                             unsigned int stride, double dt_out) const;

//! generate path record in storage (print)
   void PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride = 1) const;

};




template<typename Coordinates, Specie specie, bool record_mag_extrema, bool record_trajectory>
struct Records: public MagExtremaRecords<record_mag_extrema>, TrajectoryRecords<Coordinates, specie, record_trajectory> {

   using MagExtremaRecords = MagExtremaRecords<record_mag_extrema>;
   using TrajectoryRecords = TrajectoryRecords<Coordinates, specie, record_trajectory>;
   // methods
   using TrajectoryRecords::Segments;
   using TrajectoryRecords::PrintTrajectory;
   using TrajectoryRecords::PrintCSV;
   using MagExtremaRecords::GetMagExtrema;
   using MagExtremaRecords::GetMagInitial;

   explicit Records(int presize);

   template <typename Fields>
   void SetStart(Fields& fields);

   template <typename CoordinatesIn, typename Fields>
   void Store(CoordinatesIn& coords_in, Fields& fields);

};

};

#include "utils_records.cc"

#endif
