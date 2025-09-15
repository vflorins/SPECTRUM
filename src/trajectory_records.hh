/*!
\file spatial_data.hh
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_EXTREMA_DATA_HH
#define SPECTRUM_EXTREMA_DATA_HH


#include <vector>
#include "vectors.hh"
#include "definitions.hh" // LocateInArray
#include "physics.hh" // Mom(), Vel()

namespace Spectrum {


struct MagExtrema {

//! Minimum magnetic field magnitude along a trajectory
   double Mag_min;

//! Maximum magnetic field magnitude along a trajectory
   double Mag_max;

};



template <bool record_mag_extrema>
struct MagExtremaRecords;


template<typename Coordinates, Specie specie, bool record_trajectory>
struct TrajectoryRecords;


template <>
struct MagExtremaRecords<false> {};


template <>
struct MagExtremaRecords<true> {

//! max and min at initial (start) time
   double Mag_initial = 1e30;

//! Maximum and minimum magnetic field magnitude along a trajectory
   double Mag_min = 1e30;
   double Mag_max = -1e30;

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/11/2025
*/
   template <typename Fields>
   void SetStart(Fields& fields) {
      if constexpr (Fields::AbsMag_found())
         Mag_initial = fields.AbsMag();
      else
         Mag_initial = fields.Mag().Norm();
      Mag_min = Mag_initial;
      Mag_max = Mag_initial;
   }

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/11/2025
*/
   template <typename Fields>
   void Store(Fields& fields)
   {
// We trust the user enough to assume Mag is in the set of fields.
      Mag_min = fmin(Mag_min, fields.Mag());
      Mag_max = fmax(Mag_max, fields.Mag());
   };

   MagExtrema GetMagExtrema() {
      return {Mag_min, Mag_max};
   }

};



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

   void Store() {
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

   void PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride) const
   {
      std::cerr << "Cannot print trajectory because it is not being recorded." << std::endl;
   }

};




template<typename Coordinates, Specie specie>
struct TrajectoryRecords<Coordinates, specie, true> {
public:

   std::vector<Coordinates> traj_coords;

   explicit TrajectoryRecords(int presize) {
      traj_coords.reserve(presize);
   }

   void SetStart() {
      traj_coords.clear();
   }

   void Store(Coordinates& coords) {
      traj_coords.push_back(coords);
   }

/*!
\author Vladimir Florinski
\date 12/03/2020
\return Largest index in the trajectory arrays
*/
//! Return the number of segments in the trajectory
   int Segments(void) const
   {
      return traj_coords.size() - 1;
   };

protected:


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/11/2025
\param[in]  t_in   Time point (use a negative value for trajectory end)
\param[out] pt     Nearest index on the small side
\param[out] weight Weight of the nearest point
*/
   void GetIdx(double t_in, int& pt, double& weight) const
   {
// For a negative "t_in", return an index of "-1" and weight of 0. The calling program must check for this to avoid a memory access error.
      if (t_in < 0.0) {
         pt = -1;
         weight = 0.0;
      }

// For a very large "t_in", return the next to last point, so the interpolator will use only the last point with a weight of 1. If the trajectory has only the starting point in it, the index returned will be "-1", and the calling program must check for this to avoid a memory access error.
      else if (t_in >= traj_coords.back().Time()) {
         pt = traj_coords.size() - 2;
         weight = 0.0;
      }

// For additional safety, use the last argument to cap the output of "LocateInArray()".
      else {
         pt = LocateInArray(0, traj_coords.size() - 1, traj_coords.data(), t_in, true);
         weight = (traj_coords[pt + 1].Time() - t_in) / (traj_coords[pt + 1].Time() - traj_coords[pt].Time());
      };
   };



/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Position
*/
   GeoVector GetPosition(double t_in) const
   {
      int pt;
      double weight;
      GetIdx(t_in, pt, weight);
      if (pt < 0) return traj_coords.Pos()[0];
      else return weight * traj_coords.Pos()[pt] + (1.0 - weight) * traj_coords.Pos()[pt + 1];
   };


/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Velocity
*/
   GeoVector GetVelocity(double t_in) const
   {
      int pt;
      double weight, mom1, vel1, mom2, vel2;
      GetIdx(t_in, pt, weight);
      if (pt < 0) {
         mom1 = traj_coords[0].Mom().Norm();
         vel1 = Vel<specie>(mom1);
         return (vel1 / mom1) * traj_coords[0].Mom();
      }
   // Linear interpolation between nearest frames
      else {
         auto mom1_ = traj_coords[pt].Mom();
         auto mom2_ = traj_coords[pt + 1].Mom();
         mom1 = mom1_.Norm();
         vel1 = Vel<specie>(mom1);
         mom2 = traj_coords[pt + 1].Mom().Norm();
         vel2 = Vel<specie>(mom2);
         return weight * (vel1 / mom1) * mom1_ + (1.0 - weight) * (vel2 / mom2) * mom2_;
      };
   };

/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Kinetic energy
*/
   double GetEnergy(double t_in) const
   {
      int pt;
      double weight;
      GetIdx(t_in, pt, weight);
      if (pt < 0)
         return EnrKin<specie>(traj_coords[0].Mom().Norm());
      else
         return weight * EnrKin<specie>(traj_coords[pt].Mom().Norm()) + (1.0 - weight) * EnrKin<specie>(traj_coords[pt + 1].Mom().Norm());
   };


/*!
\author Vladimir Florinski
\date 07/14/2020
\param[in] t_in Time point (use a negative value for trajectory end)
\return Integral along trajectory
*/
   double GetDistance(double t_in) const
   {
      int pt, ipt;
      double weight, length = 0.0;
      GeoVector pos_final;

      GetIdx(t_in, pt, weight);
      if (pt >= 0) {
         for (ipt = 0; ipt < pt; ipt++) length += (traj_coords[ipt + 1].Pos() - traj_coords[ipt].Pos()).Norm();
         pos_final = weight * traj_coords[pt].Pos() + (1.0 - weight) * traj_coords[pt + 1].Pos();
         length += (pos_final - traj_coords[pt].Pos()).Norm();
      };
      return length;
   };


public:

/*!
\author Vladimir Florinski
\date 07/13/2020
\param[in] traj_name  File name
\param[in] phys_units Use physical units for output
\param[in] output     Which coordinates to print
\param[in] stride     Distance between points in the output (optional). If stride = 0, output based on dt_out.
\param[in] dt_out     Time increment at which to output quantities when stride = 0
*/
   void PrintTrajectory(const std::string traj_name, bool phys_units, unsigned int output,
                                                             unsigned int stride, double dt_out) const
   {
      unsigned int pt, iter_out = 0, max_out = 1000000;
      double mom_mag, vm_ratio, t_out = 0.0, engkin_t;
      GeoVector pos_t, vel_t;
      std::ofstream trajfile;

      trajfile.open(traj_name.c_str());

// Generate multiple column output
      trajfile << std::setprecision(12);

      if (stride) {
         for (pt = 0; pt < traj_coords.size(); pt += stride) {
//FIXME: This computation of momentum magnitude is not guaranteed to work for focused transport. It is only approximately correct when magnitude (_mom[0]) >> pitch angle cosine (_mom[1]).
            Coordinates coords = traj_coords[pt];
            mom_mag = coords.Mom().Norm();
            vm_ratio = Vel<specie>(mom_mag) / mom_mag;

            if (output & 0x01) trajfile << std::setw(20) << coords.Time() * (phys_units ? unit_time_fluid : 1.0);
            if (output & 0x02) trajfile << std::setw(20) << coords.Pos()[0] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x04) trajfile << std::setw(20) << coords.Pos()[1] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x08) trajfile << std::setw(20) << coords.Pos()[2] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x10) trajfile << std::setw(20) << vm_ratio * coords.Mom()[0] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x20) trajfile << std::setw(20) << vm_ratio * coords.Mom()[1] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x40) trajfile << std::setw(20) << vm_ratio * coords.Mom()[2] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x80) trajfile << std::setw(20) << EnrKin<specie>(mom_mag) * (phys_units ? unit_energy_particle : 1.0);
            trajfile << std::endl;
         };
      }
      else {
         while (t_out < traj_coords.back().Time() && iter_out < max_out) {
            pos_t = GetPosition(t_out);
            vel_t = GetVelocity(t_out);
            engkin_t = GetEnergy(t_out);

            if (output & 0x01) trajfile << std::setw(20) << t_out * (phys_units ? unit_time_fluid : 1.0);
            if (output & 0x02) trajfile << std::setw(20) << pos_t[0] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x04) trajfile << std::setw(20) << pos_t[1] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x08) trajfile << std::setw(20) << pos_t[2] * (phys_units ? unit_length_fluid : 1.0);
            if (output & 0x10) trajfile << std::setw(20) << vel_t[0] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x20) trajfile << std::setw(20) << vel_t[1] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x40) trajfile << std::setw(20) << vel_t[2] * (phys_units ? unit_velocity_fluid : 1.0);
            if (output & 0x80) trajfile << std::setw(20) << engkin_t * (phys_units ? unit_energy_particle : 1.0);
            trajfile << std::endl;

            t_out += dt_out;
         };
      };
      trajfile.close();
   };

/*!
\author Vladimir Florinski
\date 10/15/2020
\param[in] traj_name  File name
\param[in] phys_units Use physical units for output
\param[in] stride     Distance between points in the output (optional)
*/
   void PrintCSV(const std::string traj_name, bool phys_units, unsigned int stride) const
   {
      unsigned int pt;
      std::ofstream trajfile;

      trajfile.open(traj_name.c_str());

// Generate CSV output
      trajfile << std::setprecision(12);
      for (pt = 0; pt < traj_coords.size(); pt += stride) {
         Coordinates coords = traj_coords[pt];
         trajfile << std::setw(20) << coords.Pos()[0] * (phys_units ? unit_length_fluid : 1.0);
         trajfile << ",";
         trajfile << std::setw(20) << coords.Pos()[1] * (phys_units ? unit_length_fluid : 1.0);
         trajfile << ",";
         trajfile << std::setw(20) << coords.Pos()[2] * (phys_units ? unit_length_fluid : 1.0);
         trajfile << std::endl;
      };
      trajfile.close();
   };

};







template<typename Coordinates, Specie specie, bool record_mag_extrema, bool record_trajectory>
struct Records: public MagExtremaRecords<record_mag_extrema>, TrajectoryRecords<Coordinates, specie, record_trajectory> {

   using MagExtremaRecords = MagExtremaRecords<record_mag_extrema>;
   using TrajectoryRecords = TrajectoryRecords<Coordinates, specie, record_trajectory>;
   // methods
   using TrajectoryRecords::Segments;
   using TrajectoryRecords::PrintTrajectory;
   using TrajectoryRecords::PrintCSV;

   explicit Records(int presize):
      MagExtremaRecords(),
      TrajectoryRecords(presize)
   {}

   template <typename Fields>
   void SetStart(Fields& fields) {
      if constexpr (record_mag_extrema)
         MagExtremaRecords::SetStart(fields);
      TrajectoryRecords::SetStart();
   }

   template <typename Fields>
   void Store(Coordinates& coords, Fields& fields) {
      if constexpr (record_mag_extrema)
         MagExtremaRecords::Store(fields);
      TrajectoryRecords::Store(coords);
   }

   MagExtrema GetMagExtrema() {
      if constexpr (record_mag_extrema)
         return MagExtremaRecords::GetMagExtrema();
      else
         return {0.0,0.0};
   }

};

};


#endif
