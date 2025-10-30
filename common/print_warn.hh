/*!
\file print_warn.hh
\brief Defines simple error and warning message printing
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_PRINT_WARN_HH
#define SPECTRUM_PRINT_WARN_HH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "common/definitions.hh"

namespace Spectrum {

//! Debug level
#define GEO_DEBUG_LEVEL 3

//! Comment character
#define comment_char '#'

//! Maximum width of a line in a text file
#define line_width 1000

//! Color of error messages (red)
inline constexpr const std::string err_color = "\033[31m";

//! Color of informational messages (green)
inline constexpr const std::string inf_color = "\033[32m";

//! Color of normal messages (default)
inline constexpr const std::string std_color = "\033[0m";

/*!
\brief Print an error message in color
\author Vladimir Florinski
\date 07/22/2019
\param[in] filename Name of the source file
\param[in] line     Line in the file
\param[in] message  Error message to print
\param[in] do_print Print the message or not (based on process role)
*/
inline void PrintError(const char* filename, int line, const std::string& message, bool do_print)
{
   if (!do_print) return;
   std::cerr << err_color;
   std::cerr << filename << ":" << line << ": error: " << message << std::endl;
   std::cerr << std_color;
};

/*!
\brief Print an information message in color
\author Vladimir Florinski
\date 10/06/2021
\param[in] filename Name of the source file
\param[in] line     Line in the file
\param[in] message  Information message to print
\param[in] do_print Print the message or not (based on process role)
*/
inline void PrintMessage(const char* filename, int line, const std::string& message, bool do_print)
{
   if (!do_print) return;
   std::cerr << inf_color;
   std::cerr << filename << ":" << line << ": info message: " << message << std::endl;
   std::cerr << std_color;
};

/*!
\brief Print a pass/fail test results by comparing two numbers
\author Vladimir Florinski
\date 10/23/2025
\param[in] message Description of the test
\param[in] res1    First result
\param[in] res2    Second result
*/
inline void PrintPassFail(const std::string& message, double res1, double res2)
{
   bool cmp = fabs(res2 - res1) / (fabs(res2) + fabs(res1)) < sp_tiny;
   std::cerr << message + ": " << (cmp ? inf_color + "passed" : err_color + "failed");
   std::cerr << std_color;
   std::cerr << " (" << res1 << " / " << res2 << ")\n";
};

/*!
\brief Skip to the next non comment line
\author Vladimir Florinski
\date 07/30/2024
\param[in,out] inpfile Input stream
*/
inline void SkipComments(std::ifstream& inpfile)
{
   if (!inpfile.is_open()) return;

   char first_char;

   do {
      inpfile >> std::ws;
      first_char = inpfile.peek();
      if (first_char == comment_char) inpfile.ignore(line_width, '\n');
   } while (first_char == comment_char);
};

/*!
\brief Returns the number of records in an input file stream
\author Vladimir Florinski
\date 07/26/2024
\param[in] inpfile Input stream
\return Number of lines in the file not counting the comment lines
*/
inline unsigned long CountRecords(std::ifstream& inpfile)
{
   if (!inpfile.is_open()) return 0;

   unsigned long count = 0;
   char line[line_width + 1];

   inpfile.seekg(inpfile.beg);
   while (inpfile.peek() != EOF) {
      inpfile >> std::ws;
      inpfile.getline(line, line_width);
      if (line[0] != comment_char) count++;
   };
   return count;
};

/*!
\brief Print a general connectivity table
\author Vladimir Florinski
\date 07/22/2019
\param[in] n_nodes Total nodes
\param[in] n_nbrs  Number of neighbors per node
\param[in] n_sing  Number of singular nodes
\param[in] conn    Connectivity list
*/
inline void PrintConnectivity(int n_nodes, int n_nbrs, int n_sing, const int* const* conn)
{
   int n_nbrs_actual;

   for (auto node = 0; node < n_nodes; node++) {
      std::cerr << std::setw(8) << node << " │ ";
      n_nbrs_actual = (node < n_sing ? n_nbrs - 1 : n_nbrs);

      for (auto i = 0; i < n_nbrs_actual; i++) {
         std::cerr << std::setw(8) << conn[node][i];
      };
      std::cerr << std::endl;
   };
};

};

#endif
