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

namespace Spectrum {

//! Color of error messages (red)
const std::string err_color = "\033[31m";

//! Color of informational messages (yellow)
const std::string inf_color = "\033[33m";

//! Color of normal messages (default)
const std::string std_color = "\033[0m";

//! Maximum width of a line in a text file
const int line_width = 255;

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
   if(!do_print) return;
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
   if(!do_print) return;
   std::cerr << inf_color;
   std::cerr << filename << ":" << line << ": info message: " << message << std::endl;
   std::cerr << std_color;
};

/*!
\brief Returns the number of lines in a file stream
\author Vladimir Florinski
\date 07/24/2023
\param[in] inpfile Input stream
\return Number of lines in the file
*/
inline unsigned int CountLines(std::ifstream& inpfile)
{
   unsigned int count = 0;
   char line[line_width + 1];

   while(inpfile.peek() != EOF) {
      inpfile.getline(line, line_width);
      count++;
   };
   return count;
};

};

#endif
