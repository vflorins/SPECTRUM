/*!
\file server_exceptions.hh
\brief Defines all exception classes thrown by server objects
\author Juan G Alonso Guzman
*/

#ifndef SPECTRUM_SERVER_EXCEPTIONS_HH
#define SPECTRUM_SERVER_EXCEPTIONS_HH

#include <exception>

namespace Spectrum {

/*!
\brief Exception if server functions failed
\author Juan G Alonso Guzman
*/
class ExServerError : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 23/08/2023
\return Text describing the error
*/
inline const char* ExServerError::what(void) const noexcept
{
   return "Server error";
};

};

#endif
