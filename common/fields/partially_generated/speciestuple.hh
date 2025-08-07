

#include "fields.hh"

namespace Spectrum {

/*!
\author Lucius Schoenbaum
A Fields type can hold a Species type, but this is perhaps unintuitive, so we rename Fields for this use case.
A Fields type can also hold a mixture of Species and Field, and this can be useful in some rarer cases.
For example, if a SpeciesTuple needs to hold an indicator variable, then that is possible.
*/
   template<typename ... Ts>
   using SpeciesTuple = Fields<Ts...>;

}

