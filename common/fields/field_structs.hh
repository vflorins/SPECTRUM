/*!
\file variable_structs.hh
\author Vladimir Florinski
\author Lucius Schoenbaum
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_FIELD_STRUCTS_HH
#define SPECTRUM_FIELD_STRUCTS_HH

#include "generated/field_lists.hh"
#include <common/matrix.hh>

namespace Spectrum {

/*!
\brief Scalar with a formatted name, i.e., Scalar field
\author Lucius Schoenbaum
\date 03/25/2025
*/
template <FieldId nameid, int reconstructible_ = 0, int solvable_ = 0>
struct ScalarField {

//! static fields
   static constexpr const std::string_view name = FieldNames[nameid];
    static const int reconstructible = reconstructible_;
    static const int solvable = solvable_;

//! data field
   double value;
   
   ScalarField(double value = 0):
      value(value)
   {}

   operator double&() {
      return value;
   }

   operator double() const {
      return value;
   }

//! print string
   std::string str() const {
      std::string tmp = std::string(FieldNames[nameid]) + "[" + std::to_string(value) + "]";
      return tmp;
   }

   int scalar_size() const {
      return 1;
   }

};


/*!
\brief GeoVector with a formatted name
\author Lucius Schoenbaum
\date 03/25/2025
*/
template <FieldId nameid, int reconstructible_ = 0, int solvable_ = 0>
struct VectorField : public GeoVector {

//! static fields
   static constexpr const std::string_view name = FieldNames[nameid];
   static const int reconstructible = reconstructible_;
   static const int solvable = solvable_;

   VectorField() {}

   // todo review
   VectorField(GeoVector vec):
      GeoVector(vec)
   {}

//! print string
   std::string str() const {
      std::string tmp = std::string(FieldNames[nameid]) + "[";
      // quick lazy impl
      tmp += std::to_string(x) + " ";
      tmp += std::to_string(y) + " ";
      tmp += std::to_string(z) + "]";
      return tmp;
   }

   int scalar_size() const {
      return 3;
   }

};



template <FieldId nameid, int R, int S, FieldId nameid2, int R2, int S2>
inline GeoVector operator*(const ScalarField<nameid2, R2, S2>& lhs, const VectorField<nameid, R, S>& rhs) {
   return lhs.value*static_cast<const GeoVector&>(rhs);
}

template <FieldId nameid, int R, int S, FieldId nameid2, int R2, int S2>
inline GeoVector operator*(const VectorField<nameid, R, S>& lhs, const ScalarField<nameid2, R2, S2>& rhs) {
   return rhs.value*static_cast<const GeoVector&>(lhs);
}

template <FieldId nameid, int R, int S, FieldId nameid2, int R2, int S2>
inline GeoVector operator/(const VectorField<nameid, R, S>& lhs, const ScalarField<nameid2, R2, S2>& rhs) {
   return (1.0/rhs.value)*static_cast<const GeoVector&>(lhs);
}


/*!
\brief GeoMatrix with a formatted name
\author Lucius Schoenbaum
\date 03/25/2025
*/
template <FieldId nameid, int reconstructible_ = 0, int solvable_ = 0>
struct MatrixField : public GeoMatrix {

//! static fields
   static constexpr const std::string_view name = FieldNames[nameid];
    static const int reconstructible = reconstructible_;
    static const int solvable = solvable_;

   MatrixField() {}

   // todo review
   MatrixField(GeoMatrix mat):
      GeoMatrix(mat)
   {}

//! print string
   std::string str() const {
      // quick lazy impl
      std::string tmp = std::string(FieldNames[nameid]) + "[ ";
      for (int i = 0; i < 3; ++i) {
         tmp += "[";
         for (int j = 0; j < 3; ++j) {
            tmp += std::to_string(row[i][j]) + " ";
         }
         tmp += "] ";
      }
      tmp += "]";
      return tmp;
   }

};

}

#endif
