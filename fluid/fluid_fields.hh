/*!
\file species.hh
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPECIES_HH
#define SPECTRUM_SPECIES_HH

#include "../common/compiletime_lists.hh"
#include "../common/fields.hh"


namespace Spectrum {

enum class Model {
   GasDyn,
   MHD,
   MHDGLM,
   MHDE,
   CGL,
   CGLE,
};

enum class TurbulenceModel {
   Zank6eq,
};

enum class Form {
   primitive,
   conserved,
   flux,
};

enum class Passivity {
   active,
   passive,
};

template<Model model, Form form>
class ModelFields;

template<>
class ModelFields<Model::GasDyn, Form::primitive> : public Fields<Den_t, Vel_t, Prs_t> {};

template<>
class ModelFields<Model::GasDyn, Form::conserved> : public Fields<Den_t, Mom_t, Enr_t> {};

template<>
class ModelFields<Model::GasDyn, Form::flux> : public Fields<FlxDen_t, FlxMom_t, FlxEnr_t> {};

//

template<>
class ModelFields<Model::MHD, Form::primitive> : public Fields<Den_t, Vel_t, Prs_t, Mag_t> {};

template<>
class ModelFields<Model::MHD, Form::conserved> : public Fields<Den_t, Mom_t, Enr_t, Mag_t> {};

template<>
class ModelFields<Model::MHD, Form::flux> : public Fields<FlxDen_t, FlxMom_t, FlxEnr_t, FlxMag_t> {};

//

template<>
class ModelFields<Model::MHDGLM, Form::primitive> : public Fields<Den_t, Vel_t, Prs_t, Mag_t, Glm_t> {};

template<>
class ModelFields<Model::MHDGLM, Form::conserved> : public Fields<Den_t, Mom_t, Enr_t, Mag_t, Glm_t> {};

template<>
class ModelFields<Model::MHDGLM, Form::flux> : public Fields<FlxDen_t, FlxMom_t, FlxEnr_t, FlxMag_t, FlxGlm_t> {};


/*!
\brief General purpose class storing per-fluid-species physical data defined at a spatial location.
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/19/2025
*/
template <Model model_, Form form_, Specie specie_, Passivity passivity_ = Passivity::active>
class FluidSpecie: public ModelFields<model_, form_> {

   static constexpr Model model = model_;
   static constexpr Form form = form_;
   static constexpr Specie specie = specie_;
   static constexpr Passivity passivity = passivity_;


   template <Form form2>
   int Convert() const {
      if constexpr (form2 == form) {
         // no conversion
         return 123;
      }
      else {
         if constexpr (form2 == Form::conserved) {
            // conversion
            return 234;
         }
         else if constexpr (form2 == Form::primitive) {
            // conversion
            return 345;
         }
         else if constexpr (form2 == Form::flux) {
            // conversion
            return 456;
         }
         else {
            // conversion
            return 567;
         }
      }
   }

};


/*!
\brief Provides a container combining all conserved fluid variables into single contiguous storage.
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/26/2025
\note A Fields type can hold a FluidSpecie type, but this is perhaps unintuitive, so we rename Fields for this use case.
A Fields type can also hold a mixture of FluidSpecie and Field.
For example, if a FluidFields needs to hold an indicator variable applying to
the set of species as a group, then that is possible.
*/
template<typename ... Ts>
using FluidFields = Fields<Ts...>;

}

#endif