/*!
\file fluid_specie.hh
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_FLUID_SPECIE_HH
#define SPECTRUM_FLUID_SPECIE_HH

#include "model_fields.hh"

namespace Spectrum {

/*!
\brief General purpose class storing per-fluid-species physical data defined at a spatial location.
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/19/2025
*/
template <Model model_, Form form_, SpecieId specieid, Passivity passivity_ = Passivity::active>
class FluidSpecie: public ModelFields<model_, form_> {

   using ModelFields = ModelFields<model_, form_>;
   FluidSpecie_variables;

   static constexpr Model model = model_;
   static constexpr Form form = form_;
   static constexpr Specie specie = Specie<specieid>();
   static constexpr Passivity passivity = passivity_;
   using PrimitiveState = FluidSpecie<model, Form::primitive, specie.id, passivity>;
   using ConservedState = FluidSpecie<model, Form::conserved, specie.id, passivity>;
   using FluxFunction = FluidSpecie<model, Form::flux, specie.id, passivity>;


/*!
\brief Convert among equal-sized forms
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
This method implicitly assumes the model is preserved by conversion, otherwise behavior is undefined.
 */
// todo convert in place?
   template <Form form2>
   decltype(auto) Convert() const {
      if constexpr (form2 == form) {
         // no conversion
         FluidSpecie<model, form, specie.id, passivity> out(this);
         return out;
      }
      else {
         if constexpr (model == Model::GasDyn) {
            if constexpr (form == Form::conserved && form2 == Form::primitive) {
               // conserved -> primitive
               PrimitiveState prim;
               prim.Den() = Den();
               prim.Vel() = Mom() / Den();
               prim.Pre() = Pressure<specie>(Den(), prim.Vel().Norm2(), Enr());
               return prim;
            }
            else if constexpr (form == Form::conserved && form2 == Form::flux) {
               // conserved -> flux
               double Prs = Pressure<specie>(Den(), Mom().Norm2() / Sqr(Den()), Enr());
               FluxFunction flux;
               flux.FlxDen() = Mom()[0];
               flux.FlxMom() = Mom()[0] * Mom() / Den() + Prs * gv_nx;
               flux.FlxEnr() = (Enr() + Prs) * Mom()[0] / Den();
               return flux;
            }
            else if constexpr (form == Form::primitive && form2 == Form::conserved) {
               // primitive -> conserved
               ConservedState cons;
               cons.Den() = Den();
               cons.Mom() = Den() * Vel();
               cons.Enr() = Energy<specie>(Den(), Vel().Norm2(), Prs());
               return cons;
            }
            else if constexpr (form == Form::primitive && form2 == Form::flux) {
               // primitive -> flux
               double Enr = Energy<specie>(Den(), Vel().Norm2(), Prs());
               FluxFunction flux;
               flux.FlxDen() = Den() * Vel()[0];
               flux.FlxMom() = Den() * Vel()[0] * Vel() + Prs() * gv_nx;
               flux.FlxEnr() = (Enr + Prs()) * Vel()[0];
               return flux;
            }
            else if constexpr (form == Form::flux && form2 == Form::primitive) {
               // flux -> primitive
               // TODO
               return FluidSpecie<model, Form::flux, specie.id, passivity>();
            }
            else if constexpr (form == Form::flux && form2 == Form::conserved) {
               // flux -> conserved
               // TODO
               return FluidSpecie<model, Form::flux, specie.id, passivity>();
            }
         }
         else if constexpr (model == Model::MHD) {
            if constexpr (form == Form::conserved && form2 == Form::primitive) {
               // conserved -> primitive
               PrimitiveState prim;
               prim.Den() = Den();
               prim.Vel() = Mom() / Den();
               prim.Prs() = Pressure<specie>(Den(), prim.Vel().Norm2(), Enr(), Mag().Norm2());
               prim.Mag() = Mag();
               return prim;
            }
            else if constexpr (form == Form::conserved && form2 == Form::flux) {
               // conserved -> flux
               GeoVector vel = Mom() / Den();
               double pre = Pressure<specie>(Den(), vel.Norm2(), Enr(), Mag().Norm2());
               double pre_tot = pre + Mag().Norm2() / M_8PI;
               FluxFunction flux;
               flux.FlxDen() = Mom()[0];
               flux.FlxMom() = Mom()[0] * vel + pre_tot * gv_nx - (Mag()[0] / M_4PI) * Mag();
               flux.FlxEnr() = (Enr() + pre_tot) * vel[0] - (vel * Mag()) * Mag()[0] / M_4PI;
               flux.FlxMag() = vel[0] * Mag() - Mag()[0] * vel;
               return flux;
            }
            else if constexpr (form == Form::primitive && form2 == Form::conserved) {
               // primitive -> conserved
               ConservedState cons;
               cons.Den() = Den();
               cons.Mom() = Vel() * Den();
               cons.Enr() = Energy<specie>(Den(), Vel().Norm2(), Prs(), Mag().Norm2());
               cons.Mag() = Mag();
               return cons;
            }
            else if constexpr (form == Form::primitive && form2 == Form::flux) {
               // primitive -> flux
               double Enr = Energy<specie>(Den(), Vel().Norm2(), Prs(), Mag().Norm2());
               double Prs_tot = Prs() + Mag().Norm2() / M_8PI;
               FluxFunction flux;
               flux.FlxDen() = Den() * Vel()[0];
               flux.FlxMom() = Den() * Vel()[0] * Vel() + Prs_tot * gv_nx - (Mag()[0] / M_4PI) * Mag();
               flux.FlxEnr() = (Enr + Prs_tot) * Vel()[0] - (Vel() * Mag()) * Mag()[0] / M_4PI;
               flux.FlxMag() = Vel()[0] * Mag() - Mag()[0] * Vel();
               return flux;
            }
            else if constexpr (form == Form::flux && form2 == Form::primitive) {
               // flux -> primitive
               // TODO
               return FluidSpecie<model, Form::flux, specie.id, passivity>();
            }
            else if constexpr (form == Form::flux && form2 == Form::conserved) {
               // flux -> conserved
               // TODO
               return FluidSpecie<model, Form::flux, specie.id, passivity>();
            }
         }
         else if constexpr (model == Model::MHDGLM) {
            if constexpr (form == Form::conserved && form2 == Form::primitive) {
               // conserved -> primitive
               PrimitiveState prim;
               prim.Den() = Den();
               prim.Vel() = Mom() / Den();
               prim.Prs() = Pressure<specie>(Den(), prim.Vel().Norm2(), Enr(), Mag().Norm2());
               prim.Mag() = Mag();
               prim.Glm() = Glm();
               return prim;
            }
            else if constexpr (form == Form::conserved && form2 == Form::flux) {
               // conserved -> flux
               GeoVector Vel = Mom() / Den();
               double Prs = Pressure<specie>(Den(), Vel.Norm2(), Enr(), Mag().Norm2());
               double Prs_tot = Prs + Mag().Norm2() / M_8PI;
               FluxFunction flux;
               flux.FlxDen() = Mom()[0];
               flux.FlxMom() = Mom()[0] * Vel + Prs_tot * gv_nx - (Mag()[0] / M_4PI) * Mag();
               flux.FlxEnr() = (Enr() + Prs_tot) * Vel[0] - (Vel * Mag()) * Mag()[0] / M_4PI;
               flux.FlxMag() = Vel[0] * Mag() - Mag()[0] * Vel;
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
               flux.FlxMag() += Glm() * gv_nx;
// FIXME should this flux be calculated in the RS?
         //   flux.glmf() = Sqr(CL_MHD_GLMSAFETY * FastestWaveNormal()) * Mag()[0];
               return flux;
            }
            else if constexpr (form == Form::primitive && form2 == Form::conserved) {
               // primitive -> conserved
               ConservedState cons;
               cons.Den() = Den();
               cons.Mom() = Vel() * Den();
               cons.Enr() = Energy<specie>(Den(), Vel().Norm2(), Prs(), Mag().Norm2());
               cons.Mag() = Mag();
               cons.Glm() = Glm();
               return cons;
            }
            else if constexpr (form == Form::primitive && form2 == Form::flux) {
               // primitive -> flux
               double enr = Energy<specie>(Den(), Vel().Norm2(), Prs(), Mag().Norm2());
               double pre_tot = Prs() + Mag().Norm2() / M_8PI;
               FluxFunction flux;
               flux.FlxDen() = Den() * Vel()[0];
               flux.FlxMom() = Den() * Vel()[0] * Vel() + pre_tot * gv_nx - (Mag()[0] / M_4PI) * Mag();
               flux.FlxEnr() = (enr + pre_tot) * Vel()[0] - (Vel() * Mag()) * Mag()[0] / M_4PI;
               flux.FlxMag() = Vel()[0] * Mag() - Mag()[0] * Vel();
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
               flux.FlxMag() += Glm() * gv_nx;
// FIXME should this flux be calculated in the RS?
         //   flux.glmf() = Sqr(CL_MHD_GLMSAFETY * FastestWaveNormal()) * Mag()[0];
               return flux;
            }
            else if constexpr (form == Form::flux && form2 == Form::primitive) {
               // flux -> primitive
               // TODO
               return FluidSpecie<model, Form::flux, specie.id, passivity>();
            }
            else if constexpr (form == Form::flux && form2 == Form::conserved) {
               // flux -> conserved
               // TODO
               return FluidSpecie<model, Form::flux, specie.id, passivity>();
            }
         }
         // etc.
         return 0;
      }
   }



/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
\return Fastest wave speed (sound)
*/
   inline double FastestWave() const {
      if constexpr (model == Model::GasDyn) {
         return SoundSpeed<specie>(Den(), Prs());
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
\return Fastest normal wave speed (sound)
*/
   inline double FastestWaveNormal() const {
      if constexpr (model == Model::GasDyn) {
         return 123;
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
\return Vector of primitive variables in a sonic state
*/
   void ToSonic() {
      if constexpr (model == Model::GasDyn) {
         return 123;
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }




/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
*/
   void Invert() {
      if constexpr (model == Model::GasDyn) {
         return 123;
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
   void FixEnergy() {
      // todo (HLLC)
      // resv_cons.FixEnergy(
      // prime_flux_rght, SR, SC); // U4
      if constexpr (model == Model::GasDyn) {
         return 123;
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }



/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
   void FixMomentum() {
      // todo (HLLC)
      if constexpr (model == Model::GasDyn) {
         return 123;
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/20/2025
\return Gas pressure
*/
   void GetPressure() {
      // todo (HLLC) - review
      if constexpr (model == Model::GasDyn) {
         return 123;
      }
      else if constexpr (model == Model::MHD) {
         return 123;
      }
      else if constexpr (model == Model::MHDGLM) {
         return 123;
      }
      return -1;
   }


};

}

#endif
