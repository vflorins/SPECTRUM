/*!
\file model_fields.hh
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MODEL_FIELDS_HH
#define SPECTRUM_MODEL_FIELDS_HH

#include "../common/compiletime_lists.hh"
#include "../common/fields.hh"

namespace Spectrum {

template<Model model, Form form>
struct ModelFields;

#define ModelFields_variables \
using Fields::Array;       \
using Fields::Den;         \
using Fields::Vel;         \
using Fields::Prs;        \
using Fields::Mag;        \
using Fields::Mom;         \
using Fields::Enr;         \
using Fields::Glm;       \
using Fields::FlxDen;     \
using Fields::FlxMom;     \
using Fields::FlxEnr;       \
using Fields::FlxMag;     \
using Fields::FlxGlm;

#define FluidSpecie_variables \
using ModelFields::Array;       \
using ModelFields::Den;         \
using ModelFields::Vel;         \
using ModelFields::Prs;        \
using ModelFields::Mag;        \
using ModelFields::Mom;         \
using ModelFields::Enr;         \
using ModelFields::Glm;       \
using ModelFields::FlxDen;     \
using ModelFields::FlxMom;     \
using ModelFields::FlxEnr;       \
using ModelFields::FlxMag;     \
using ModelFields::FlxGlm;

template<>
struct ModelFields<Model::GasDyn, Form::primitive> : public Fields<Den_t, Vel_t, Prs_t> {ModelFields_variables};

template<>
struct ModelFields<Model::GasDyn, Form::conserved> : public Fields<Den_t, Mom_t, Enr_t> {ModelFields_variables};

template<>
struct ModelFields<Model::GasDyn, Form::flux> : public Fields<FlxDen_t, FlxMom_t, FlxEnr_t> {ModelFields_variables};

//

template<>
class ModelFields<Model::MHD, Form::primitive> : public Fields<Den_t, Vel_t, Prs_t, Mag_t> {ModelFields_variables};

template<>
class ModelFields<Model::MHD, Form::conserved> : public Fields<Den_t, Mom_t, Enr_t, Mag_t> {ModelFields_variables};

template<>
class ModelFields<Model::MHD, Form::flux> : public Fields<FlxDen_t, FlxMom_t, FlxEnr_t, FlxMag_t> {ModelFields_variables};

//

template<>
class ModelFields<Model::MHDGLM, Form::primitive> : public Fields<Den_t, Vel_t, Prs_t, Mag_t, Glm_t> {ModelFields_variables};

template<>
class ModelFields<Model::MHDGLM, Form::conserved> : public Fields<Den_t, Mom_t, Enr_t, Mag_t, Glm_t> {ModelFields_variables};

template<>
class ModelFields<Model::MHDGLM, Form::flux> : public Fields<FlxDen_t, FlxMom_t, FlxEnr_t, FlxMag_t, FlxGlm_t> {ModelFields_variables};

//

// etc.

}


#endif
