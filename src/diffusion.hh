//
// Created by Lucius Schoenbaum on 10/1/25.
//

#ifndef SPECTRUM_DIFFUSION_HH
#define SPECTRUM_DIFFUSION_HH


#include "common/compiletime_lists.hh"
#include "common/fields.hh"

#include "diffusion_other.hh"


namespace Spectrum {


template<typename HConfig>
using DiffusionList = Fields<
FConfig<>,
DiffusionNone<HConfig>,
DiffusionIsotropicConstant<HConfig>,
DiffusionParaConstant<HConfig>,
DiffusionPerpConstant<HConfig>,
DiffusionQLTConstant<HConfig>,
DiffusionWNLTConstant<HConfig>,
DiffusionWNLTRampVLISM<HConfig>,
DiffusionFlowMomentumPowerLaw<HConfig>,
DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>,
DiffusionRigidityMagneticFieldPowerLaw<HConfig>,
DiffusionStraussEtAl2013<HConfig>,
DiffusionGuoEtAl2014<HConfig>,
DiffusionPotgieterEtAl2015<HConfig>,
DiffusionEmpiricalSOQLTandUNLT<HConfig>
>;


template<typename HConfig>
using Diffusion = FieldOps::Nth<DiffusionList<HConfig>, reinterpret_cast<int>(HConfig::DiffusionConfig::diffusionid)>;


}



#endif


