# File config.baseline.py Created by Lucius Schoenbaum October 20, 2025






backgrounds = [
    "CylindricalObstacle",
    "Dipole",
    "Discontinuity",
    "MagnetizedCylinder",
    "MagnetizedSphere",
    "Server",
    "ServerBATL",
    "ServerCartesian",
    "Shock",
    "SmoothDiscontinuity",
    "SmoothShock",
    "SolarWind",
    "SolarWindTermShock",
    "SphericalObstacle",
    "Uniform",
    "VLISMBochum",
    "Waves",
]

servers = {
    "Server": 0,
    "ServerBATL": 0,
    "ServerCartesian": 0,
}

trajectories = [
    # todo
]

diffusions = [
    "TimeUniform",
    "QLTConstant",
    # todo
]



class ParameterInfo:
    def __init__(self, name, description, default, possible_values = None, parameter_type = int):
        self.name = name
        self.description = description
        self.possible_values = possible_values
        self.default = default
        self.parameter_type = parameter_type

    def str(self, cpp_comment = True):
        out = "//! " if cpp_comment else ""
        nl = "\n// " if cpp_comment else "\n"
        out += f"Name: {self.name}"
        out += f"{nl}Description: {self.description}"
        if self.possible_values:
            out += f"{nl}Options: " + " | ".join(self.possible_values)
        out += f"{nl}Default: {self.default}"
        return out




parameters_general = {
    'trajectory': ParameterInfo(
        name = 'trajectory',
        description = "The trajectory type used for simulation",
        default = "Lorentz",
        possible_values=trajectories,
    ),
    'background': ParameterInfo(
        name = 'background',
        description = "The background type used for simulation",
        default = "Dipole",
        possible_values=backgrounds,
    ),
    'diffusion': ParameterInfo(
        name = 'diffusion',
        description = "The diffusion type used for simulation",
        default = "Uniform",
        possible_values=diffusions,
    ),
    # ----- General ----- #
    'build_mode': ParameterInfo(
        name = 'build_mode',
        description = "The build mode used for simulation",
        default = "release",
        possible_values=["release", "debug"],
        parameter_type="BuildMode",
    ),
    'specie': ParameterInfo(
        name = 'specie',
        description = "The specie (particle species) used for simulation",
        default = "proton_core",
        possible_values=["proton_core", "electron_core"], # todo
    ),
    'num_trajectories': ParameterInfo(
        name = 'num_trajectories',
        description = "The number of trajectories solved during simulation",
        default = 1,
    ),
    'batch_size': ParameterInfo(
        name = 'batch_size',
        description = "The number of trajectories in a single batch during simulation",
        default = 1,
    ),
}


# entries must match list-wise with those in background.config.hh
parameters_background = {
    'derivative_method': ParameterInfo(
        name = 'derivative_method',
        description = "The method used to evaluate derivatives of spatially located field quantities",
        default = "numeric",
        possible_values=["analytic", "numeric"],
        parameter_type="BackgroundOptions::DerivativeMethod",
    ),
    'num_numeric_grad_evals': ParameterInfo(
        name = 'num_numeric_grad_evals',
        description = "The number of derivative evaluations applied in the derivative-averaging method",
        default = 1,
    ),
    'incr_dmax_ratio': ParameterInfo(
        name = 'incr_dmax_ratio',
        description = "What fraction of _dmax to use to calculate the field increment",
        default = 0.0001,
        parameter_type = float,
    ),
    'server_interpolation_order': ParameterInfo(
        name = 'server_interpolation_order',
        description = "server_interpolation_order",
        default = 1,
    ),
    'smooth_discontinuity_order': ParameterInfo(
        name = 'smooth_discontinuity_order',
        description = "Parameter controlling smoothness of discontinuity/shock",
        default = 4,
        possible_values=["0: not continuous", "1: differentiable", "2: twice differentiable", "3: thrice differentiable", "4, 5, >5: smooth"],
    ),
    'server_num_ghost_cells': ParameterInfo(
        name = 'server_num_ghost_cells',
        description = "number of ghost cells (server parameter)",
        default = 2,
    ),
    # [test_parker_spiral.cc] Make sure that "SOLARWIND_CURRENT_SHEET" and SOLARWIND_POLAR_CORRECTION are (#)defined as 0 in src/background_solarwind.hh.
    'solarwind_current_sheet': ParameterInfo(
        name = 'solarwind_current_sheet',
        description = "Heliospheric current sheet",
        default = 'disabled',
        possible_values=["disabled", "flat", "wavy_static: wavy (Jokipii-Thomas 1981) and static", "wavy_time_dependent: wavy and time-dependent"],
        parameter_type="BackgroundOptions::CurrentSheet",
    ),
    'solarwind_sectored_region': ParameterInfo(
        name = 'solarwind_sectored_region',
        description = "Magnetic topology region",
        default = 'nowhere',
        possible_values=["nowhere", "HCS"],
        parameter_type="BackgroundOptions::SectoredRegion",
    ),
    'solarwind_polar_correction': ParameterInfo(
        name = 'solarwind_polar_correction',
        description = "Correction to Parker Spiral, mainly for polar regions",
        default = 'none',
        possible_values=["none", "Smith_Bieber: Smith-Bieber 1991", "Zurbuchen_etal: Zurbuchen et al. 1997", "Schwadron_McComas: Schwadron-McComas 2003"],
        parameter_type="BackgroundOptions::PolarCorrection",
    ),
    'solarwind_speed_latitude_profile': ParameterInfo(
        name = 'solarwind_speed_latitude_profile',
        description = "Latitudinal profile for bulk speed",
        default = 'constant',
        possible_values=["constant", "linear_step", "smooth_step"],
        parameter_type="BackgroundOptions::SpeedLatitudeProfile",
    ),
    'solarwind_termshock_speed_exponent': ParameterInfo(
        name = 'solarwind_termshock_speed_exponent',
        description = "Integer exponent of decrease of solar wind speed beyond the termination shock",
        default = 'square',
        possible_values=["zero", "one", "square", "cube"],
        parameter_type="BackgroundOptions::TermShockSpeedExponent",
    ),
    'mod_type': ParameterInfo(
        name = 'mod_type',
        description = "What function to use within 'get_ampfactor'",
        default = 'scaled',
        possible_values=["none", "zero", "constant", "scaled"],
        parameter_type="BackgroundOptions::ModType",
    ),
    'mod_rpos': ParameterInfo(
        name = 'mod_rpos',
        description = "Whether to scale relative to s=0 or s=+inf",
        default = 'scale_rel_zero',
        possible_values=["scale_rel_zero", "scale_rel_inf"],
        parameter_type="BackgroundOptions::ModRPos",
    ),
}


# entries must match list-wise with those in trajectory.config.hh
parameters_trajectory = {
    'TrajectoryFields': ParameterInfo(
        name = 'TrajectoryFields',
        description = "The fields needed to advance the trajectory during simulation.",
        default = "Fields<FConfig<>>",
    ),
    'RecordCoordinates': ParameterInfo(
        name = 'RecordCoordinates',
        description = "The coordinates (and format spec) used to record the trajectory progress.",
        default = "Fields<FConfig<>, Pos_t, Time_t>",
    ),
    'time_flow': ParameterInfo(
        name = 'time_flow',
        description = "The time flow direction used for simulation",
        default = "forward",
        possible_values=["forward", "backward"],
        parameter_type="TrajectoryOptions::TimeFlow",
    ),
    'rk_integrator': ParameterInfo(
        name = 'rk_integrator',
        description = "The discretization scheme used for simulation, a Runge-Kutta scheme",
        default = "DormandPrince_54E",
        # todo extend list
        possible_values=["Euler1E", "Euler1I", "Midpoint_2I", "Midpoint_2E", "RungeKutta_4E", "Kutta38_4E", "RungeKuttaFehlberg_54E", "CashKarp_54E", "DormandPrince_54E", "RungeKuttaFehlberg65E", "RungeKuttaFehlberg76E", "RungeKuttaFehlberg_87E"],
        parameter_type="RKIntegrator",
    ),
    'record_mag_extrema': ParameterInfo(
        name = 'record_mag_extrema',
        description = "vcx",
        default = "vcx",
    ),
    'record_trajectory': ParameterInfo(
        name = 'record_trajectory',
        description = "vcx",
        default = "vcx",
    ),
    'record_trajectory_segment_presize': ParameterInfo(
        name = 'record_trajectory_segment_presize',
        description = "vcx",
        default = "vcx",
    ),
    'trajectory_adv_safety_level': ParameterInfo(
        name = 'trajectory_adv_safety_level',
        description = "vcx",
        default = "vcx",
    ),
    'max_trajectory_steps': ParameterInfo(
        name = 'max_trajectory_steps',
        description = "vcx",
        default = "vcx",
    ),
    'max_time_adaptations': ParameterInfo(
        name = 'max_time_adaptations',
        description = "vcx",
        default = "vcx",
    ),
    'n_max_calls': ParameterInfo(
        name = 'n_max_calls',
        description = "vcx",
        default = "vcx",
    ),
    'cfl_advection': ParameterInfo(
        name = 'cfl_advection',
        description = "vcx",
        default = "vcx",
    ),
    'cfl_diffusion': ParameterInfo(
        name = 'cfl_diffusion',
        description = "vcx",
        default = "vcx",
    ),
    'cfl_acceleration': ParameterInfo(
        name = 'cfl_acceleration',
        description = "vcx",
        default = "vcx",
    ),
    'cfl_pitchangle': ParameterInfo(
        name = 'cfl_pitchangle',
        description = "vcx",
        default = "vcx",
    ),
    'drift_safety': ParameterInfo(
        name = 'drift_safety',
        description = "vcx",
        default = "vcx",
    ),
    'mirror_threshold': ParameterInfo(
        name = 'mirror_threshold',
        description = "vcx",
        default = "vcx",
    ),
    'pperp_method': ParameterInfo(
        name = 'pperp_method',
        description = "vcx",
        default = "vcx",
    ),
    'use_B_drifts': ParameterInfo(
        name = 'use_B_drifts',
        description = "vcx",
        default = "vcx",
    ),
    'stochastic_method': ParameterInfo(
        name = 'stochastic_method',
        description = "vcx",
        default = "vcx",
    ),
    'stochastic_method_mu': ParameterInfo(
        name = 'stochastic_method_mu',
        description = "vcx",
        default = "vcx",
    ),
    'stochastic_method_perp': ParameterInfo(
        name = 'stochastic_method_perp',
        description = "vcx",
        default = "vcx",
    ),
    'split_scatt': ParameterInfo(
        name = 'split_scatt',
        description = "vcx",
        default = "vcx",
    ),
    'split_scatt_fraction': ParameterInfo(
        name = 'split_scatt_fraction',
        description = "vcx",
        default = "vcx",
    ),
    'const_dmumax': ParameterInfo(
        name = 'const_dmumax',
        description = "vcx",
        default = "vcx",
    ),
    'steps_per_orbit': ParameterInfo(
        name = 'steps_per_orbit',
        description = "vcx",
        default = "vcx",
    ),
    'divk_method': ParameterInfo(
        name = 'divk_method',
        description = "vcx",
        default = "vcx",
    ),
    'dlnp_max': ParameterInfo(
        name = 'dlnp_max',
        description = "vcx",
        default = "vcx",
    ),
}

# entries must match list-wise with those in diffusion.config.hh
parameters_diffusion = {
    'use_qlt_scatt': ParameterInfo(
        name = 'use_qlt_scatt',
        description = "vcx",
        default = False,
    ),
    'D0': ParameterInfo(
        name = 'D0',
        description = "vcx",
        default = "vcx",
    ),
    'Dperp': ParameterInfo(
        name = 'Dperp',
        description = "vcx",
        default = "vcx",
    ),
    'Dpara': ParameterInfo(
        name = 'Dpara',
        description = "vcx",
        default = "vcx",
    ),
    'A2A': ParameterInfo(
        name = 'A2A',
        description = "vcx",
        default = "vcx",
    ),
    'l_max': ParameterInfo(
        name = 'l_max',
        description = "vcx",
        default = "vcx",
    ),
    'k_min': ParameterInfo(
        name = 'k_min',
        description = "vcx",
        default = "vcx",
    ),
    'ps_index': ParameterInfo(
        name = 'ps_index',
        description = "vcx",
        default = "vcx",
    ),
    'ps_minus': ParameterInfo(
        name = 'ps_minus',
        description = "vcx",
        default = "vcx",
    ),
    'A2T': ParameterInfo(
        name = 'A2T',
        description = "vcx",
        default = "vcx",
    ),
    'A2L': ParameterInfo(
        name = 'A2L',
        description = "vcx",
        default = "vcx",
    ),
    'ps_plus': ParameterInfo(
        name = 'ps_plus',
        description = "vcx",
        default = "vcx",
    ),
    'l_max_HP': ParameterInfo(
        name = 'l_max_HP',
        description = "vcx",
        default = "vcx",
    ),
    'dl_max': ParameterInfo(
        name = 'dl_max',
        description = "vcx",
        default = "vcx",
    ),
    'nose_z_nose': ParameterInfo(
        name = 'nose_z_nose',
        description = "vcx",
        default = "vcx",
    ),
    'nose_z_sheath': ParameterInfo(
        name = 'nose_z_sheath',
        description = "vcx",
        default = "vcx",
    ),
    'nose_z_dz': ParameterInfo(
        name = 'nose_z_dz',
        description = "vcx",
        default = "vcx",
    ),
    'kappa0': ParameterInfo(
        name = 'kappa0',
        description = "vcx",
        default = "vcx",
    ),
    'kappa_ratio': ParameterInfo(
        name = 'kappa_ratio',
        description = "vcx",
        default = "vcx",
    ),
    'U0': ParameterInfo(
        name = 'U0',
        description = "vcx",
        default = "vcx",
    ),
    'p0': ParameterInfo(
        name = 'p0',
        description = "vcx",
        default = "vcx",
    ),
    'T0': ParameterInfo(
        name = 'T0',
        description = "vcx",
        default = "vcx",
    ),
    'r0': ParameterInfo(
        name = 'r0',
        description = "vcx",
        default = "vcx",
    ),
    'lam0': ParameterInfo(
        name = 'lam0',
        description = "vcx",
        default = "vcx",
    ),
    'R0': ParameterInfo(
        name = 'R0',
        description = "vcx",
        default = "vcx",
    ),
    'B0': ParameterInfo(
        name = 'B0',
        description = "vcx",
        default = "vcx",
    ),
    'pow_law_U': ParameterInfo(
        name = 'pow_law_U',
        description = "vcx",
        default = "vcx",
    ),
    'pow_law_p': ParameterInfo(
        name = 'pow_law_p',
        description = "vcx",
        default = "vcx",
    ),
    'pow_law_T': ParameterInfo(
        name = 'pow_law_T',
        description = "vcx",
        default = "vcx",
    ),
    'pow_law_r': ParameterInfo(
        name = 'pow_law_r',
        description = "vcx",
        default = "vcx",
    ),
    'pow_law_R': ParameterInfo(
        name = 'pow_law_R',
        description = "vcx",
        default = "vcx",
    ),
    'pow_law_B': ParameterInfo(
        name = 'pow_law_B',
        description = "vcx",
        default = "vcx",
    ),
    'LISM_idx': ParameterInfo(
        name = 'LISM_idx',
        description = "vcx",
        default = "vcx",
    ),
    'LISM_ind': ParameterInfo(
        name = 'LISM_ind',
        description = "vcx",
        default = "vcx",
    ),
    'kappa_ratio_inner': ParameterInfo(
        name = 'kappa_ratio_inner',
        description = "vcx",
        default = "vcx",
    ),
    'kappa_ratio_outer': ParameterInfo(
        name = 'kappa_ratio_outer',
        description = "vcx",
        default = "vcx",
    ),
    'lam_inner': ParameterInfo(
        name = 'lam_inner',
        description = "vcx",
        default = "vcx",
    ),
    'lam_outer': ParameterInfo(
        name = 'lam_outer',
        description = "vcx",
        default = "vcx",
    ),
    'kappa_inner': ParameterInfo(
        name = 'kappa_inner',
        description = "vcx",
        default = "vcx",
    ),
    'kappa_outer': ParameterInfo(
        name = 'kappa_outer',
        description = "vcx",
        default = "vcx",
    ),
    'lam_para': ParameterInfo(
        name = 'lam_para',
        description = "vcx",
        default = "vcx",
    ),
    'lam_perp': ParameterInfo(
        name = 'lam_perp',
        description = "vcx",
        default = "vcx",
    ),
    'kappa_ratio_red': ParameterInfo(
        name = 'kappa_ratio_red',
        description = "vcx",
        default = "vcx",
    ),
    'radial_limit_perp_red': ParameterInfo(
        name = 'radial_limit_perp_red',
        description = "vcx",
        default = "vcx",
    ),
    'solar_cycle_idx': ParameterInfo(
        name = 'solar_cycle_idx',
        description = "vcx",
        default = "vcx",
    ),
    'solar_cycle_effect': ParameterInfo(
        name = 'solar_cycle_effect',
        description = "vcx",
        default = "vcx",
    ),
    'Bmix_idx': ParameterInfo(
        name = 'Bmix_idx',
        description = "vcx",
        default = "vcx",
    ),
    'Bmix_ind': ParameterInfo(
        name = 'Bmix_ind',
        description = "vcx",
        default = "vcx",
    ),
}


parameters = parameters_general | parameters_background | parameters_trajectory | parameters_diffusion


