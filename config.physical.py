# File config.physical.py Created by Lucius Schoenbaum October 20, 2025

"""
The following dict of dicts define a reasonable set of defaults
for particular trajectory/background/diffusion types.
These defaults are intended to be more useful
than the fallthrough defaults. This information makes this script
the location in the code where these default values are maintained.

To aid both the maintainer/developer and user,
only the *relevant* parameters are included
in each list. The maintainer/developer can therefore
focus on a particular class whose default values are being intelligently
chosen based on physics and numerics. Meanwhile, the
lists and dicts in `config.py` are more software-oriented,
and less physics-oriented.
"""


physical_defaults = {
    'background':{
        'CylindricalObstacle': {
            'derivative_method': 'analytic',
        },
        'Dipole': {
        },
        'Discontinuity': {
        },
        'MagnetizedCylinder': {
        },
        'MagnetizedSphere': {
        },
        'Server': {
        },
        'ServerBATL': {
        },
        'ServerCartesian': {
        },
        'Shock': {
        },
        'SmoothDiscontinuity': {
        },
        'SmoothShock': {
        },
        'SolarWind': {
            'solarwind_current_sheet': 'disabled',
            'solarwind_sectored_region': 'nowhere',
            'solarwind_polar_correction': 'none',
            'solarwind_speed_latitude_profile': 'constant',
            'solarwind_termshock_speed_exponent': 'square',
        },
        'SolarWindTermShock': {
        },
        'SphericalObstacle': {
        },
        'Uniform': {
        },
        'VLISMBochum': {
            'mod_type': 'scaled',
            'mod_rpos': 'scale_rel_zero',
        },
        'Waves': {
        },
    },
    'trajectory': {
        'Fieldline': {
            # nothing at this time.
        },
        'Focused': {
            'pperp_method': "scheme",
            'use_B_drifts': "none",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'Guiding': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'GuidingDiff': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'GuidingScatt': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'split_scatt': False,
            'split_scatt_fraction': 0.5,
            'const_dmumax': "constant_dthetamax",
            'stochastic_method': "Euler",
            'cfl_pitchangle': 0.5,
        },
        'GuidingDiffScatt': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'split_scatt': False,
            'split_scatt_fraction': 0.5,
            'const_dmumax': "constant_dthetamax",
            'stochastic_method_mu': "Euler",
            'stochastic_method_perp': "Euler",
            'cfl_pitchangle': 0.5,
        },
        'Lorentz': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 100000,
            'cfl_advection': 0.1,
            'drift_safety': 0.5,
            'mirror_threshold': 300,
            'steps_per_orbit': 100,
        },
        'Parker': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 100000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'stochastic_method': "Euler",
            'use_B_drifts': "none",
            'divk_method': "direct",
            'cfl_diffusion': 0.5,
            'cfl_acceleration': 0.5,
            'dlnp_max': 0.01,
        },
    },
    'diffusion': {
        'IsotropicConstant': {
            'D0': 1234,
        },
        'QLTConstant': {
            'A2A': 1234,
            'l_max': 1234,
            'k_min': 1234,
            'ps_index': 1234,
            'ps_minus': 1234,
        },
        'WNLTConstant': {
            'use_qlt_scatt': False,
            'A2A': 1234,
            'l_max': 1234,
            'k_min': 1234,
            'ps_index': 1234,
            'ps_minus': 1234,
            'A2T': 1234,
            'A2L': 1234,
            'ps_plus': 1234,
        },
        'WNLTRampVLISM': {
            'A2A': 1234,
            'l_max': 1234,
            'k_min': 1234,
            'ps_index': 1234,
            'ps_minus': 1234,
            'l_max_HP': 1234,
            'dl_max': 1234,
            'nose_z_nose': 1234,
            'nose_z_sheath': 1234,
            'nose_dz': 1234,
            'kappa0': 1234,
            'kappa_ratio': 1234,
        },
        'ParaConstant': {
            'D0': 1234,
        },
        'PerpConstant': {
            'D0': 1234,
        },
        'FullConstant': {
            'Dperp': 1234,
            'Dpara': 1234,
        },
        'FlowMomentumPowerLaw': {
            'kappa0': 1234,
            'kappa_ratio': 1234,
            'U0': 1234,
            'p0': 1234,
            'pow_law_U': 1234,
            'pow_law_p': 1234,
        },
        'KineticEnergyRadialDistancePowerLaw': {
            'kappa0': 1234,
            'kappa_ratio': 1234,
            'T0': 1234,
            'r0': 1234,
            'pow_law_T': 1234,
            'pow_law_r': 1234,
        },
        'RigidityMagneticFieldPowerLaw': {
            'kappa_ratio': 1234,
            'lam0': 1234,
            'R0': 1234,
            'B0': 1234,
            'pow_law_R': 1234,
            'pow_law_B': 1234,
        },
        'StraussEtAl2013': {
            'R0': 1234,
            'B0': 1234,
            'LISM_idx': 1234,
            'LISM_ind': 1234,
            'kappa_ratio_inner': 1234,
            'kappa_ratio_outer': 1234,
            'lam_inner': 1234,
            'lam_outer': 1234,
        },
        'GuoEtAl2014': {
            'R0': 1234,
            'B0': 1234,
            'LISM_idx': 1234,
            'LISM_ind': 1234,
            'kappa_ratio_inner': 1234,
            'kappa_ratio_outer': 1234,
            'lam_inner': 1234,
            'lam_outer': 1234,
        },
        'PotgieterEtAl2015': {
            'R0': 1234,
            'B0': 1234,
            'LISM_idx': 1234,
            'LISM_ind': 1234,
            'kappa_ratio_inner': 1234,
            'kappa_ratio_outer': 1234,
            'kappa_inner': 1234,
            'kappa_outer': 1234,
        },
        'EmpiricalSOQLTandUNLT': {
            'R0': 1234,
            'B0': 1234,
            'lam_para': 1234,
            'lam_perp': 1234,
            'kappa_ratio_red': 1234,
            'radial_limit_perp_red': 1234,
            'solar_cycle_idx': 1234,
            'solar_cycle_effect': 1234,
            'Bmix_idx': 1234,
            'Bmix_ind': 1234,
        },
    },
}




