# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
PSyclone Script for applying OpenMP transformations to lw kernel mod
"""

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import (
    Call, Loop)
from transmute_psytrans.transmute_functions import (
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
    replace_n_threads,
    first_priv_red_init)

options = {
    "node-type-check": False,
    "ignore_dependencies_for": [
        "aer_lw_absorption",
        "aer_lw_asymmetry",
        "aer_lw_scattering",
        "aer_mix_ratio",
        "ch4",
        "cloud_drop_no_conc",
        "co",
        "co2",
        "conv_frozen_fraction",
        "conv_frozen_mmr",
        "conv_frozen_number",
        "conv_liquid_fraction",
        "conv_liquid_mmr",
        "cs",
        "d_mass",
        "frozen_fraction",
        "h2",
        "h2o",
        "hcn",
        "he",
        "layer_heat_capacity",
        "li",
        "liquid_fraction",
        "mcf",
        "mcl",
        "n2",
        "n2o",
        "n_ice",
        "na",
        "nh3",
        "o2",
        "o3",
        "potassium",
        "pressure_in_wth",
        "profile_list",
        "radiative_cloud_fraction",
        "radiative_conv_fraction",
        "rand_seed",
        "rb",
        "rho_in_wth",
        "sigma_mi",
        "sigma_ml",
        "so2",
        "sulphuric",
        "t_layer_boundaries",
        "temperature_in_wth",
        "tile_fraction",
        "tile_lw_albedo",
        "tile_temperature",
        "tio",
        "vo",
    ],
}

# Associate ll index with segmentation_ll tool type
Loop.set_loop_type_inference_rules({"segmentation_ll": {"variable": "ll"}})


def trans(psyir):
    """
    Entry point for OpenMP transformations for `lw_kernel_mod.F90`.
    """

    # Replace max_threads = 1
    replace_n_threads(psyir, "max_threads")

    # Enable pure for runes call to allow PSyclone to add OMP around it
    for call in psyir.walk(Call):
        if call.routine.name == "runes":
            call.routine.symbol.is_pure = True

    # Walk the loops of the psyir obj
    for loop in psyir.walk(Loop):
        # if the loop is of the set type above, segmentation_ll
        if loop.loop_type == "segmentation_ll":
            # first private workaround
            first_priv_red_init(loop, ["n_profile_list_seg", "seg_end"])
            # Apply the transformation
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop, options)
            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
