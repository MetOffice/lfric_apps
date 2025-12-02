import pytest
from ..extract_checker import *

# Generate some inputs and outputs that will be used in unit tests below
list_of_strings = ["first string", "second string", "third string"]
extract_source_input = ["extract.location{diff}[um] = $um_sources",
 "extract.path-excl[um] = src/atmosphere/convection/comorph \\",
 "                        / # everything",
 "extract.path-incl[um] = src/atmosphere/AC_assimilation/iau_mod.F90 \\",
 "                        src/atmosphere/aerosols \\",
 "                        src/atmosphere/atmosphere_service \\"]

set_of_files = {
"science/src/shumlib/shum_byteswap/src/c_shum_byteswap.c",
"science/src/um/src/control/top_level/init_corner_pr.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/invert_mod.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/acssrw_mod.F90",
"science/src/jules/src/science/soil/vmc_from_head_mod.F90",
"science/src/casim/src/dust_hack.F90",
"science/src/shumlib/shum_number_tools/src/f_shum_is_denormal.F90",
"science/src/jules/src/science/vegetation/emerge.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/acscos_mod.F90",
"science/src/shumlib/shum_number_tools/src/f_shum_is_nan.F90",
"science/src/um/src/control/coupling/oasis3_geto2a.F90",
"science/src/jules/src/science/vegetation/sow.F90",
"science/src/shumlib/shum_data_conv/src/c_shum_data_conv.c",
"science/src/um/src/control/grids/calc_npmsl_redbl.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/acsno_mod.F90",
"science/src/jules/src/control/shared/jules_soil_ecosse_mod.F90",
"science/src/um/src/control/misc/address_check.F90",
"science/src/casim/src/ship_tracks.F90",
}
set_of_um_files = {
"science/src/um/src/control/top_level/init_corner_pr.F90",
"science/src/um/src/control/coupling/oasis3_geto2a.F90",
"science/src/um/src/control/grids/calc_npmsl_redbl.F90",
"science/src/um/src/control/misc/address_check.F90",
}
set_of_shumlib_files = {
"science/src/shumlib/shum_byteswap/src/c_shum_byteswap.c",
"science/src/shumlib/shum_number_tools/src/f_shum_is_denormal.F90",
"science/src/shumlib/shum_number_tools/src/f_shum_is_nan.F90",
"science/src/shumlib/shum_data_conv/src/c_shum_data_conv.c",
}
set_of_ukca_files = {
"science/src/ukca/src/science/photolysis/stratospheric/photolib/invert_mod.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/acssrw_mod.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/acscos_mod.F90",
"science/src/ukca/src/science/photolysis/stratospheric/photolib/acsno_mod.F90",
}
empty_set = set([])

# Test subsetting_by_pattern()
@pytest.mark.parametrize(
    ("Set_In", "pattern", "expected"),
    [
        (set_of_files, "ukca", set_of_ukca_files),
        (set_of_files, "um", set_of_um_files),
        (set_of_files, "shumlib", set_of_shumlib_files),
        (set_of_files, "", set_of_files),
        (set_of_files, "no-such-repos", empty_set),
    ]
    )
def test_subsetting_by_pattern(Set_In, pattern, expected):
    assert subsetting_by_pattern(Set_In, "science/src/"+pattern) == expected


# Test extract_paths()
@pytest.mark.parametrize(
    ("inlist", "out1", "out2"),
    [
        (["21:42:24 *Compile* algorithm/iterative_solvers/pressure_precon_alg_mod.f90",
          "21:32:59 *Analysing* science/src/um/src/control/ukca_interface/ukca_option_mod.F90",
          "src/jules/src/science/vegetation/woodprod.F90",
          "src/jules/src/science_cable/",
          "21:17:42 *Copying source ../gungho/source/kernel/transport/ffsl/cosmic_flux_mod.F90",
          "sending incremental file list",
          "./"],
         {'algorithm/iterative_solvers/pressure_precon_alg_mod.f90',
          'science/src/um/src/control/ukca_interface/ukca_option_mod.F90',
          'src/jules/src/science/vegetation/woodprod.F90',
          '../gungho/source/kernel/transport/ffsl/cosmic_flux_mod.F90'},
         {'src/jules/src/science_cable/',
          './'}
         )
    ]
)
def test_extract_paths(inlist, out1, out2):
    assert extract_paths(inlist) == (out1, out2)


# Test pattern_search()
@pytest.mark.parametrize(
    ("inlist", "pattern", "expected"),
    [
        (list_of_strings, "first", {"first string": 1}),
        (list_of_strings, "string",{"first string": 1, "second string": 1, "third string": 1}),
        (list_of_strings + list_of_strings, "second", {"second string": 2})
    ]
)
def test_pattern_search(inlist, pattern, expected):
    assert pattern_search(inlist, pattern) == expected


# Test pattern_search_start()
@pytest.mark.parametrize(
    ("inlist", "pattern", "expected"),
    [
        (list_of_strings, "first", {"first string": 1}),
        (list_of_strings, "string", {}),
        (list_of_strings + list_of_strings, "second", {"second string": 2})
    ]
)
def test_pattern_search_start(inlist, pattern, expected):
    assert pattern_search(inlist, pattern, "start") == expected


# Test things_in_A_but_not_B()
@pytest.mark.parametrize(
    ("set_A", "set_B", "expected"),
    [
        (set([1,2,3,4,5]), set([1,2,4,5]), set([3])),
        (set([1,2,3,4,5]), set([1,2,3,4,5,6,7]), set([])),
    ]
)
def test_things_in_A_but_not_B(set_A, set_B, expected):
    assert things_in_A_but_not_B(set_A, set_B) == expected


# Test things_only_in_A_OR_B()
@pytest.mark.parametrize(
    ("set_A", "set_B", "expected"),
    [
        (set([1,2,3,4,5]), set([1,2,4,5]), set([3])),
        (set([1,2,3,4,5]), set([1,2,3,4,5,6,7]), set([6,7])),
        (set([1,2,3,4,5,6,7]), set([1,2,3,4,5]), set([6,7])),
        (set([1,2,3]), set([3,4,5]), set([1,2,4,5])),
    ]
)
def test_things_only_in_A_OR_B(set_A, set_B, expected):
    assert things_only_in_A_OR_B(set_A, set_B) == expected


# Test things_only_in_A_AND_B()
@pytest.mark.parametrize(
    ("set_A", "set_B", "expected"),
    [
        (set([1,2,3,4,5]), set([1,2,4,5]), set([1,2,4,5])),
        (set([1,2,3,4,5]), set([5,6,7]), set([5])),
        (set([1,3,5,6,7]), set([1,2,3,4,5]), set([1,3,5])),
        (set([1,2,3]), set([4,5]), set([])),
    ]
    )
def test_things_only_in_A_AND_B(set_A, set_B, expected):
    assert things_only_in_A_AND_B(set_A, set_B) == expected

