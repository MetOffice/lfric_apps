#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This module contains the default Fab configuration class.
'''

from typing import cast, Optional

from fab.api import (BuildConfig, Category, Compiler, ContainFlags,
                     ToolRepository)

from site_specific.default.config import Config as DefaultSiteConfig


class Config(DefaultSiteConfig):
    '''
    This class is the default application-specific Configuration object
    for LFRic ATM. It inherits from the LFRic site-specific configuration,
    and modifies compiler flags for certain files (e.g. UM physics).
    '''

    @staticmethod
    def get_fortran_compiler(name: str) -> Optional[Compiler]:
        """
        A small method to return a Fortran compiler with a given
        name if it is available. This function will test for the
        given compiler name, but also for mpif90 wrapper. This is
        required since in certain environments e.g. gfortran might
        not be in PATH, but the mpif90 wrapper around gfortran
        is. In this case, the Fab mpif90-gfortran compiler wrapper
        is returned.

        :param name: name of the compiler.
        :returns: the compiler if it is available, or None otherwise.
        """
        tool_repo = ToolRepository()
        compiler = tool_repo.get_tool(Category.FORTRAN_COMPILER, name)

        if not compiler.is_available:
            try:
                # There might not be a mpif90 tool, e.g.
                # for Cray compilers
                compiler = tool_repo.get_tool(Category.FORTRAN_COMPILER,
                                              f"mpif90-{name}")
            except KeyError:
                return None
            if not compiler.is_available:
                return None

        compiler = cast(Compiler, compiler)
        return compiler

    @staticmethod
    def set_um_physics_flags(compiler: Compiler,
                             flags: list[str]) -> None:
        """Set specific flags for UM physics files (e.g. typically
        at least 8-byte default reals, ...). This is a convenience
        function to keep the list of directories in only one place.

        :param compiler: the compiler for which to set the flags.
        :Param flags: the list of flags to set.
        """
        for pattern in ["/AC_assimilation/", "/aerosols/",
                        "/atmosphere_service/", "/boundary_layer/",
                        "/carbon/", "/convection/", "/legacy/",
                        "/diffusion_and_filtering/", "/dynamics/",
                        "/dynamics_advection/", "/electric/",
                        "/free_tracers/", "/gravity_wave_drag/",
                        "/idealised/", "/large_scale_cloud/",
                        "/large_scale_precipitation/",
                        "/physics_diagnostics/", "/radiation_control/",
                        "/stochastic_physics/", "/tracer_advection/",
                        '/casim/', "/jules/", "/socrates/", "/ukca/"]:
            compiler.add_flags(ContainFlags(pattern, flags), "base")

    def setup_cray(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Cray compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        super().setup_cray(build_config)
        ftn = self.get_fortran_compiler("crayftn-ftn")
        if not ftn:
            return

        self.set_um_physics_flags(ftn, ["-s", "real64"])

        # Fast-debug
        ftn.add_flags(ContainFlags("/parcel_ascent_5a.", "-hvector0"), "base")

    def setup_gnu(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Gnu compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        super().setup_gnu(build_config)
        gfortran = self.get_fortran_compiler("gfortran")
        if not gfortran:
            return

        um_physics = ["-fdefault-real-8", "-Wno-error=conversion"]
        self.set_um_physics_flags(gfortran, um_physics)

        # Buggy source code:
        # bl_ : bl_type_ind(map_bl(1,i)+0) = bl_type_1(i,1)  ...
        # iau_: iau_ts_end = iau_ts_start + iau_ts_num - 1.0_r_def
        #       timestep_index = nint( iau_time / dt ) + 1.0_i_def
        for pattern in ["/bl_exp_kernel_mod.f90", "/iau_time_control_mod.f90"]:
            gfortran.add_flags(ContainFlags(pattern, "-Wno-error=conversion"),
                               "base")

    def setup_intel_classic(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Intel classic compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        super().setup_intel_classic(build_config)
        ifort = self.get_fortran_compiler("ifort")
        if not ifort:
            return

        self.set_um_physics_flags(ifort, ["-r8"])

        # Some SOCRATES functions do not currently declare interfaces.
        # Flag was introduced in Intel Fortran v19.1.0 according to
        # Intel release notes.
        if (19, 1) <= ifort.get_version():
            no_externals = ["-warn", "noexternals"]
            for pattern in ["socrates/src/radiance_core",
                            "socrates/src/interface_core"]:
                ifort.add_flags(ContainFlags(pattern, no_externals),
                                "base")

    def setup_intel_llvm(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Intel LLVM compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
        super().setup_intel_llvm(build_config)
        ifx = self.get_fortran_compiler("ifx")
        if not ifx:
            return

        self.set_um_physics_flags(ifx, ["-r8"])

        # Some SOCRATES functions do not currently declare interfaces.
        # Flag was introduced in Intel Fortran v19.1.0 according to
        # Intel release notes.
        if (19, 1) <= ifx.get_version():
            no_externals = ["-warn", "noexternals"]
            for pattern in ["socrates/src/radiance_core",
                            "socrates/src/interface_core"]:
                ifx.add_flags(ContainFlags(pattern, no_externals),
                              "base")

        # Disabling OpenMP due to a compiler bug for intel-compiler newer
        # than 2020.3.304 - see Ticket 3853)
        if (2020, 3, 304) < ifx.get_version():
            ifx.add_flags(ContainFlags(pattern, "-qno-openmp"),
                          "base")

    def setup_nvidia(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Nvidia compiler and linker flags.

        :param build_config: the Fab build configuration instance
        '''
