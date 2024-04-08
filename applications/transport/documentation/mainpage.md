@mainpage Overview

@section transport_miniapp Transport Miniapp

This miniapp forms part of the LFRic and Gungho project. It is a miniapp that allows transport-only tests with any of the transport methods to be run without having to use the full model.

A specified wind field is used to transport all fields. A number of fields are transported:
1. A W3 field in a conservative (flux) form equation (equivalent to density in the full model)
2. A W3 field in a non-conservative (advective) form equation (equivalent to the wind components in the full model)
3. A Wtheta field in a non-conservative (advective) form equation (equivalent to the potential temperature in the full model)
4. A Wtheta field in a consevative (flux) form equation as a mixing ratio (equivalent to moisture in the full model)

@section lfric_gungho_sec LFRic and GungHo

The documentation related to LFRic and GungHo projects is hosted at the <a href="https://code.metoffice.gov.uk/trac/home">Met Office Science Repository Service (MOSRS)</a>.

- <a href="https://code.metoffice.gov.uk/trac/lfric/wiki">LFRic Project Space</a>,
- <a href="https://code.metoffice.gov.uk/trac/lfric/wiki/GungHo">GungHo Project Space</a>.

More information about checking out and running the code can be found at <a href="https://code.metoffice.gov.uk/trac/lfric/wiki/LFRicTechnical/QuickStart">LFRic QuickStart page</a>.

@section psyclone_in_lfric_sec PSyclone in LFRic

The LFRic code uses <a href="https://code.metoffice.gov.uk/trac/lfric/wiki/PsycloneTool">PSyclone tool</a>, which is hosted at the
<a href="https://github.com/stfc/PSyclone">PSyclone GitHub repository</a>.

One of the PSyclone features used very explicitly in the LFRic code are <a href="https://github.com/stfc/PSyclone/blob/master/doc/built_ins.rst">Built-ins</a>:
operations which can be specified within an invoke call in the algorithm layer but do not require an associated kernel to be implemented as they are provided
directly by the infrastructure.

For more up-to-date information about the <b>LFRic-specific Built-ins</b> functionality (e.g. names, argument order) please refer to the
<a href="https://github.com/stfc/PSyclone/blob/master/doc/dynamo0p3.rst#built-ins">dynamo 0.3 API Built-ins documentation</a>.


