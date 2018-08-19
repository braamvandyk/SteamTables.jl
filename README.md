# SteamTables

A Julia implementation of the IAPWS-IF97 properties of water and steam.
Provides the Gibbs and Helmholtz free energies, enthalpy, entropy, Cp, Cv and sonic velocity.

Inputs are either P&T, P&h or P&s.

[![SteamTables](http://pkg.julialang.org/badges/SteamTables_0.6.svg)](http://pkg.julialang.org/?pkg=SteamTables)

Linux and macOS: [![Build Status](https://travis-ci.org/braamvandyk/SteamTables.jl.svg?branch=master)](https://travis-ci.org/braamvandyk/SteamTables.jl)

Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/braamvandyk/SteamTables.jl?branch=master&svg=true)](https://ci.appveyor.com/project/braamvandyk/SteamTables-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/braamvandyk/SteamTables.jl/badge.svg?branch=master)](https://coveralls.io/r/braamvandyk/SteamTables.jl?branch=master)
[![codecov.io](http://codecov.io/github/braamvandyk/SteamTables.jl/coverage.svg?branch=master)](http://codecov.io/github/braamvandyk/SteamTables.jl?branch=master)

## Exported functions

### Single input

  Psat(T) and Tsat(P) returns the saturation linearindices

### Two inputs

#### P and T

  SpecificG, SpecificF, SpecificV, SpecificU, SpecificS, SpecificH, SpecificCP, SpecificCV, SpeedOfSound

#### P and h

  SpecificG_Ph, SpecificF_Ph, SpecificV_Ph, SpecificU_Ph, SpecificS_Ph, SpecificH_Ph, SpecificCP_Ph, SpecificCV_Ph, SpeedOfSound_Ph

#### P and s

  SpecificG_Ps, SpecificF_Ps, SpecificV_Ps, SpecificU_Ps, SpecificS_Ps, SpecificH_Ps, SpecificCP_Ps, SpecificCV_Ps, SpeedOfSound_Ps

    SpecificG     [kJ/kg]   Specific Gibbs free energy

    SpecificF     [kJ/kg]   Specific Helmholtz free energy

    SpecificV     [m3/kg]   Specific volume

    SpecificU     [kJ/kg]   Specific internal energy

    SpecificS     [kJ/kgK]  Specific entropy

    SpecificH     [kJ/kg]   Specific enthalpy

    SpecificCp    [kJ/kgK]  Specific isobaric heat capacity

    SpecificCv    [kJ/kgK]  Specific isochoric heat capacity

    SpeedOfSound  [m/s]     Sonic velocity

    Temperatures in K, Pressures in MPa

## Exported constants

    R  = 0.461526   [kJ/kg/K] Universal gas constant

    Tc = 647.096    [K]       Critical temperature of water

    Pc = 22.064     [kg/m3]   Critical density of water

    T3 = 273.16     [K]       Triple point temperature of water

    P3 = 611.657E-6 [MPa]     Triple point pressure of water

    Mr = 18.01528   [kg/kmol] Molecular weight of water
  