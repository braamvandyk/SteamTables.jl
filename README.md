# SteamTables

A Julia implementation of the IAPWS-IF97 properties of water and steam.
Provides the Gibbs and Helmholtz free energies, enthalpy, entropy, Cp, Cv and sonic velocity.

Optional use of physical units via Unitful.jl

Inputs are either P&T, P&h or P&s.

[![Build status (Github Actions)](https://github.com/sylvaticus/MyAwesomePackage.jl/workflows/CI/badge.svg)](https://github.com/sylvaticus/MyAwesomePackage.jl/actions)
[![codecov.io](http://codecov.io/github/sylvaticus/MyAwesomePackage.jl/coverage.svg?branch=main)](http://codecov.io/github/sylvaticus/MyAwesomePackage.jl?branch=main)

## Exported functions

### Single input

Name         |Units  |Properties returned
-------------|-------|-------------------
Psat(T)      |MPa    |Saturation pressure
Tsat(P)      |K      |Saturation temperature
SatDensL(T)  |kg/m3  |Saturated liquid density
SatDensV(T)  |kg/m3  |Saturated vapour density
SatHL(T)     |J/kg   |Saturated liquid enthalpy
SatHV(T)     |J/kg   |Saturated vapour enthalpy
SatSL(T)     |J/kgK  |Saturated liquid entropy
SatSV(T)     |J/kgK  |Saturated vapour entropy
DeltaHvap    |J/kg   |Latent heat of vaporisation  

### Two inputs

#### P and T

    SpecificG, SpecificF, SpecificV, SpecificU, SpecificS, SpecificH, SpecificCP, SpecificCV, SpeedOfSound

#### P and h

    SpecificG_Ph, SpecificF_Ph, SpecificV_Ph, SpecificU_Ph, SpecificS_Ph, SpecificH_Ph, SpecificCP_Ph, SpecificCV_Ph, SpeedOfSound_Ph, Quality_Ph, Temperature_Ph

#### P and s

    SpecificG_Ps, SpecificF_Ps, SpecificV_Ps, SpecificU_Ps, SpecificS_Ps, SpecificH_Ps, SpecificCP_Ps, SpecificCV_Ps, SpeedOfSound_Ps, Quality_Ps, Temperature_Ps

#### T and h
    Quality_Th

#### T and s
    Quality_Ts


Name         |Units  |Properties returned
-------------|-------|-------------------
SpecificG    |kJ/kg  |Specific Gibbs free energy
SpecificF    |kJ/kg  |Specific Helmholtz free energy
SpecificV    |m3/kg  |Specific volume
SpecificU    |kJ/kg  |Specific internal energy
SpecificS    |kJ/kgK |Specific entropy
SpecificH    |kJ/kg  |Specific enthalpy
SpecificCp   |kJ/kgK |Specific isobaric heat capacity
SpecificCv   |kJ/kgK |Specific isochoric heat capacity
SpeedOfSound |m/s    |Sonic velocity
Quality      |       |Vapour quality  

Temperatures in K, Pressures in MPa


## Exported constants

Name |Value      |Units   |Physical Constant
-----|-----------|--------|-------------------
R    |0.461526   |kJ/kg/K |Universal gas constant
Tc   |647.096    |K       |Critical temperature of water
Pc   |22.064     |kg/m3   |Critical density of water
T3   |273.16     |K       |Triple point temperature of water
P3   |611.657E-6 |MPa     |Triple point pressure of water
Mr   |18.01528   |kg/kmol |Molecular weight of water

## Use of physical units

If unitless values are passed, the physical units listed above are assumed and the returned values will also be in these units, e.g. Kelvin for temperature, MPa for pressure.

Optionally, however, you may pass Unitful.jl Quantity values, e.g. 100u"°C". The returned values will still be in the specified units, but may be easily converted via Unitful's `uconvert()`

### Example
```
julia> using SteamTables, Unitful

julia> Psat(400)
0.24575318630408327

julia> Psat(400u"K")
0.24575318630408327 MPa

julia> Psat(300u"°C")
8.587708329557278 MPa

julia> uconvert(u"psi", Psat(400u"K"))
35.643486181534875 psi
```

For more details on the use of Unitful, see the package documentation.
