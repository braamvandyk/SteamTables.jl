# SteamTables

A Julia implementation of the IAPWS-IF97 properties of water and steam. 
Provides the Gibbs and Helmholtz free energies, enthalpy, entropy Cp, Cv and sonic velocity.
Inputs are eith P&T, P&h or P&s.

## Exported functions:

### Single input:

  Psat(T) and Tsat(P) returns the saturation linearindices

### Two inputs:

P and T
  SpecificG,    SpecificF,     SpecificV,     SpecificU,      SpecificS,
  SpecificH,    SpecificCP,    SpecificCV,    SpeedOfSound

P and h
  SpecificG_Ph, SpecificF_Ph,  SpecificV_Ph,  SpecificU_Ph,   SpecificS_Ph,
  SpecificH_Ph, SpecificCP_Ph, SpecificCV_Ph, SpeedOfSound_Ph

P and s
  SpecificG_Ps, SpecificF_Ps,  SpecificV_Ps,  SpecificU_Ps,   SpecificS_Ps,
  SpecificH_Ps, SpecificCP_Ps, SpecificCV_Ps, SpeedOfSound_Ps
    

SpecificG    [kJ/kg]  Specific Gibbs free energy 
SpecificF    [kJ/kg]  Specific Helmholtz free energy 
SpecificV    [m3/kg]  Specific volume 
SpecificU    [kJ/kg]  Specific internal energy
SpecificS    [kJ/kgK] Specific entropy 
SpecificH    [kJ/kg]  Specific enthalpy 
SpecificCp   [kJ/kgK] Specific isobaric heat capacity 
SpecificCv   [kJ/kgK] Specific isochoric heat capacity 
SpeedOfSound [m/s]    Sonic velocity 

Temperatures in K, Pressures in MPa

## Exported constants
  R  = 0.461526   [kJ/kg/K] Universal gas constant
  Tc = 647.096    [K]       Critical temperature of water
  Pc = 22.064     [kg/m3]   Critical density of water
  T3 = 273.16     [K]       Triple point temperature of water
  P3 = 611.657E-6 [MPa]     Triple point pressure of water
  Mr = 18.01528   [kg/kmol] Molecular weight of water
