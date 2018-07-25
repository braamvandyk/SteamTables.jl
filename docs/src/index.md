#SteamTables.jl
 
This is a Julia implementation of the International Association for the Properties of Water and Steam Industrial Formulation (1997) models (IAPWS IF-97).
 
##Exported functions:

###Single input:

  Psat(T) and Tsat(P) returns the saturation linearindices

###Two inputs:

####P and T
  SpecificG,    SpecificF,     SpecificV,     SpecificU,      SpecificS,
  SpecificH,    SpecificCP,    SpecificCV,    SpeedOfSound

#####P and h
  SpecificG_Ph, SpecificF_Ph,  SpecificV_Ph,  SpecificU_Ph,   SpecificS_Ph,
  SpecificH_Ph, SpecificCP_Ph, SpecificCV_Ph, SpeedOfSound_Ph

####P and s
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

##Exported constants
  R  = 0.461526   [kJ/kg/K] Universal gas constant
  Tc = 647.096    [K]       Critical temperature of water
  Pc = 22.064     [kg/m3]   Critical density of water
  T3 = 273.16     [K]       Triple point temperature of water
  P3 = 611.657E-6 [MPa]     Triple point pressure of water
  Mr = 18.01528   [kg/kmol] Molecular weight of water


 
## Index
 
```@docs
SteamTable.B23(InputType::Symbol, InputValue)
SteamTable.B2bc(InputType::Symbol, InputValue)
SteamTable.Region1(Output::Symbol, P, T)
SteamTable.Region1_TPh(P, h)
SteamTable.Region1_TPs(P, s)
SteamTable.Region2(Output::Symbol, P, 
SteamTable.Region2meta(Output::Symbol, P, 
SteamTable.Region2a_TPh(P, h)
SteamTable.Region2b_TPh(P, h)
SteamTable.Region2c_TPh(P, h)
SteamTable.Region2a_TPs(P, s)
SteamTable.Region2b_TPs(P, s)
SteamTable.Region2c_TPs(P, s)
SteamTable.Region3_ρ(Output::Symbol, ρ, T) 
SteamTable.Region3(Output::Symbol, P, T) 
SteamTable.Region4(InputType::Symbol, InputValue)
SteamTable.Region5(Output::Symbol, P, T)
SteamTable.RegionID(P, T)::Symbol
SteamTable.RegionID_Ph(P, h)::Symbol
SteamTable.RegionID_Ps(P, s)::Symbol
SteamTable.Psat(T)
SteamTable.Tsat(P)
SteamTable.SpecificG(P, T)
SteamTable.SpecificG_Ph(P, h)
SteamTable.SpecificG_Ps(P, s)
SteamTable.SpecificF(P, T)
SteamTable.SpecificF_Ph(P, h)
SteamTable.SpecificF_Ps(P, s)
SteamTable.SpecificV(P, T)
SteamTable.SpecificV_Ph(P, h)
SteamTable.SpecificV_Ps(P, s)
SteamTable.SpecificU(P, T)
SteamTable.SpecificU_Ph(P, h)
SteamTable.SpecificU_Ps(P, s)
SteamTable.SpecificS(P, T)
SteamTable.SpecificS_Ph(P, h)
SteamTable.SpecificH(P, T)
SteamTable.SpecificH_Ps(P, s)
SteamTable.SpecificCP(P, T)
SteamTable.SpecificCP_Ph(P, h)
SteamTable.SpecificCP_Ps(P, s)
SteamTable.SpecificCV(P, T)
SteamTable.SpecificCV_Ph(P, h)
SteamTable.SpecificCV_Ps(P, s)
SteamTable.SpeedOfSound(P, T)
SteamTable.SpeedOfSound_Ph(P, h)
SteamTable.SpeedOfSound_Ps(P, s)
```