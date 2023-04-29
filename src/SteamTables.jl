"""
# SteamTables

IAPWS-97 Industrial Formulation Properties of water and steam.
Provides the Gibbs and Helmholtz free energies, enthalpy, entropy Cp, Cv and sonic velocity given inputs that are eith P&T, P&h or P&s.

## Exported functions:

### Single input:

  Psat(T) and Tsat(P) returns the saturation temperature or pressure
  SatDensL(T) and SatDensV(T) return the saturated densities
  SatHL(T) and SatHV(T) return the saturated enthalpies
  SatSL(T) and SatSV(T) return the saturated entropies
  DeltaHvap(T) returns the latent heat of vaporisation

### Two inputs:

P and T
  SpecificG,    SpecificF,     SpecificV,     SpecificU,      SpecificS,
  SpecificH,    SpecificCP,    SpecificCV,    SpeedOfSound

P and h
  SpecificG_Ph, SpecificF_Ph,  SpecificV_Ph,  SpecificU_Ph,   SpecificS_Ph,
  SpecificH_Ph, SpecificCP_Ph, SpecificCV_Ph, SpeedOfSound_Ph
  Quality_Ph

P and s
  SpecificG_Ps, SpecificF_Ps,  SpecificV_Ps,  SpecificU_Ps,   SpecificS_Ps,
  SpecificH_Ps, SpecificCP_Ps, SpecificCV_Ps, SpeedOfSound_Ps
  Quality_Ps

T and h
  Quality_Th

T and s
  Quality_Ts

SpecificG    [kJ/kg]  Specific Gibbs free energy
SpecificF    [kJ/kg]  Specific Helmholtz free energy
SpecificV    [m3/kg]  Specific volume
SpecificU    [kJ/kg]  Specific internal energy
SpecificS    [kJ/kgK] Specific entropy
SpecificH    [kJ/kg]  Specific enthalpy
SpecificCp   [kJ/kgK] Specific isobaric heat capacity
SpecificCv   [kJ/kgK] Specific isochoric heat capacity
SpeedOfSound [m/s]    Sonic velocity

By default temperatures in K, pressures in MPa

If units are associated with inputs via Unitful.jl, the return values will also
have default units associated with them. These can be converted to preferred units
in the calling function.

##Exported constants
  R  = 0.461526   [kJ/kg/K] Universal gas constant
  Tc = 647.096    [K]       Critical temperature of water
  Pc = 22.064     [MPa]     Critical pressure of water
  ρc = 322.       [kg/m3]   Critical density of water
  T3 = 273.16     [K]       Triple point temperature of water
  P3 = 611.657E-6 [MPa]     Triple point pressure of water
  Mr = 18.01528   [kg/kmol] Molecular weight of water

  Exported constants do not have associated units. Units can be easily added,
  e.g.

  Ru = R * 1.0u"kJ/kg/K"

"""
module SteamTables

using PrecompileTools

const R = 0.461526      #kJ/kg/K    Universal gas constant
const Tc = 647.096      #K          Critical temperature of water
const Pc = 22.064       #Pa         Critical pressure of water
const ρc = 322.         #kg/m3      Critical density of water
const T3 = 273.16       #K          Triple point temperature of water
const P3 = 611.657E-6   #MPa        Triple point pressure of water
const Mr = 18.01528     #kg/kmol    Molecular weight of water

export Psat, Tsat,
       SatDensL, SatDensV,
       SatHL, SatHV, SatSL, SatSV,
       DeltaHvap,
       Quality_Ph, Quality_Th, Quality_Ps, Quality_Ts,
       SpecificG,     SpecificF,     SpecificV,    SpecificU,    SpecificS,
       SpecificH,     SpecificCP,    SpecificCV,   SpeedOfSound,
       SpecificG_Ph,  SpecificF_Ph,  SpecificV_Ph, SpecificU_Ph, SpecificS_Ph,
       SpecificCP_Ph, SpecificCV_Ph, SpeedOfSound_Ph,
       SpecificG_Ps,  SpecificF_Ps,  SpecificV_Ps, SpecificU_Ps, SpecificH_Ps,
       SpecificCP_Ps, SpecificCV_Ps, SpeedOfSound_Ps,
       R, Tc, Pc, ρc, T3, P3, Mr

using Roots
using Unitful
using ForwardDiff

struct UnitsError <: Exception
    var
    message::String
end
Base.showerror(io::IO, e::UnitsError) = print(io, "UnitsError with ",e.var, ": ", e.message)


"""
    B23(InputType::Symbol, InputValue)

Returns the boundary between regions 2 and 3.
InputType is either :T or :P to indicate that InputValue is temperature [K] or pressure [MPa]. The complimentary value is returned.
"""
function B23(InputType::Symbol, InputValue)
    n = [0.348_051_856_289_69E3,
        -0.116_718_598_799_75E1,
         0.101_929_700_393_26E-2,
         0.572_544_598_627_46E3,
         0.139_188_397_788_70E2]

    if InputType == :T
        return n[1] + n[2] * InputValue + n[3] * InputValue^2
    elseif InputType == :P
        return n[4] + √((InputValue - n[5])/n[3])
    else
        throw(DomainError(InputType, "Unknown input. Expecting :T or :P"))
    end
end


"""
    B2bc

    Returns the boundary between regions 2b and 2c, approx s=5.85kJ/kgK
    InputType is either :h or :P to indicate that InputValue is enthalpy [kJ/kg]
    or pressure [MPa]. The complimentary value is returned.
    Valid from saturation line at 554.485K & 6.54670MPa to 1019.32K & 100MPa
"""
function B2bc(InputType::Symbol, InputValue)
    n = [0.905_842_785_147_23E3,
        -0.679_557_863_992_41,
         0.128_090_027_301_36E-3,
         0.265_265_719_084_28E4,
         0.452_575_789_059_48E1]

    if InputType == :h
        return n[1] + n[2] * InputValue + n[3] * InputValue^2
    elseif InputType == :P
        return n[4] + √((InputValue - n[5])/n[3])
    else
        throw(DomainError(InputType, "Unknown input. Expecting :P or :h."))
    end
end

"""
    Region1

    Returns all the property values in region 1.
    273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa
    Pressures in MPa and temperature in [K]
    The property to be returned is specified in Output.
        :SpecificG        kJ/kg
        :SpecificF        kJ/kg
        :SpecificV        m3/kg
        :SpecificU        kJ/kg
        :SpecificS        kJ/kgK
        :SpecificH        kJ/kg
        :SpecificCP       kJ/kgK
        :SpecificCV       kJ/kgK
        :SpeedOfSound     m/s
"""
function Region1(Output::Symbol, P, T)
    Pstar = 16.53   #MPa
    Tstar = 1386.0  #K

    n = [0.146_329_712_131_67,
        -0.845_481_871_691_14,
        -0.375_636_036_720_40E1,
         0.338_551_691_683_85E1,
        -0.957_919_633_878_72,
         0.157_720_385_132_28,
        -0.166_164_171_995_01E-1,
         0.812_146_299_835_68E-3,
         0.283_190_801_238_04E-3,
        -0.607_063_015_658_74E-3,
        -0.189_900_682_184_19E-1,
        -0.325_297_487_705_05E-1,
        -0.218_417_171_754_14E-1,
        -0.528_383_579_699_30E-4,
        -0.471_843_210_732_67E-3,
        -0.300_017_807_930_26E-3,
         0.476_613_939_069_87E-4,
        -0.441_418_453_308_46E-5,
        -0.726_949_962_975_94E-15,
        -0.316_796_448_450_54E-4,
        -0.282_707_979_853_12E-5,
        -0.852_051_281_201_03E-9,
        -0.224_252_819_080_00E-5,
        -0.651_712_228_956_01E-6,
        -0.143_417_299_379_24E-12,
        -0.405_169_968_601_17E-6,
        -0.127_343_017_416_41E-8,
        -0.174_248_712_306_34E-9,
        -0.687_621_312_955_31E-18,
         0.144_783_078_285_21E-19,
         0.263_357_816_627_95E-22,
        -0.119_476_226_400_71E-22,
         0.182_280_945_814_04E-23,
        -0.935_370_872_924_58E-25]

    I = [0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         1,
         1,
         1,
         1,
         1,
         1,
         2,
         2,
         2,
         2,
         2,
         3,
         3,
         3,
         4,
         4,
         4,
         5,
         8,
         8,
         21,
         23,
         29,
         30,
         31,
         32]

    J = [-2,
         -1,
         0,
         1,
         2,
         3,
         4,
         5,
         -9,
         -7,
         -1,
         0,
         1,
         3,
         -3,
         0,
         1,
         3,
         17,
         -4,
         0,
         6,
         -5,
         -2,
         10,
         -8,
         -11,
         -6,
         -29,
         -31,
         -38,
         -39,
         -40,
         -41]

    π = P / Pstar
    τ = Tstar / T
    ππ = 7.1 - π
    ττ = τ - 1.222

    γ    = sum([n[i]*(ππ^I[i])*(ττ^J[i]) for i=1:34])
    γ_π  = sum([-n[i]*I[i]*(ππ^(I[i]-1))*(ττ^J[i]) for i=1:34])
    γ_ππ = sum([n[i]*I[i]*(I[i]-1)*(ππ^(I[i]-2))*(ττ^J[i]) for i=1:34])
    γ_τ  = sum([n[i]*(ππ^I[i])*J[i]*(ττ^(J[i]-1)) for i=1:34])
    γ_ττ = sum([n[i]*(ππ^I[i])*J[i]*(J[i]-1)*(ττ^(J[i]-2)) for i=1:34])
    γ_πτ = sum([-n[i]*I[i]*(ππ^(I[i]-1))*J[i]*(ττ^(J[i]-1)) for i=1:34])

    if Output == :SpecificG            #kJ/kg
        return R*T*γ
    elseif Output == :SpecificF        #kJ/kg
        return R*T*(γ - π*γ_π/1000)
    elseif Output == :SpecificV        #m3/kg
        return R*T*π*γ_π/P/1000
    elseif Output == :SpecificU        #kJ/kg
        return R*T*(τ*γ_τ - π*γ_π)
    elseif Output == :SpecificS        #kJ/kgK
        return R*(τ*γ_τ - γ)
    elseif Output == :SpecificH        #kJ/kg
        return R*T*τ*γ_τ
    elseif Output == :SpecificCP       #kJ/kgK
        return -R*τ^2*γ_ττ
    elseif Output == :SpecificCV       #kJ/kgK
        return R*((-τ^2)*γ_ττ + (γ_π - τ*γ_πτ)^2/γ_ππ)
    elseif Output == :SpeedOfSound     #m/s
        return √(1000*R*T*(γ_π)^2/((γ_π - τ*γ_πτ)^2/τ^2/γ_ττ - γ_ππ))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end

"""
    Region1_TPh

    Returns T [K] from P[MPa] and h[kJ/kg] in Region 1.
    273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa
    The inconsistency with the main Region 1 model is ~25mK.
    If you need better accuracy (why???), start with this value and iterate.
"""
function Region1_TPh(P, h)
    hstar = 2500 #kJ/kg

    n = [-0.238_724_899_245_21E3,
          0.404_211_886_379_45E3,
          0.113_497_468_817_18E3,
         -0.584_576_160_480_39E1,
         -0.152_854_824_131_40E-3,
         -0.108_667_076_953_77E-5,
         -0.133_917_448_726_02E2,
          0.432_110_391_835_59E2,
         -0.540_100_671_705_06E2,
          0.305_358_922_039_16E2,
         -0.659_647_494_236_38E1,
          0.939_654_008_783_63E-2,
          0.115_736_475_053_40E-6,
         -0.258_586_412_820_73E-4,
         -0.406_443_630_847_99E-8,
          0.664_561_861_916_35E-7,
          0.806_707_341_030_27E-10,
         -0.934_777_712_139_47E-12,
          0.582_654_420_206_01E-14,
         -0.150_201_859_535_03E-16]

    I = [0,
         0,
         0,
         0,
         0,
         0,
         1,
         1,
         1,
         1,
         1,
         1,
         1,
         2,
         2,
         3,
         3,
         4,
         5,
         6]

    J = [0,
         1,
         2,
         6,
         22,
         32,
         0,
         1,
         2,
         3,
         4,
         10,
         32,
         10,
         32,
         10,
         32,
         32,
         32,
         32]

    η = h / hstar
    return sum([n[i]*P^I[i]*(η+1)^J[i] for i=1:20])
end

"""
    Region1_TPs

    Returns T [K] from P[MPa] and s[kJ/kgK] in Region 1.
    273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa
    The inconsistency with the main Region 1 model is ~25mK.
    If you need better accuracy (why???), start with this value and iterate.
"""
function Region1_TPs(P, s)
    n = [0.174_782_680_583_07E3,
         0.348_069_308_928_73E2,
         0.652_925_849_784_55E1,
         0.330_399_817_754_89,
        -0.192_813_829_231_96E-6,
        -0.249_091_972_445_73E-22,
        -0.261_076_364_893_32,
         0.225_929_659_815_86,
        -0.642_564_633_952_26E-1,
         0.788_762_892_705_26E-2,
         0.356_721_106_073_66E-9,
         0.173_324_969_948_95E-23,
         0.566_089_006_548_37E-3,
        -0.326_354_831_397_17E-3,
         0.447_782_866_906_32E-4,
        -0.513_221_569_085_07E-9,
        -0.425_226_570_422_07E-25,
         0.264_004_413_606_89E-12,
         0.781_246_004_597_23E-28,
        -0.307_321_999_036_68E-30]

    I = [0,
         0,
         0,
         0,
         0,
         0,
         1,
         1,
         1,
         1,
         1,
         1,
         2,
         2,
         2,
         2,
         2,
         3,
         3,
         4]

    J = [0,
         1,
         2,
         3,
         11,
         31,
         0,
         1,
         2,
         3,
         12,
         31,
         0,
         1,
         2,
         9,
         31,
         10,
         32,
         32]

    return sum([n[i]*P^I[i]*(s+2)^J[i] for i=1:20])
end

"""
    Region2

    Returns all the property values in region 2.
    273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
    623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
    863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa
    Accuracy in the metastable region is reasonable above 10MPa, only
    Pressures in MPa and temperature in [K]
    The property to be returned is specified in Output.
        :SpecificG        kJ/kg
        :SpecificF        kJ/kg
        :SpecificV        m3/kg
        :SpecificU        kJ/kg
        :SpecificS        kJ/kgK
        :SpecificH        kJ/kg
        :SpecificCP       kJ/kgK
        :SpecificCV       kJ/kgK
        :SpeedOfSound     m/s
"""
function Region2(Output::Symbol, P, T)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K

    no = [-0.969_276_865_002_17E1,
           0.100_866_559_680_18E2,
          -0.560_879_112_830_20E-2,
           0.714_527_380_814_55E-1,
          -0.407_104_982_239_28,
           0.142_408_191_714_44E1,
          -0.438_395_113_194_50E1,
          -0.284_086_324_607_72,
           0.212_684_637_533_07E-1]

    Jo = [0,
          1,
         -5,
         -4,
         -3,
         -2,
         -1,
          2,
          3]

    nr = [-0.177_317_424_732_13E-2,
          -0.178_348_622_923_58E-1,
          -0.459_960_136_963_65E-1,
          -0.575_812_590_834_32E-1,
          -0.503_252_787_279_30E-1,
          -0.330_326_416_702_03E-4,
          -0.189_489_875_163_15E-3,
          -0.393_927_772_433_55E-2,
          -0.437_972_956_505_73E-1,
          -0.266_745_479_140_87E-4,
           0.204_817_376_923_09E-7,
           0.438_706_672_844_35E-6,
          -0.322_776_772_385_70E-4,
          -0.150_339_245_421_48E-2,
          -0.406_682_535_626_49E-1,
          -0.788_473_095_593_67E-9,
           0.127_907_178_522_85E-7,
           0.482_253_727_185_07E-6,
           0.229_220_763_376_61E-5,
          -0.167_147_664_510_61E-10,
          -0.211_714_723_213_55E-2,
          -0.238_957_419_341_04E2,
          -0.590_595_643_242_70E-17,
          -0.126_218_088_991_01E-5,
          -0.389_468_424_357_39E-1,
           0.112_562_113_604_59E-10,
          -0.823_113_408_979_98E1,
           0.198_097_128_020_88E-7,
           0.104_069_652_101_74E-18,
          -0.102_347_470_959_29E-12,
          -0.100_181_793_795_11E-8,
          -0.808_829_086_469_85E-10,
           0.106_930_318_794_09,
          -0.336_622_505_741_71,
           0.891_858_453_554_21E-24,
           0.306_293_168_762_32E-12,
          -0.420_024_676_982_08E-5,
          -0.590_560_296_856_39E-25,
           0.378_269_476_134_57E-5,
          -0.127_686_089_346_81E-14,
           0.730_876_105_950_61E-28,
           0.554_147_153_507_78E-16,
          -0.943_697_072_412_10E-6]

    Ir = [1,
          1,
          1,
          1,
          1,
          2,
          2,
          2,
          2,
          2,
          3,
          3,
          3,
          3,
          3,
          4,
          4,
          4,
          5,
          6,
          6,
          6,
          7,
          7,
          7,
          8,
          8,
          9,
          10,
          10,
          10,
          16,
          16,
          18,
          20,
          20,
          20,
          21,
          22,
          23,
          24,
          24,
          24]

    Jr = [0,
          1,
          2,
          3,
          6,
          1,
          2,
          4,
          7,
          36,
          0,
          1,
          3,
          6,
          35,
          1,
          2,
          3,
          7,
          3,
          16,
          35,
          0,
          11,
          25,
          8,
          36,
          13,
          4,
          10,
          14,
          29,
          50,
          57,
          20,
          35,
          48,
          21,
          53,
          39,
          26,
          40,
          58]

    π = P / Pstar
    τ = Tstar / T

    γo    = log(π) + sum([no[i]*τ^Jo[i] for i=1:9])
    γo_π  = 1/π
    γo_ππ = -1/(π^2)
    γo_τ  = sum([no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9])
    γo_ττ = sum([no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9])
    γo_πτ = 0

    γr    = sum([nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:43])
    γr_π  = sum([nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:43])
    γr_ππ = sum([nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:43])
    γr_τ  = sum([nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:43])
    γr_ττ = sum([nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:43])
    γr_πτ = sum([nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:43])


    if Output == :SpecificG            #kJ/kg
        return R*T*(γo + γr)
    elseif Output == :SpecificF        #kJ/kg
        return R*T*(γo + γr - π*(γo_π + γr_π)/1000)
    elseif Output == :SpecificV        #m3/kg
        return R*T*π*(γo_π + γr_π)/P/1000
    elseif Output == :SpecificU        #kJ/kg
        return R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))
    elseif Output == :SpecificS        #kJ/kgK
        return R*(τ*(γo_τ + γr_τ) - (γo + γr))
    elseif Output == :SpecificH        #kJ/kg
        return R*T*τ*(γo_τ + γr_τ)
    elseif Output == :SpecificCP       #kJ/kgK
        return -R*τ^2*(γo_ττ + γr_ττ)
    elseif Output == :SpecificCV       #kJ/kgK
        return R*(-τ^2*(γo_ττ+γr_ττ) - (1+π*γr_π-τ*π*γr_πτ)^2/(1-π^2*γr_ππ))
    elseif Output == :SpeedOfSound     #m/s
        return sqrt((1000*R*T)*(1+2*π*γr_π+π^2*γr_π^2)/((1-π^2*γr_ππ)+(1+π*γr_π-τ*π*γr_πτ)^2/τ^2/(γo_ττ + γr_ττ)))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end

"""
    Region2meta

    Returns all the property values in region 2 in the metastable region,
    from the saturated-vapour line to the 5% equlibrium moisture line, a.k.a the
    practical Wilson line:
        %equilibrium moisture = [h - h_liq(P)] / [h_vap(P) - h_liq(P)]
    for pressures from P3 up to 10MPa.
    Pressures in MPa and temperature in [K]
    The property to be returned is specified in Output.
        :SpecificG        kJ/kg
        :SpecificF        kJ/kg
        :SpecificV        m3/kg
        :SpecificU        kJ/kg
        :SpecificS        kJ/kgK
        :SpecificH        kJ/kg
        :SpecificCP       kJ/kgK
        :SpecificCV       kJ/kgK
        :SpeedOfSound     m/s
"""
function Region2meta(Output::Symbol, P, T)
    Pstar = 1.0    #MPa
    Tstar = 540.0  #K

    no = [-0.969_372_683_930_49E1, #Updated from Region_2
           0.100_872_759_700_06E2, #Updated from Region_2
          -0.560_879_112_830_20E-2,
           0.714_527_380_814_55E-1,
          -0.407_104_982_239_28,
           0.142_408_191_714_44E1,
          -0.438_395_113_194_50E1,
          -0.284_086_324_607_72,
           0.212_684_637_533_07E-1]

    Jo = [ 0,
           1,
          -5,
          -4,
          -3,
          -2,
          -1,
           2,
           3]

    nr = [-0.733_622_601_865_06E-2,
          -0.882_238_319_431_46E-1,
          -0.723_345_552_132_45E-1,
          -0.408_131_785_344_55E-2,
           0.200_978_033_802_07E-2,
          -0.530_459_218_986_42E-1,
          -0.761_904_090_869_70E-2,
          -0.634_980_376_573_13E-2,
          -0.860_430_930_285_88E-1,
           0.753_215_815_227_70E-2,
          -0.792_383_754_461_39E-2,
          -0.228_881_607_784_47E-3,
          -0.264_565_014_828_10E-2]

    Ir = [1,
          1,
          1,
          1,
          2,
          2,
          2,
          3,
          3,
          4,
          4,
          5,
          5]

    Jr = [0,
          2,
          5,
          11,
          1,
          7,
          16,
          4,
          16,
          7,
          10,
          9,
          10]

    π = P / Pstar
    τ = Tstar / T

    γo    = log(π) + sum([no[i]*τ^Jo[i] for i=1:9])
    γo_π  = 1/π
    γo_ππ = -1/(π^2)
    γo_τ  = sum([no[i]*Jo[i]*τ^(Jo[i]-1) for i=1:9])
    γo_ττ = sum([no[i]*Jo[i]*(Jo[i]-1)*τ^(Jo[i]-2) for i=1:9])
    γo_πτ = 0

    γr    = sum([nr[i]*π^Ir[i]*(τ-0.5)^Jr[i] for i = 1:13])
    γr_π  = sum([nr[i]*Ir[i]*π^(Ir[i]-1)*(τ-0.5)^Jr[i] for i=1:13])
    γr_ππ = sum([nr[i]*Ir[i]*(Ir[i]-1)*π^(Ir[i]-2)*(τ-0.5)^Jr[i] for i=1:13])
    γr_τ  = sum([nr[i]*π^Ir[i]*Jr[i]*(τ-0.5)^(Jr[i]-1) for i=1:13])
    γr_ττ = sum([nr[i]*π^Ir[i]*Jr[i]*(Jr[i]-1)*(τ-0.5)^(Jr[i]-2) for i=1:13])
    γr_πτ = sum([nr[i]*Ir[i]*π^(Ir[i]-1)*Jr[i]*(τ-0.5)^(Jr[i]-1) for i = 1:13])

    if Output == :SpecificG            #kJ/kg
        return R*T*(γo + γr)
    elseif Output == :SpecificF        #kJ/kg
        return R*T*(γo + γr - π*(γo_π + γr_π)/1000)
    elseif Output == :SpecificV        #m3/kg
        return R*T*π*(γo_π + γr_π)/P/1000
    elseif Output == :SpecificU        #kJ/kg
        return R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))
    elseif Output == :SpecificS        #kJ/kgK
        return R*(τ*(γo_τ + γr_τ) - (γo + γr))
    elseif Output == :SpecificH        #kJ/kg
        return R*T*τ*(γo_τ + γr_τ)
    elseif Output == :SpecificCP       #kJ/kgK
        return -R*τ^2*(γo_ττ + γr_ττ)
    elseif Output == :SpecificCV       #kJ/kgK
        return R*(-τ^2*(γo_ττ+γr_ττ) - (1+π*γr_π-τ*π*γr_πτ)^2/(1-π^2*γr_ππ))
    elseif Output == :SpeedOfSound     #m/s
        return √((1000*R*T)*(1+2*π*γr_π+π^2*γr_π^2)/((1-π^2*γr_ππ)+(1+π*γr_π-τ*π*γr_πτ)^2/τ^2/(γo_ττ + γr_ττ)))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end

"""
    Region2a_TPh

    Returns T [K] from P[MPa] and h[kJ/kg] in Region 2a.
    Pressures in MPa and temperature in [K]
"""
function Region2a_TPh(P, h)
    hstar = 2000.0
    Pstar = 1.0

    n = [0.108_989_523_182_88E4,
         0.849_516_544_955_35E3,
        -0.107_817_480_918_26E3,
         0.331_536_548_012_63E2,
        -0.742_320_167_902_48E1,
         0.117_650_487_243_56E2,
         0.184_457_493_557_90E1,
        -0.417_927_005_496_24E1,
         0.624_781_969_358_12E1,
        -0.173_445_631_081_14E2,
        -0.200_581_768_620_96E3,
         0.271_960_654_737_96E3,
        -0.455_113_182_858_18E3,
         0.309_196_886_047_55E4,
         0.252_266_403_578_72E6,
        -0.617_074_228_683_39E-2,
        -0.310_780_466_295_83,
         0.116_708_730_771_07E2,
         0.128_127_984_040_46E9,
        -0.985_549_096_232_76E9,
         0.282_245_469_730_02E10,
        -0.359_489_714_107_03E10,
         0.172_273_499_131_97E10,
        -0.135_513_342_407_75E5,
         0.128_487_346_646_50E8,
         0.138_657_242_832_26E1,
         0.235_988_325_565_14E6,
        -0.131_052_365_450_54E8,
         0.739_998_354_747_66E4,
        -0.551_966_970_300_60E6,
         0.371_540_859_962_33E7,
         0.191_277_292_396_60E5,
        -0.415_351_648_356_34E6,
        -0.624_598_551_925_07E2]

    I = [0,
         0,
         0,
         0,
         0,
         0,
         1,
         1,
         1,
         1,
         1,
         1,
         1,
         1,
         1,
         2,
         2,
         2,
         2,
         2,
         2,
         2,
         2,
         3,
         3,
         4,
         4,
         4,
         5,
         5,
         5,
         6,
         6,
         7]

    J = [0,
         1,
         2,
         3,
         7,
         20,
         0,
         1,
         2,
         3,
         7,
         9,
         11,
         18,
         44,
         0,
         2,
         7,
         36,
         38,
         40,
         42,
         44,
         24,
         44,
         12,
         32,
         44,
         32,
         36,
         42,
         34,
         44,
         28]

    π = P / Pstar
    ηη = h / hstar - 2.1

    return sum([n[i]*π^I[i]*ηη^J[i] for i=1:34])
end

"""
    Region2b_TPh

    Returns T [K] from P[MPa] and h[kJ/kg] in Region 2b.
    Pressures in MPa and temperature in [K]
"""
function Region2b_TPh(P, h)
    hstar = 2000.0
    Pstar = 1.0

    n = [0.148_950_410_795_16E4,
         0.743_077_983_140_34E3,
        -0.977_083_187_978_37E2,
         0.247_424_647_056_74E1,
        -0.632_813_200_160_26,
         0.113_859_521_296_58E1,
        -0.478_118_636_486_25,
         0.852_081_234_315_44E-2,
         0.937_471_473_779_32,
         0.335_931_186_049_16E1,
         0.338_093_556_014_54E1,
         0.168_445_396_719_04,
         0.738_757_452_366_95,
        -0.471_287_374_361_86,
         0.150_202_731_397_07,
        -0.217_641_142_197_50E-2,
        -0.218_107_553_247_61E-1,
        -0.108_297_844_036_77,
        -0.463_333_246_358_12E-1,
         0.712_803_519_595_51E-4,
         0.110_328_317_899_99E-3,
         0.189_552_483_879_02E-3,
         0.308_915_411_605_37E-2,
         0.135_555_045_549_49E-2,
         0.286_402_374_774_56E-6,
        -0.107_798_573_575_12E-4,
        -0.764_627_124_548_14E-4,
         0.140_523_928_183_16E-4,
        -0.310_838_143_314_34E-4,
        -0.103_027_382_121_03E-5,
         0.282_172_816_350_40E-6,
         0.127_049_022_719_45E-5,
         0.738_033_534_682_92E-7,
        -0.110_301_392_389_09E-7,
        -0.814_563_652_078_33E-13,
        -0.251_805_456_829_62E-10,
        -0.175_652_339_694_07E-17,
         0.869_341_563_441_63E-14]

    I = [0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         1,
         1,
         1,
         1,
         1,
         1,
         1,
         1,
         2,
         2,
         2,
         2,
         3,
         3,
         3,
         3,
         4,
         4,
         4,
         4,
         4,
         4,
         5,
         5,
         5,
         6,
         7,
         7,
         9,
         9]

    J = [0,
         1,
         2,
         12,
         18,
         24,
         28,
         40,
         0,
         2,
         6,
         12,
         18,
         24,
         28,
         40,
         2,
         8,
         18,
         40,
         1,
         2,
         12,
         24,
         2,
         12,
         18,
         24,
         28,
         40,
         18,
         24,
         40,
         28,
         2,
         28,
         1,
         40]

    ππ = P / Pstar - 2
    ηη = h / hstar - 2.6

    return sum([n[i]*ππ^I[i]*ηη^J[i] for i=1:38])
end

"""
    Region2c_TPh

    Returns T [K] from P[MPa] and h[kJ/kg] in Region 2b.
    Pressures in MPa and temperature in [K]
"""
function Region2c_TPh(P, h)
    hstar = 2000.0
    Pstar = 1.0

    n = [-0.323_683_985_552_42E13,
          0.732_633_509_021_81E13,
          0.358_250_899_454_47E12,
         -0.583_401_318_515_90E12,
         -0.107_830_682_174_70E11,
          0.208_255_445_631_71E11,
          0.610_747_835_645_16E6,
          0.859_777_225_355_80E6,
         -0.257_457_236_041_70E5,
          0.310_810_884_227_14E5,
          0.120_823_158_659_36E4,
          0.482_197_551_092_55E3,
          0.379_660_012_724_86E1,
         -0.108_429_848_800_77E2,
         -0.453_641_726_766_60E-1,
          0.145_591_156_586_98E-12,
          0.112_615_974_072_30E-11,
         -0.178_049_822_406_86E-10,
          0.123_245_796_908_32E-6,
         -0.116_069_211_309_84E-5,
          0.278_463_670_885_54E-4,
         -0.592_700_384_741_76E-3,
          0.129_185_829_918_78E-2]

    I = [-7,
         -7,
         -6,
         -6,
         -5,
         -5,
         -2,
         -2,
         -1,
         -1,
          0,
          0,
          1,
          1,
          2,
          6,
          6,
          6,
          6,
          6,
          6,
          6,
          6]

    J = [0,
         4,
         0,
         2,
         0,
         2,
         0,
         1,
         0,
         2,
         0,
         1,
         4,
         8,
         4,
         0,
         1,
         4,
         10,
         12,
         16,
         20,
         22]

    ππ = P / Pstar + 25
    ηη = h / hstar - 1.8

    return sum([n[i]*ππ^I[i]*ηη^J[i] for i=1:23])
end

"""
    Region2a_TPs

    Returns T [K] from P[MPa] and s[kJ/kgK] in Region 2a.
    Pressures in MPa and temperature in [K]
"""
function Region2a_TPs(P, s)
    sstar = 2.0
    Pstar = 1.0

    n = [-0.392_359_838_619_84E6,
          0.515_265_738_272_70E6,
          0.404_824_431_610_48E5,
         -0.321_937_909_239_02E3,
          0.969_614_242_186_94E2,
         -0.228_678_463_717_73E2,
         -0.449_429_141_243_57E6,
         -0.501_183_360_201_66E4,
          0.356_844_635_600_15,
          0.442_353_358_481_90E5,
         -0.136_733_888_117_08E5,
          0.421_632_602_078_64E6,
          0.225_169_258_374_75E5,
          0.474_421_448_656_46E3,
         -0.149_311_307_976_47E3,
         -0.197_811_263_204_52E6,
         -0.235_543_994_707_60E5,
         -0.190_706_163_020_76E5,
          0.553_756_698_831_64E5,
          0.382_936_914_373_63E4,
         -0.603_918_605_805_67E3,
          0.193_631_026_203_31E4,
          0.426_606_436_986_10E4,
         -0.597_806_388_727_18E4,
         -0.704_014_639_268_62E3,
          0.338_367_841_075_53E3,
          0.208_627_866_351_87E2,
          0.338_341_726_561_96E-1,
         -0.431_244_284_148_93E-4,
          0.166_537_913_564_12E3,
         -0.139_862_920_558_98E3,
         -0.788_495_479_998_72,
          0.721_324_117_538_72E-1,
         -0.597_548_393_982_83E-2,
         -0.121_413_589_539_04E-4,
          0.232_270_967_338_71E-6,
         -0.105_384_635_661_94E2,
          0.207_189_254_965_02E1,
         -0.721_931_552_604_27E-1,
          0.207_498_870_811_20E-6,
         -0.183_406_579_113_79E-1,
          0.290_362_723_486_96E-6,
          0.210_375_278_936_19,
          0.256_812_397_299_99E-3,
         -0.127_990_029_337_81E-1,
         -0.821_981_026_520_18E-5]

    I = [-1.5,
         -1.5,
         -1.5,
         -1.5,
         -1.5,
         -1.5,
         -1.25,
         -1.25,
         -1.25,
         -1.0,
         -1.0,
         -1.0,
         -1.0,
         -1.0,
         -1.0,
         -0.75,
         -0.75,
         -0.5,
         -0.5,
         -0.5,
         -0.5,
         -0.25,
         -0.25,
         -0.25,
         -0.25,
         0.25,
         0.25,
         0.25,
         0.25,
         0.5,
         0.5,
         0.5,
         0.5,
         0.5,
         0.5,
         0.5,
         0.75,
         0.75,
         0.75,
         0.75,
         1.0,
         1.0,
         1.25,
         1.25,
         1.5,
         1.5]

    J = [-24,
         -23,
         -19,
         -13,
         -11,
         -10,
         -19,
         -15,
         -6,
         -26,
         -21,
         -17,
         -16,
         -9,
         -8,
         -15,
         -14,
         -26,
         -13,
         -9,
         -7,
         -27,
         -25,
         -11,
         -6,
         1,
         4,
         8,
         11,
         0,
         1,
         5,
         6,
         10,
         14,
         16,
         0,
         4,
         9,
         17,
         7,
         18,
         3,
         15,
         5,
         18]

    π = P / Pstar
    σσ = s / sstar - 2.0

    return sum([n[i]*π^I[i]*σσ^J[i] for i=1:46])
end

"""
    Region2b_TPs

    Returns T [K] from P[MPa] and s[kJ/kgK] in Region 2a.
    Pressures in MPa and temperature in [K]
"""
function Region2b_TPs(P, s)
    sstar = 0.7853
    Pstar = 1.0

    n = [0.316_876_650_834_97E6,
         0.208_641_758_818_58E2,
        -0.398_593_998_035_99E6,
        -0.218_160_585_188_77E2,
         0.223_697_851_942_42E6,
        -0.278_417_034_458_17E4,
         0.992_074_360_714_80E1,
        -0.751_975_122_991_57E5,
         0.297_086_059_511_58E4,
        -0.344_068_785_485_26E1,
         0.388_155_642_491_15,
         0.175_112_950_857_50E5,
        -0.142_371_128_544_49E4,
         0.109_438_033_641_67E1,
         0.899_716_193_084_95,
        -0.337_597_400_989_58E4,
         0.471_628_858_183_55E3,
        -0.191_882_419_936_79E1,
         0.410_785_804_921_96,
        -0.334_653_781_720_97,
         0.138_700_347_775_05E4,
        -0.406_633_261_958_38E3,
         0.417_273_471_596_10E2,
         0.219_325_494_345_32E1,
        -0.103_200_500_090_77E1,
         0.358_829_435_167_03,
         0.525_114_537_260_66E-2,
         0.128_389_164_507_05E2,
        -0.286_424_372_193_81E1,
         0.569_126_836_648_55,
        -0.999_629_545_849_31E-1,
        -0.326_320_377_784_59E-2,
         0.233_209_225_767_23E-3,
        -0.153_348_098_574_50,
         0.290_722_882_399_02E-1,
         0.375_347_027_411_67E-3,
         0.172_966_917_024_11E-2,
        -0.385_560_508_445_04E-3,
        -0.350_177_122_926_08E-4,
        -0.145_663_936_314_92E-4,
         0.564_208_572_672_69E-5,
         0.412_861_500_746_05E-7,
        -0.206_846_711_188_24E-7,
         0.164_093_936_747_25E-8]

    I = [-6,
         -6,
         -5,
         -5,
         -4,
         -4,
         -4,
         -3,
         -3,
         -3,
         -3,
         -2,
         -2,
         -2,
         -2,
         -1,
         -1,
         -1,
         -1,
         -1,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          1,
          1,
          1,
          1,
          1,
          1,
          2,
          2,
          2,
          3,
          3,
          3,
          4,
          4,
          5,
          5,
          5]

    J = [0,
         11,
         0,
         11,
         0,
         1,
         11,
         0,
         1,
         11,
         12,
         0,
         1,
         6,
         10,
         0,
         1,
         5,
         8,
         9,
         0,
         1,
         2,
         4,
         5,
         6,
         9,
         0,
         1,
         2,
         3,
         7,
         8,
         0,
         1,
         5,
         0,
         1,
         3,
         0,
         1,
         0,
         1,
         2]

    π = P / Pstar
    σσ = 10.0 - s / sstar

    return sum([n[i]*π^I[i]*σσ^J[i] for i=1:44])
end

"""
    Region2c_TPs

    Returns T [K] from P[MPa] and s[kJ/kgK] in Region 2a.
    Pressures in MPa and temperature in [K]
"""
function Region2c_TPs(P, s)
    sstar = 2.9251
    Pstar = 1.0

    n = [0.909_685_010_053_65E3,
         0.240_456_670_884_20E4,
        -0.591_623_263_871_30E3,
         0.541_454_041_280_74E3,
        -0.270_983_084_111_92E3,
         0.979_765_250_979_26E3,
        -0.469_667_729_594_35E3,
         0.143_992_746_047_23E2,
        -0.191_042_042_304_29E2,
         0.532_991_671_119_71E1,
        -0.212_529_753_759_34E2,
        -0.311_473_344_137_60,
         0.603_348_408_946_23,
        -0.427_648_397_025_09E-1,
         0.581_855_972_552_59E-2,
        -0.145_970_082_847_53E-1,
         0.566_311_756_310_27E-2,
        -0.761_558_645_845_77E-4,
         0.224_403_429_193_32E-3,
        -0.125_610_950_134_13E-4,
         0.633_231_326_609_34E-6,
        -0.205_419_896_753_75E-5,
         0.364_053_703_900_82E-7,
        -0.297_598_977_892_15E-8,
         0.101_366_185_297_63E-7,
         0.599_257_196_923_51E-11,
        -0.206_778_701_051_64E-10,
        -0.208_742_781_818_86E-10,
         0.101_621_668_250_89E-9,
        -0.164_298_282_813_47E-9]

    I = [-2,
         -2,
         -1,
          0,
          0,
          0,
          0,
          1,
          1,
          1,
          1,
          2,
          2,
          2,
          3,
          3,
          3,
          4,
          4,
          4,
          5,
          5,
          5,
          6,
          6,
          7,
          7,
          7,
          7,
          7]

    J = [0,
         1,
         0,
         0,
         1,
         2,
         3,
         0,
         1,
         3,
         4,
         0,
         1,
         2,
         0,
         1,
         5,
         0,
         1,
         4,
         0,
         1,
         2,
         0,
         1,
         0,
         1,
         3,
         4,
         5]

    π = P / Pstar
    σσ = 2.0 - s / sstar

    return sum([n[i]*π^I[i]*σσ^J[i] for i=1:30])
end


"""
    Region3_ρ

    Returns all the property values in region 3.
    623.15K ≤ T ≤ T(P) from B23-model P(T) from B23-model ≤ P ≤ 100MPa
    Pressures in MPa and temperature in [K]
    The property to be returned is specified in Output.
        :SpecificG        kJ/kg
        :SpecificF        kJ/kg
        :Pressure         MPa
        :SpecificU        kJ/kg
        :SpecificS        kJ/kgK
        :SpecificH        kJ/kg
        :SpecificCP       kJ/kgK
        :SpecificCV       kJ/kgK
        :SpeedOfSound     m/s
"""
function Region3_ρ(Output::Symbol, ρ, T)
    ρstar = ρc      #kg/m3
    Tstar = Tc      #K

    n = [0.106_580_700_285_13E1,
        -0.157_328_452_902_39E2,
         0.209_443_969_743_07E2,
        -0.768_677_078_787_16E1,
         0.261_859_477_879_54E1,
        -0.280_807_811_486_20E1,
         0.120_533_696_965_17E1,
        -0.845_668_128_125_02E-2,
        -0.126_543_154_777_14E1,
        -0.115_244_078_066_81E1,
         0.885_210_439_843_18,
        -0.642_077_651_816_07,
         0.384_934_601_866_71,
        -0.852_147_088_242_06,
         0.489_722_815_418_77E1,
        -0.305_026_172_569_65E1,
         0.394_205_368_791_54E-1,
         0.125_584_084_243_08,
        -0.279_993_296_987_10,
         0.138_997_995_694_60E1,
        -0.201_899_150_235_70E1,
        -0.821_476_371_739_63E-2,
        -0.475_960_357_349_23,
         0.439_840_744_735_00E-1,
        -0.444_764_354_287_39,
         0.905_720_707_197_33,
         0.705_224_500_879_67,
         0.107_705_126_263_32,
        -0.329_136_232_589_54,
        -0.508_710_620_411_58,
        -0.221_754_008_730_96E-1,
         0.942_607_516_650_92E-1,
         0.164_362_784_479_61,
        -0.135_033_722_413_48E-1,
        -0.148_343_453_524_72E-1,
         0.579_229_536_280_84E-3,
         0.323_089_047_037_11E-2,
         0.809_648_029_962_15E-4,
        -0.165_576_797_950_37E-3,
        -0.449_238_990_618_15E-4]

    I = [0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         1,
         1,
         1,
         1,
         2,
         2,
         2,
         2,
         2,
         2,
         3,
         3,
         3,
         3,
         3,
         4,
         4,
         4,
         4,
         5,
         5,
         5,
         6,
         6,
         6,
         7,
         8,
         9,
         9,
         10,
         10,
         11]

    J = [0,
         0,
         1,
         2,
         7,
         10,
         12,
         23,
         2,
         6,
         15,
         17,
         0,
         2,
         6,
         7,
         22,
         26,
         0,
         2,
         4,
         16,
         26,
         0,
         2,
         4,
         26,
         1,
         3,
         26,
         0,
         2,
         26,
         2,
         26,
         2,
         26,
         0,
         1,
         26]

    δ = ρ / ρstar
    τ = Tstar / T

    ϕ    =  n[1]*log(δ) + sum([n[i]*(δ^I[i])*(τ^J[i]) for i=2:40])
    ϕ_δ  =  n[1]/δ + sum([n[i]*I[i]*(δ^(I[i]-1))*(τ^J[i]) for i=2:40])
    ϕ_δδ = -n[1]/(δ^2) + sum([n[i]*I[i]*(I[i]-1)*(δ^(I[i]-2))*(τ^J[i]) for i=2:40])
    ϕ_τ  =  sum([n[i]*(δ^I[i])*J[i]*(τ^(J[i]-1)) for i=2:40])
    ϕ_ττ =  sum([n[i]*(δ^I[i])*J[i]*(J[i]-1)*(τ^(J[i]-2)) for i=2:40])
    ϕ_δτ =  sum([n[i]*I[i]*(δ^(I[i]-1))*J[i]*(τ^(J[i]-1)) for i=2:40])

    if Output == :SpecificG            #kJ/kg
        return R*T*(δ*ϕ_δ + ϕ)
    elseif Output == :SpecificF        #kJ/kg
        return R*T*ϕ
    elseif Output == :Pressure         #MPa
        return ρ*R*T*δ*ϕ_δ/1000
    elseif Output == :SpecificU        #kJ/kg
        return R*T*τ*ϕ_τ
    elseif Output == :SpecificS        #kJ/kgK
        return R*(τ*ϕ_τ - ϕ)
    elseif Output == :SpecificH        #kJ/kg
        return R*T*(τ*ϕ_τ + δ*ϕ_δ)
    elseif Output == :SpecificCP       #kJ/kgK
        return R*(-τ^2*ϕ_ττ+((δ*ϕ_δ - δ*τ*ϕ_δτ)^2)/(2*δ*ϕ_δ + δ^2*ϕ_δδ))
    elseif Output == :SpecificCV       #kJ/kgK
        return R*(-τ^2*ϕ_ττ)
    elseif Output == :SpeedOfSound     #m/s
        return sqrt(1000*R*T*(2*δ*ϕ_δ + δ^2*ϕ_δδ - ((δ*ϕ_δ - δ*τ*ϕ_δτ)^2)/(τ^2*ϕ_ττ)))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end


"""
    Initialise from ideal gas law, then use root finder to calculate P by iterating on Region3ρ.
    Pass through properties from Region3ρ
"""
function Region3(Output::Symbol, P, T)
    # Start with a higher density to fix the problems in convergence.
    # TODO Replace this with a density from Peng-Robinson or similar EoS
    # ρ0 = 1000*P/(R*T) #Starting value from ideal gas
    ρ0 = 500.0

    f(ρ) = Region3_ρ(:Pressure, ρ, T) - P
    ρ = Roots.find_zero(f, ρ0)

    if Output == :SpecificV
        return 1.0 / ρ
    else
        return Region3_ρ(Output, ρ, T)
    end
end


"""
    Region4

    Returns the vapour-liquid phase boundary (saturation line).
    Valid from triple point to critical point
    273.15K ≤ T ≤ 647.096K
    InputType is either :T or :P to indicate that InputValue is temperature [K]
    or pressure [MPa]. The complimentary value is returned.
"""
function Region4(InputType::Symbol, InputValue)
    n = [0.116_705_214_527_67E4,
        -0.724_213_167_032_06E6,
        -0.170_738_469_400_92E2,
         0.120_208_247_024_70E5,
        -0.323_255_503_223_33E7,
         0.149_151_086_135_30E2,
        -0.482_326_573_615_91E4,
         0.405_113_405_420_57E6,
        -0.238_555_575_678_49,
         0.650_175_348_447_98E3]

    if InputType == :T
        Θ = InputValue + n[9]/(InputValue - n[10])
        A =      Θ^2 + n[1]*Θ + n[2]
        B = n[3]*Θ^2 + n[4]*Θ + n[5]
        C = n[6]*Θ^2 + n[7]*Θ + n[8]

        P = (2C / (-B + √(B^2 - 4*A*C)))^4
        return P
    elseif InputType == :P
        β = InputValue^0.25
        E =      β^2 + n[3]*β + n[6]
        F = n[1]*β^2 + n[4]*β + n[7]
        G = n[2]*β^2 + n[5]*β + n[8]
        D = 2G / (-F - √(F^2 - 4*E*G))

        T = (n[10]+D-√((n[10]+D)^2 - 4(n[9]+n[10]*D)))/2.0
        return T
    else
        throw(DomainError(InputType, "Unknown input. Expecting :T or :P"))
    end
end


"""
    Region5

    Returns all the property values in region 5.
    1073.15K ≤ T ≤ 2273.15K  0 ≤ P ≤ 50MPa
    Pressures in MPa and temperature in [K]

    The property to be returned is specified in Output.
        :SpecificG        kJ/kg
        :SpecificF        kJ/kg
        :SpecificV        m3/kg
        :SpecificU        kJ/kg
        :SpecificS        kJ/kgK
        :SpecificH        kJ/kg
        :SpecificCP       kJ/kgK
        :SpecificCV       kJ/kgK
        :SpeedOfSound     m/s
"""
function Region5(Output::Symbol, P, T)
    Pstar = 1.0     #MPa
    Tstar = 1000.0  #K

    no = [-0.131_799_836_742_01E2,
           0.685_408_416_344_34E1,
          -0.248_051_489_334_66E-1,
           0.369_015_349_803_33,
          -0.311_613_182_139_25E1,
          -0.329_616_265_389_17]

    Jo = [0,
          1,
         -3,
         -2,
         -1,
          2]

    nr = [0.157_364_048_552_59E-2,
          0.901_537_616_739_44E-3,
         -0.502_700_776_776_48E-2,
          0.224_400_374_094_85E-5,
         -0.411_632_754_534_71E-5,
          0.379_194_548_229_55E-7]

    Ir = [1,
          1,
          1,
          2,
          2,
          3]

    Jr = [1,
          2,
          3,
          3,
          9,
          7]

    π = P / Pstar
    τ = Tstar / T

    γo    =  log(π) + sum([no[i]*(τ^Jo[i]) for i=1:6])
    γo_π  =  1/π
    γo_ππ = -1/(π^2)
    γo_τ  =          sum([no[i]*Jo[i]*(τ^(Jo[i]-1)) for i=1:6])
    γo_ττ =          sum([no[i]*Jo[i]*(Jo[i]-1)*(τ^(Jo[i]-2)) for i=1:6])
    γo_πτ = 0

    γr    = sum([nr[i]*(π^Ir[i])*(τ^Jr[i]) for i=1:6])
    γr_π  = sum([nr[i]*Ir[i]*(π^(Ir[i]-1))*(τ^Jr[i]) for i=1:6])
    γr_ππ = sum([nr[i]*Ir[i]*(Ir[i]-1)*(π^(Ir[i]-2))*(τ^Jr[i]) for i=1:6])
    γr_τ  = sum([nr[i]*(π^Ir[i])*Jr[i]*(τ^(Jr[i]-1)) for i=1:6])
    γr_ττ = sum([nr[i]*(π^Ir[i])*Jr[i]*(Jr[i]-1)*(τ^(Jr[i]-2)) for i=1:6])
    γr_πτ = sum([nr[i]*Ir[i]*(π^(Ir[i]-1))*Jr[i]*(τ^(Jr[i]-1)) for i=1:6])

    if Output == :SpecificG            #kJ/kg
        return R*T*(γo + γr)
    elseif Output == :SpecificF            #kJ/kg
        return R*T*((γo + γr) - π*(γo_π + γr_π)/1000)
    elseif Output == :SpecificV        #m3/kg
        return R*T*π*(γo_π + γr_π)/P/1000
    elseif Output == :SpecificU        #kJ/kg
        return R*T*(τ*(γo_τ + γr_τ) - π*(γo_π + γr_π))
    elseif Output == :SpecificS        #kJ/kgK
        return R*(τ*(γo_τ + γr_τ) - (γo + γr))
    elseif Output == :SpecificH        #kJ/kg
        return R*T*τ*(γo_τ + γr_τ)
    elseif Output == :SpecificCP       #kJ/kgK
        return -R*τ^2*(γo_ττ + γr_ττ)
    elseif Output == :SpecificCV       #kJ/kgK
        return R*((-τ^2)*(γo_ττ + γr_ττ) - (1 + π*γr_π - τ*π*γr_πτ)^2/(1-(π^2)*γr_ππ))
    elseif Output == :SpeedOfSound     #m/s
        return √(1000*R*T*(1 + 2*π*γr_π + (π^2)*(γr_π^2))/((1 - (π^2)*(γr_ππ))+((1+π*γr_π-τ*π*γr_πτ)^2)/(τ^2*(γo_ττ + γr_ττ))))
    else
        throw(DomainError(Output, "Unknown value requested."))
    end
end


"""
    RegionID

    Identifies the applicable region, based on specified T and P.
    This allows the correct function to be called to retrieve properties.
"""

function RegionID(P, T)::Symbol
    #=
    Region 1:
        273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa

    Region 2:
        273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
        623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
        863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa

    Region 2meta:
    From the saturated-vapour line to the 5% equlibrium moisture line, a.k.a the
    practical Wilson line:
        %equilibrium moisture = [h - h_liq(P)] / [h_vap(P) - h_liq(P)]
    for pressures from the triple point (273.16K, 611.657MPa) up to 10MPa.
    These equations only used on demand.

    Region 3:
        623.15K ≤ T ≤ T(P) from B23-model P(T) from B23-model ≤ P ≤ 100MPa

    Region 4:
    Returns the vapour-liquid phase boundary (saturation line).
    Valid from triple point to critical point
    273.15K ≤ T ≤ 647.096K
    Since this is a single floating point number, Region 4 is never returned and only used on demand.

    Region 5:emperature in [K]
    =#

    #Check region 5 first:
    if 1073.15 ≤ T ≤ 2273.15 && 0 ≤ P ≤ 50
        return :Region5
    elseif T < 273.15 || T > 1073.15 || P < 0 || P > 100
        throw(DomainError((P, T), "Pressure/Temperature outside valid ranges.")) #Outside the regions where equations are available.
    elseif 273.15 ≤ T ≤ 623.15
        if Psat(T) ≤ P ≤ 100
            return :Region1
        else
            return :Region2
        end
    elseif 623.15 <  T ≤ 863.15
        if 0 <  P ≤ B23(:T, T)
            return :Region2
        else
            return :Region3
        end
    else
        return :Region2
    end
end

"""
    RegionID_Ph

    Identifies the applicable region, based on specified P [MPa] and h [kJ/kg].
    This allows the correct function to be called to retrieve properties.
    Backwards equations are only available for Regions 1 and 2 and the vapour-liquid
    equilibrium line, which is only used when explicitly called.

    As there are no explicit domains defined for (T, h), the method is iterative.
    The reverse equation for T = f(P,h) is called for each region and the result
    is checked against the bounds for the region.
"""
function RegionID_Ph(P, h)::Symbol

    #=
        Region 1:
        273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa

        Region 2:
        273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
        623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
        863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa
        For reverse:
            2a: P ≤ 4Mpa
            2b: s ≥ 5.85 kJ/kgK (or use function B2bc if P,h specificed)
            2c: s ≥ 5.85 kJ/kgK (or use function B2bc if P,h specificed)
    =#

    # Check overall region first
    if P > 100
        throw(DomainError(P, "Pressure not in valid ranges."))
    end

    T = Region1_TPh(P, h)
    if 273.15 ≤ T ≤ 623.15
        # could be Region 1
        if P ≥ Psat(T)
            return :Region1 #, Region1(:SpecificH, P, T) # Return forward h for consistency check
        end
    end

    if P ≤ 4.0
        # could be Region 2a
        T = Region2a_TPh(P, h)
        if Tsat(P) ≤ T ≤ 1073.15
            return :Region2a #, Region2(:SpecificH, P, T) # Return forward h for consistency check
        else
            throw(DomainError((P, T), "Pressure/Temperature outside valid ranges.")) # Only other region with backwards mdoels for P ≤ 4.0 is Region 1, which is already eliminated
        end
    else
        htest = B2bc(:P, P)
        if h < htest
            # could be Region 2c
            T = Region2c_TPh(P, h)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :Region2c #, Region2(:SpecificH, P, T) # Return forward h for consistency check
                end
            elseif B23(:P, P) ≤ T ≤ 1073.15
                return :Region2c #, Region2(:SpecificH, P, T) # Return forward h for consistency check
            else
                throw(DomainError((P, T), "Pressure/Temperature outside valid ranges."))
            end
        else
            # could be Region 2b
            T = Region2b_TPh(P, h)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :Region2b #, Region2(:SpecificH, P, T) # Return forward h for consistency check
                end
            elseif B23(:P, P) ≤ T ≤ 1073.15
                return :Region2b #, Region2(:SpecificH, P, T) # Return forward h for consistency check
            else
                throw(DomainError((P, T), "Pressure/Temperature outside valid ranges."))
            end
        end
    end
    throw(DomainError((P, T), "Pressure/Temperature outside valid ranges."))
end


"""
    RegionID_Ps

    Identifies the applicable region, based on specified P [MPa] and h [kJ/kg].
    This allows the correct function to be called to retrieve properties.
    Backwards equations are only available for Regions 1 and 2 and the vapour-liquid
    equilibrium line, which is only used when explicitly called.

    As there are no explicit domains defined for (T, s), the method is iterative.
    The reverse equation for T = f(P,s) is called for each region and the result
    is checked against the bounds for the region.
"""
function RegionID_Ps(P, s)::Symbol

    #=
        Region 1:
        273.15K ≤ T ≤ 623.15K  Psat(T) ≤ P ≤ 100MPa

        Region 2:
        273.15K ≤ T ≤ 623.15K  0 ≤ P ≤ Psat(T)
        623.15K <  T ≤ 863.15K  0 <  P ≤ P(T) from B23-model
        863.15K <  T ≤ 1073.15K 0 <  P ≤ 100MPa
        For reverse:
            2a: P ≤ 4Mpa
            2c: s ≥ 5.85 kJ/kgK
            2b: s ≥ 5.85 kJ/kgK
    =#

    # Check overall region first
    if P > 100
        throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))
    end

    T = Region1_TPs(P, s)
    if 273.15 ≤ T ≤ 623.15
        # could be Region 1
        if P ≥ Psat(T)
            return :Region1 #, Region1(:SpecificS, P, T) # Return forward s for consistency check
        end
    end

    if P ≤ 4.0
        # could be Region 2a
        T = Region2a_TPs(P, s)
        if Tsat(P) ≤ T ≤ 1073.15
            return :Region2a #, Region2(:SpecificS, P, T) # Return forward s for consistency check
        else
            throw(DomainError((P, s), "Pressure/entropy not in valid ranges.")) # Only other region with backwards mdoels for P ≤ 4.0 is Region 1, which is already eliminated
        end
    else
        if s < 5.85
            # could be Region 2c
            T = Region2c_TPs(P, s)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :Region2 #c, Region2(:SpecificS, P, T) # Return forward s for consistency check
                end
            elseif B23(:P, P) ≤ T ≤ 1073.15
                return :Region2c #, Region2(:SpecificS, P, T) # Return forward s for consistency check
            else
                throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))
            end
        else
            # could be Region 2b
            T = Region2b_TPs(P, s)
            if P ≤ Psat(623.15)
                if Tsat(P) ≤ T ≤ 1073.15
                    return :Region2b #, Region2(:SpecificS, P, T) # Return forward s for consistency check
                end
            elseif B23(:P, P) ≤ T ≤ 1073.15
                return :Region2b #, Region2(:SpecificS, P, T) # Return forward s for consistency check
            else
                throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))

            end
        end
    end
    throw(DomainError((P, s), "Pressure/entropy not in valid ranges."))
end


"""
    Psat2

    Do not call this function, but rather use Psat().
    Returns the saturated vapour pressure according to IAPWS SR1-86(1992)
    with T in K and P in MPa.
    This replaces the Region 4 equations for Psat and is used to calcualate
    the saturated enthalpy and entropy to match the values in IAPWS SR1-86(1992)
    It is not used in Tsat(), as the difference in vapour pressure is in the
    8th decimal at the boiling point.
"""
function Psat2(T)

    a = [ -7.85951783000,
           1.84408259000,
          -11.78664970000,
          22.68074110000,
         -15.96187190000,
           1.80122502000
        ]

    θ = T/Tc
    τ = 1 - θ

    return Pc*exp(Tc/T*(a[1]*τ + a[2]*τ^1.5 + a[3]*τ^3 + a[4]*τ^3.5
                      + a[5]*τ^4 + a[6]*τ^7.5))
end


#============================================================================
=============================================================================
                            Exported functions:
                            ~~~~~~~~~~~~~~~~~~
============================================================================
============================================================================#


"""
    Psat

    Utility function that returns the vapour pressure at temperature T [K].
    If inputs have associated units, the value is returned with associated units
    of MPa via Uniful.jl.
"""
function Psat(T)
    if T3 ≤ T ≤ Tc
        return Region4(:T, T)
    else
        throw(DomainError(T, "Temperature not between tripple and critical points."))
    end
end


function Psat(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        return Region4(:T, T.val)*u"MPa"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points."))
    end
end


"""
    Tsat

    Utility function that returns the vapour pressure at pressure P [MPa].
    If inputs have associated units, the value is returned with associated
    units of K via Uniful.jl.
"""
function Tsat(P)
    if 0.611_213E-3 ≤ P ≤ Pc # Note that botton limit not exactly P3
        return Region4(:P, P)
    else
        throw(DomainError(P, "Pressure not between tripple and critical points"))
    end
end


function Tsat(P::Q) where Q <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end

    if 0.611_213E-3 ≤ P.val ≤ Pc
        return Region4(:P, P.val)*u"K"
    else
        throw(DomainError(P, "Pressure not between tripple and critical points"))
    end
end


"""
    SpecificG

    Utility function that returns the Gibbs free energy [kJ/kgK] from P [MPa]
    and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificG(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificG, P, T)
    elseif Region == :Region2
        return Region2(:SpecificG, P, T)
    elseif Region == :Region3
        return Region3(:SpecificG, P, T)
    elseif Region == :Region5
        return Region5(:SpecificG, P, T)
    end
end


function SpecificG(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
        T = 1.0*uconvert(u"K", T)
        Region = RegionID(P.val, T.val)

        if Region == :Region1
            return Region1(:SpecificG, P.val, T.val)*u"kJ/kg"
        elseif Region == :Region2
            return Region2(:SpecificG, P.val, T.val)*u"kJ/kg"
        elseif Region == :Region3
            return Region3(:SpecificG, P.val, T.val)*u"kJ/kg"
        elseif Region == :Region5
            return Region5(:SpecificG, P.val, T.val)*u"kJ/kg"
        end
    catch
        throw(UnitsError((P,T), "Invalid input units."))
    end

end


"""
    SpecificG_Ph

    Utility function that returns the Gibbs free energy [kJ/kgK] from P [MPa] and h [kJ/kg]
    The explicit backwards equations are only available in regions 1 and 2. Input outside these
    will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificG_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificG, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificG, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificG, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificG, P, T)
    end
end


function SpecificG_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificG, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificG, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificG, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificG, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificG_Ps

    Utility function that returns the Gibbs free energy [kJ/kgK] from P [MPa]
    and s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificG_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificG, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificG, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificG, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificG, P, T)
    end
end


function SpecificG_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificG, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificG, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificG, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificG, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificF

    Utility function that returns the Helmholtz free energy [kJ/kgK] from
    P [MPa] and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificF(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificF, P, T)
    elseif Region == :Region2
        return Region2(:SpecificF, P, T)
    elseif Region == :Region3
        return Region3(:SpecificF, P, T)
    elseif Region == :Region5
        return Region5(:SpecificF, P, T)
    end
end


function SpecificF(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificF, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region2
        return Region2(:SpecificF, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region3
        return Region3(:SpecificF, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region5
        return Region5(:SpecificF, P.val, T.val)*u"kJ/kg"
    end
end


"""
    SpecificF_Ph

    Utility function that returns the Helmholtz free energy [kJ/kgK] from
    P [MPa] and h [kJ/kg].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificF_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificF, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificF, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificF, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificF, P, T)
    end
end


function SpecificF_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificF, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificF, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificF, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificF, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificF_Ps

    Utility function that returns the Helmholtz free energy [kJ/kgK] from
    P [MPa] and s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificF_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificF, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificF, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificF, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificF, P, T)
    end
end


function SpecificF_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificF, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificF, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificF, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificF, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificV

    Utility function that returns the specific volume [m3/kg] from P [MPa] and T [K].
    If inputs have associated units, the value is returned with associated
    units of m3/kg via Uniful.jl.
"""
function SpecificV(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificV, P, T)
    elseif Region == :Region2
        return Region2(:SpecificV, P, T)
    elseif Region == :Region3
        return Region3(:SpecificV, P, T)
    elseif Region == :Region5
        return Region5(:SpecificV, P, T)
    end
end


function SpecificV(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificV, P.val, T.val)*u"m^3/kg"
    elseif Region == :Region2
        return Region2(:SpecificV, P.val, T.val)*u"m^3/kg"
    elseif Region == :Region3
        return Region3(:SpecificV, P.val, T.val)*u"m^3/kg"
    elseif Region == :Region5
        return Region5(:SpecificV, P.val, T.val)*u"m^3/kg"
    end
end


"""
    SpecificV_Ph

    Utility function that returns the specific volume [m3/kg] from P [MPa] and h [kJ/kg]
    The explicit backwards equations are only available in regions 1 and 2. Input outside these
    will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of m3/kg via Uniful.jl.
"""
function SpecificV_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificV, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificV, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificV, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificV, P, T)
    end
end


function SpecificV_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificV, P.val, T)*u"m^3/kg"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificV, P.val, T)*u"m^3/kg"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificV, P.val, T)*u"m^3/kg"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificV, P.val, T)*u"m^3/kg"
    end
end


"""
    SpecificV_Ps

    Utility function that returns the specific volume [m3/kg] from P [MPa] and
    s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of m3/kg via Uniful.jl.
"""
function SpecificV_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificV, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificV, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificV, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificV, P, T)
    end
end


function SpecificV_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificV, P.val, T)*u"m^3/kg"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificV, P.val, T)*u"m^3/kg"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificV, P.val, T)*u"m^3/kg"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificV, P.val, T)*u"m^3/kg"
    end
end


"""
    SpecificU

    Utility function that returns the specific internal energy [kJ/kg] from
    P [MPa] and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificU(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificU, P, T)
    elseif Region == :Region2
        return Region2(:SpecificU, P, T)
    elseif Region == :Region3
        return Region3(:SpecificU, P, T)
    elseif Region == :Region5
        return Region5(:SpecificU, P, T)
    end
end


function SpecificU(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificU, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region2
        return Region2(:SpecificU, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region3
        return Region3(:SpecificU, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region5
        return Region5(:SpecificU, P.val, T.val)*u"kJ/kg"
    end
end


"""
    SpecificU_Ph

    Utility function that returns the internal energy [kJ/kgK] from P [MPa] and
    h [kJ/kg].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificU_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificU, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificU, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificU, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificU, P, T)
    end
end


function SpecificU_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificU, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificU, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificU, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificU, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificU_Ps

    Utility function that returns the internal energy [kJ/kgK] from P [MPa] and
    s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.
"""
function SpecificU_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificU, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificU, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificU, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificU, P, T)
    end
end


function SpecificU_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificU, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificU, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificU, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificU, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificS

    Utility function that returns the specific entropy [kJ/kgK] from P [MPa]
    and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificS(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificS, P, T)
    elseif Region == :Region2
        return Region2(:SpecificS, P, T)
    elseif Region == :Region3
        return Region3(:SpecificS, P, T)
    elseif Region == :Region5
        return Region5(:SpecificS, P, T)
    end
end


function SpecificS(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificS, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region2
        return Region2(:SpecificS, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region3
        return Region3(:SpecificS, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region5
        return Region5(:SpecificS, P.val, T.val)*u"kJ/kg/K"
    end
end


"""
    SpecificS_Ph

    Utility function that returns the specific entropy [kJ/kgK] from P [MPa]
    and h [kJ/kg].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificS_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificS, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificS, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificS, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificS, P, T)
    end
end


function SpecificS_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificS, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificS, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificS, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificS, P.val, T)*u"kJ/kg/K"
    end
end


"""
    SpecificH

    Utility function that returns the specific enthalpy [kJ/kg] from P [MPa]
    and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.

    Note: Do not use this function to attempt to find values on the saturation
    line, as any rounding errors anywhere will result in a point on either side
    of the phase boundary with large differences in enthalpy. Use SatHL/SatHV
    for saturated phase enthalpies.
"""
function SpecificH(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificH, P, T)
    elseif Region == :Region2
        return Region2(:SpecificH, P, T)
    elseif Region == :Region3
        return Region3(:SpecificH, P, T)
    elseif Region == :Region5
        return Region5(:SpecificH, P, T)
    end
end


function SpecificH(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificH, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region2
        return Region2(:SpecificH, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region3
        return Region3(:SpecificH, P.val, T.val)*u"kJ/kg"
    elseif Region == :Region5
        return Region5(:SpecificH, P.val, T.val)*u"kJ/kg"
    end
end


"""
    SpecificH_Ps

    Utility function that returns the speciic enthalpy [kJ/kg] from P [MPa]
    and s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kg via Uniful.jl.

    Note: Do not use this function to attempt to find values on the saturation
    line, as any rounding errors anywhere will result in a point on either side
    of the phase boundary with large differences in enthalpy. Use SatHL/SatHV
    for saturated phase enthalpies,
"""
function SpecificH_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificH, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificH, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificH, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificH, P, T)
    end
end


function SpecificH_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificH, P.val, T)*u"kJ/kg"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificH, P.val, T)*u"kJ/kg"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificH, P.val, T)*u"kJ/kg"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificH, P.val, T)*u"kJ/kg"
    end
end


"""
    SpecificCP

    Utility function that returns the specific isobaric heat capacity [kJ/kgK]
    from P [MPa] and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificCP(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificCP, P, T)
    elseif Region == :Region2
        return Region2(:SpecificCP, P, T)
    elseif Region == :Region3
        return Region3(:SpecificCP, P, T)
    elseif Region == :Region5
        return Region5(:SpecificCP, P, T)
    end
end


function SpecificCP(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificCP, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region2
        return Region2(:SpecificCP, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region3
        return Region3(:SpecificCP, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region5
        return Region5(:SpecificCP, P.val, T.val)*u"kJ/kg/K"
    end
end


"""
    SpecificCP_Ph

    Utility function that returns the specific entropy [kJ/kgK] from P [MPa]
    and h [kJ/kg].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificCP_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificCP, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificCP, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificCP, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificCP, P, T)
    end
end


function SpecificCP_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificCP, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificCP, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificCP, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificCP, P.val, T)*u"kJ/kg/K"
    end
end


"""
    SpecificCP_Ps

    Utility function that returns the speciic entropy [kJ/kgK] from P [MPa]
    and s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificCP_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificCP, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificCP, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificCP, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificCP, P, T)
    end
end


function SpecificCP_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificCP, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificCP, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificCP, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificCP, P.val, T)*u"kJ/kg/K"
    end
end


"""
    SpecificCV

    Utility function that returns the specific isochoric heat capacity [kJ/kgK]
    from P [MPa] and T [K].
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificCV(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpecificCV, P, T)
    elseif Region == :Region2
        return Region2(:SpecificCV, P, T)
    elseif Region == :Region3
        return Region3(:SpecificCV, P, T)
    elseif Region == :Region5
        return Region5(:SpecificCV, P, T)
    end
end


function SpecificCV(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpecificCV, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region2
        return Region2(:SpecificCV, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region3
        return Region3(:SpecificCV, P.val, T.val)*u"kJ/kg/K"
    elseif Region == :Region5
        return Region5(:SpecificCV, P.val, T.val)*u"kJ/kg/K"
    end
end


"""
    SpecificCV_Ph

    Utility function that returns the specific entropy [kJ/kgK] from P [MPa]
    and h [kJ/kg].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificCV_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpecificCV, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpecificCV, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpecificCV, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpecificCV, P, T)
    end
end


function SpecificCV_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpecificCV, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpecificCV, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpecificCV, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpecificCV, P.val, T)*u"kJ/kg/K"
    end
end


"""
    SpecificCV_Ps

    Utility function that returns the speciic entropy [kJ/kgK] from P [MPa]
    and s [kJ/kgK].
    The explicit backwards equations are only available in regions 1 and 2.
    Input outside these will result in a DomainError exception.
    If inputs have associated units, the value is returned with associated
    units of kJ/kgK via Uniful.jl.
"""
function SpecificCV_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpecificCV, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpecificCV, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpecificCV, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpecificCV, P, T)
    end
end


function SpecificCV_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpecificCV, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpecificCV, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpecificCV, P.val, T)*u"kJ/kg/K"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpecificCV, P.val, T)*u"kJ/kg/K"
    end
end


"""
    SpeedOfSound

    Utility function that returns the sonic velocity [m/s] from P [MPa] and T [K].
    If inputs have associated units, the value is returned with associated
    units of m/s via Uniful.jl.
"""
function SpeedOfSound(P, T)
    Region = RegionID(P, T)

    if Region == :Region1
        return Region1(:SpeedOfSound, P, T)
    elseif Region == :Region2
        return Region2(:SpeedOfSound, P, T)
    elseif Region == :Region3
        return Region3(:SpeedOfSound, P, T)
    elseif Region == :Region5
        return Region5(:SpeedOfSound, P, T)
    end
end


function SpeedOfSound(P::Q1, T::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    Region = RegionID(P.val, T.val)

    if Region == :Region1
        return Region1(:SpeedOfSound, P.val, T.val)*u"m/s"
    elseif Region == :Region2
        return Region2(:SpeedOfSound, P.val, T.val)*u"m/s"
    elseif Region == :Region3
        return Region3(:SpeedOfSound, P.val, T.val)*u"m/s"
    elseif Region == :Region5
        return Region5(:SpeedOfSound, P.val, T.val)*u"m/s"
    end
end


"""
    SpeedOfSound_Ph

    Utility function that returns the sonic velocity [m/s] from P [MPa] and h [kJ/kg].
    If inputs have associated units, the value is returned with associated
    units of m/s via Uniful.jl."""
function SpeedOfSound_Ph(P, h)
    Region = RegionID_Ph(P, h)

    if Region == :Region1
        T = Region1_TPh(P, h)
        return Region1(:SpeedOfSound, P, T)
    elseif Region == :Region2a
        T = Region2a_TPh(P, h)
        return Region2(:SpeedOfSound, P, T)
    elseif Region == :Region2b
        T = Region2b_TPh(P, h)
        return Region2(:SpeedOfSound, P, T)
    elseif Region == :Region2c
        T = Region2c_TPh(P, h)
        return Region2(:SpeedOfSound, P, T)
    end
end


function SpeedOfSound_Ph(P::Q1, h::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        h = 1.0*uconvert(u"kJ/kg", h)
    catch
        throw(UnitsError(h, "Invalid input units."))
    end

    Region = RegionID_Ph(P.val, h.val)

    if Region == :Region1
        T = Region1_TPh(P.val, h.val)
        return Region1(:SpeedOfSound, P.val, T)*u"m/s"
    elseif Region == :Region2a
        T = Region2a_TPh(P.val, h.val)
        return Region2(:SpeedOfSound, P.val, T)*u"m/s"
    elseif Region == :Region2b
        T = Region2b_TPh(P.val, h.val)
        return Region2(:SpeedOfSound, P.val, T)*u"m/s"
    elseif Region == :Region2c
        T = Region2c_TPh(P.val, h.val)
        return Region2(:SpeedOfSound, P.val, T)*u"m/s"
    end
end


"""
    SpeedOfSound_Ps

    Utility function that returns the sonic velocity [m/s] from P [MPa] and s [kJ/kgK].
    If inputs have associated units, the value is returned with associated
    units of m/s via Uniful.jl.
"""
function SpeedOfSound_Ps(P, s)
    Region = RegionID_Ps(P, s)

    if Region == :Region1
        T = Region1_TPs(P, s)
        return Region1(:SpeedOfSound, P, T)
    elseif Region == :Region2a
        T = Region2a_TPs(P, s)
        return Region2(:SpeedOfSound, P, T)
    elseif Region == :Region2b
        T = Region2b_TPs(P, s)
        return Region2(:SpeedOfSound, P, T)
    elseif Region == :Region2c
        T = Region2c_TPs(P, s)
        return Region2(:SpeedOfSound, P, T)
    end
end


function SpeedOfSound_Ps(P::Q1, s::Q2) where Q1 <: Quantity where Q2 <: Quantity
    try
        P = 1.0*uconvert(u"MPa", P)
    catch
        throw(UnitsError(P, "Invalid input units."))
    end
    try
        s = 1.0*uconvert(u"kJ/kg/K", s)
    catch
        throw(UnitsError(s, "Invalid input units."))
    end

    Region = RegionID_Ps(P.val, s.val)

    if Region == :Region1
        T = Region1_TPs(P.val, s.val)
        return Region1(:SpeedOfSound, P.val, T)*u"m/s"
    elseif Region == :Region2a
        T = Region2a_TPs(P.val, s.val)
        return Region2(:SpeedOfSound, P.val, T)*u"m/s"
    elseif Region == :Region2b
        T = Region2b_TPs(P.val, s.val)
        return Region2(:SpeedOfSound, P.val, T)*u"m/s"
    elseif Region == :Region2c
        T = Region2c_TPs(P.val, s.val)
        return Region2(:SpeedOfSound, P.val, T)*u"m/s"
    end
end


"""
    SatDensL

    Returns the saturated liquid density [kg/m3] from T [K].
    If inputs have associated units, the value is returned with associated
    units of kg/m3 via Uniful.jl.
"""
function SatDensL(T)
    if T3 ≤ T ≤ Tc
        b = [   1.992_740_64,
                1.099_653_42,
               -0.510_839_303,
               -1.754_934_79,
              -45.517_035_2,
               -6.746_944_50E5]

        θ = T/Tc
        τ = 1 - θ

        return ρc*(1 + b[1]*τ^(1/3) + b[2]*τ^(2/3) + b[3]*τ^(5/3) + b[4]*τ^(16/3)
               + b[5]*τ^(43/3) + b[6]*τ^(110/3))
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end #SatDensL


function SatDensL(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        b = [   1.992_740_64,
                1.099_653_42,
               -0.510_839_303,
               -1.754_934_79,
              -45.517_035_2,
               -6.746_944_50E5]

        θ = T.val/Tc
        τ = 1 - θ

        return ρc*(1 + b[1]*τ^(1/3) + b[2]*τ^(2/3) + b[3]*τ^(5/3) + b[4]*τ^(16/3)
               + b[5]*τ^(43/3) + b[6]*τ^(110/3))*u"kg/m^3"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end #SatDensL


"""
    SatDensV

    Returns the saturated vapour density [kg/m3] from T [K].
    If inputs have associated units, the value is returned with associated
    units of kg/m3 via Uniful.jl.
"""
function SatDensV(T)
    if T3 ≤ T ≤ Tc
        c = [  -2.031_502_40,
               -2.683_029_40,
               -5.386_264_92,
              -17.299_160_5,
              -44.758_658_1,
              -63.920_106_3
        ]

        θ = T/Tc
        τ = 1 - θ

        return ρc*exp(c[1]*τ^(2/6) + c[2]*τ^(4/6) + c[3]*τ^(8/6) + c[4]*τ^(18/6)
             + c[5]*τ^(37/6) + c[6]*τ^(71/6))
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end #SatDensV


function SatDensV(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        c = [  -2.031_502_40,
               -2.683_029_40,
               -5.386_264_92,
              -17.299_160_5,
              -44.758_658_1,
              -63.920_106_3
        ]

        θ = T.val/Tc
        τ = 1 - θ

        return ρc*exp(c[1]*τ^(2/6) + c[2]*τ^(4/6) + c[3]*τ^(8/6) + c[4]*τ^(18/6)
             + c[5]*τ^(37/6) + c[6]*τ^(71/6))*u"kg/m^3"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end #SatDensV


"""
    SatHL

    Returns the saturated liquid specific enthalpy [J/kg] from T [K].
    If inputs have associated units, the value is returned with associated
    units of J/kg via Uniful.jl.
"""
function SatHL(T)
    if T3 ≤ T ≤ Tc
        α0 = 1000
        dα = -1135.905_627_715
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T/Tc
        α = α0*(dα + d[1]*θ^(-19) + d[2]*θ + d[3]*θ^4.5 + d[4]*θ^5 + d[5]*θ^54.5)
        return α + T/SatDensL(T)*1e6*ForwardDiff.derivative(Psat2, T)
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


function SatHL(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        α0 = 1000
        dα = -1135.905_627_715
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T.val/Tc
        α = α0*(dα + d[1]*θ^(-19) + d[2]*θ + d[3]*θ^4.5 + d[4]*θ^5 + d[5]*θ^54.5)
        return (α + T.val/SatDensL(T.val)*1e6*ForwardDiff.derivative(Psat2, T.val))*u"J/kg"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


"""
    SatHV

    Returns the saturated vapour specific enthalpy [J/kg] from T [K].
    If inputs have associated units, the value is returned with associated
    units of J/kg via Uniful.jl.
"""
function SatHV(T)
    if T3 ≤ T ≤ Tc
        α0 = 1000
        dα = -1135.905_627_715
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T/Tc
        α = α0*(dα + d[1]*θ^(-19) + d[2]*θ + d[3]*θ^4.5 + d[4]*θ^5 + d[5]*θ^54.5)
        return α + T/SatDensV(T)*1e6*ForwardDiff.derivative(Psat2, T)
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


function SatHV(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        α0 = 1000
        dα = -1135.905_627_715
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T.val/Tc
        α = α0*(dα + d[1]*θ^(-19) + d[2]*θ + d[3]*θ^4.5 + d[4]*θ^5 + d[5]*θ^54.5)
        return (α + T.val/SatDensV(T.val)*1e6*ForwardDiff.derivative(Psat2, T.val))*u"J/kg"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


"""
    SatSL

    Returns the saturated liquid specific entropy [J/kg] from T [K].
    If inputs have associated units, the value is returned with associated
    units of J/kgK via Uniful.jl.
"""
function SatSL(T)
    if T3 ≤ T ≤ Tc
        α0 = 1000
        ϕ0 = α0/Tc

        dϕ = 2319.5246
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T/Tc
        ϕ = ϕ0*(dϕ + 19/20*d[1]*θ^(-20) + d[2]*log(θ) + 9/7*d[3]*θ^3.5
                   + 5/4*d[4]*θ^4 + 109/107*d[5]*θ^53.5)
        return ϕ + 1/SatDensL(T)*1e6*ForwardDiff.derivative(Psat2, T)
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


function SatSL(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        α0 = 1000
        ϕ0 = α0/Tc

        dϕ = 2319.5246
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T.val/Tc
        ϕ = ϕ0*(dϕ + 19/20*d[1]*θ^(-20) + d[2]*log(θ) + 9/7*d[3]*θ^3.5
                   + 5/4*d[4]*θ^4 + 109/107*d[5]*θ^53.5)
        return (ϕ + 1/SatDensL(T.val)*1e6*ForwardDiff.derivative(Psat2, T.val))*u"J/kg"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


"""
    SatSV

    Returns the saturated liquid specific entropy [J/kg] from T [K].
    If inputs have associated units, the value is returned with associated
    units of J/kgK via Uniful.jl.
"""
function SatSV(T)
    if T3 ≤ T ≤ Tc
        α0 = 1000
        ϕ0 = α0/Tc

        dϕ = 2319.5246
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T/Tc
        ϕ = ϕ0*(dϕ + 19/20*d[1]*θ^(-20) + d[2]*log(θ) + 9/7*d[3]*θ^3.5
                   + 5/4*d[4]*θ^4 + 109/107*d[5]*θ^53.5)
        return ϕ + 1/SatDensV(T)*1e6*ForwardDiff.derivative(Psat2, T)
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


function SatSV(T::Q) where Q <: Quantity
    try
        T = 1.0*uconvert(u"K", T)
    catch
        throw(UnitsError(T, "Invalid input units."))
    end

    if T3 ≤ T.val ≤ Tc
        α0 = 1000
        ϕ0 = α0/Tc

        dϕ = 2319.5246
        d = [   -5.651_349_98E-8,
              2690.666_31,
               127.287_297,
              -135.003_439,
                 0.981_825_814
            ]

        θ = T.val/Tc
        ϕ = ϕ0*(dϕ + 19/20*d[1]*θ^(-20) + d[2]*log(θ) + 9/7*d[3]*θ^3.5
                   + 5/4*d[4]*θ^4 + 109/107*d[5]*θ^53.5)
        return (ϕ + 1/SatDensV(T.val)*1e6*ForwardDiff.derivative(Psat2, T.val))*u"J/kg"
    else
        throw(DomainError(T, "Temperature not between tripple and critical points"))
    end
end


"""
    DeltaHvap

    Returns the latent heat of vaporisation [J/kg] at T [K].
    If inputs have associated units, the value is returned with associated
    units of J/kg via Uniful.jl.
"""
function DeltaHvap(T)
    return SatHV(T) - SatHL(T)
end


"""
    Quality_Ph

    Returns vapour quality from P and h.
"""
function Quality_Ph(P, h)
    # Get the temperature from P, assuming saturated
    T = Tsat(P)
    hl = SatHL(T)
    hv = SatHV(T)

    return (h - hl)/(hv - hl)
end


"""
    Quality_Th

    Returns vapour quality from T and h.
"""
function Quality_Th(T, h)
    # Get the temperature from P, assuming saturated
    hl = SatHL(T)
    hv = SatHV(T)

    return (h - hl)/(hv - hl)
end


"""
    Quality_Ps

    Returns vapour quality from P and s.
"""
function Quality_Ps(P, s)
    # Get the temperature from P, assuming saturated
    T = Tsat(P)
    sl = SatSL(T)
    sv = SatSV(T)

    return (s - sl)/(sv - sl)
end


"""
    Quality_Ts

    Returns vapour quality from T and s.
"""
function Quality_Ts(T, s)
    # Get the temperature from P, assuming saturated
    sl = SatSL(T)
    sv = SatSV(T)

    return (s - sl)/(sv - sl)
end

@setup_workload begin
    @compile_workload begin
        include("compilefile.jl")
    end
end

end # module
