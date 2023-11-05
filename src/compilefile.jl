@inline mysignif(x, n) = round(x; sigdigits = n)

@inline function mysignif(x::Q, n::Int) where Q <: Quantity
    z = x.val
    y = Quantity(round(z, sigdigits = n), unit(x))
    return y
end

function runprecompworkload()
    dummy = Bool[]

    push!(dummy, SpecificV(3.0, 300.0) ≈ 0.100_215_168E-2)
    push!(dummy, SpecificV(3.0u"MPa", 300.0u"K") ≈ 0.100_215_168E-2u"m^3/kg")

    push!(dummy, SpecificH(3.0, 300.0) ≈ 0.115_331_273E3)
    push!(dummy, SpecificH(3.0u"MPa", 300.0u"K") ≈ 0.115_331_273E3u"kJ/kg")

    push!(dummy, SpecificU(3.0, 500.0) ≈ 0.971_934_985E3)
    push!(dummy, SpecificU(3.0u"MPa", 500.0u"K") ≈ 0.971_934_985E3u"kJ/kg")

    push!(dummy, SpecificS(3.0, 300.0) ≈ 0.392_294_792)
    push!(dummy, SpecificS(3.0u"MPa", 300.0u"K") ≈ 0.392_294_792u"kJ/kg/K")

    push!(dummy, SpecificCP(3.0, 300.0) ≈ 0.417_301_218E1)
    push!(dummy, SpecificCP(3.0u"MPa", 300.0u"K") ≈ 0.417_301_218E1u"kJ/kg/K")

    push!(dummy, SpeedOfSound(3.0, 300.0) ≈ 0.150_773_921E4)
    push!(dummy, SpeedOfSound(3.0u"MPa", 300.0u"K") ≈ 0.150_773_921E4u"m/s")

    push!(dummy, mysignif(SpecificG(3.0, 0.391_798_509E3), 6) ≈ mysignif(SpecificG_Ph(3.0, 500.0), 6))
    push!(dummy, mysignif(SpecificG(3.0u"MPa", 0.391_798_509E3u"K"), 6) ≈ mysignif(SpecificG_Ph(3.0u"MPa", 500.0u"kJ/kg"), 6))

    push!(dummy, SpecificF(3.0, 0.391_798_509E3) ≈ SpecificF_Ph(3.0, 500.0))
    push!(dummy, SpecificF(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificF_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificV(3.0, 0.391_798_509E3) ≈ SpecificV_Ph(3.0, 500.0))
    push!(dummy, SpecificV(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificV_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificU(3.0, 0.391_798_509E3) ≈ SpecificU_Ph(3.0, 500.0))
    push!(dummy, SpecificU(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificU_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificS(3.0, 0.391_798_509E3) ≈ SpecificS_Ph(3.0, 500.0))
    push!(dummy, SpecificS(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificS_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificCP(3.0, 0.391_798_509E3) ≈ SpecificCP_Ph(3.0, 500.0))
    push!(dummy, SpecificCP(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificCP_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificCV(3.0, 0.391_798_509E3) ≈ SpecificCV_Ph(3.0, 500.0))
    push!(dummy, SpecificCV(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificCV_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(3.0, 0.391_798_509E3) ≈ SpeedOfSound_Ph(3.0, 500.0))
    push!(dummy, SpeedOfSound(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, mysignif(SpecificG(3.0, 0.307_842_258E3), 6) ≈ mysignif(SpecificG_Ps(3.0, 0.5), 6))
    push!(dummy, mysignif(SpecificG(3.0u"MPa", 0.307_842_258E3u"K"), 6) ≈ mysignif(SpecificG_Ps(3.0u"MPa", 0.5u"kJ/kg/K"), 6))

    push!(dummy, mysignif(SpecificF(3.0, 0.307_842_258E3), 6) ≈ mysignif(SpecificF_Ps(3.0, 0.5), 6))
    push!(dummy, mysignif(SpecificF(3.0u"MPa", 0.307_842_258E3u"K"), 6) ≈ mysignif(SpecificF_Ps(3.0u"MPa", 0.5u"kJ/kg/K"), 6))

    push!(dummy, SpecificV(3.0, 0.307_842_258E3) ≈ SpecificV_Ps(3.0, 0.5))
    push!(dummy, SpecificV(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificV_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificU(3.0, 0.307_842_258E3) ≈ SpecificU_Ps(3.0, 0.5))
    push!(dummy, SpecificU(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificU_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificH(3.0, 0.307_842_258E3) ≈ SpecificH_Ps(3.0, 0.5))
    push!(dummy, SpecificH(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificH_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificCP(3.0, 0.307_842_258E3) ≈ SpecificCP_Ps(3.0, 0.5))
    push!(dummy, SpecificCP(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificCP_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(3.0, 0.307_842_258E3) ≈ SpeedOfSound_Ps(3.0, 0.5))
    push!(dummy, SpeedOfSound(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpeedOfSound_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificV(0.0035, 300.0) ≈ 0.394_913_866E2)
    push!(dummy, SpecificV(0.0035u"MPa", 300.0u"K") ≈ 0.394_913_866E2u"m^3/kg")

    push!(dummy, SpecificH(0.0035, 300.0) ≈ 0.254_991_145E4)
    push!(dummy, SpecificH(0.0035u"MPa", 300.0u"K") ≈ 0.254_991_145E4u"kJ/kg")

    push!(dummy, SpecificU(0.0035, 300.0) ≈ 0.241_169_160E4)
    push!(dummy, SpecificU(0.0035u"MPa", 300.0u"K") ≈ 0.241_169_160E4u"kJ/kg")

    push!(dummy, SpecificS(0.0035, 300.0) ≈ 0.852_238_967E1)
    push!(dummy, SpecificS(0.0035u"MPa", 300.0u"K") ≈ 0.852_238_967E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.0035, 300.0) ≈ 0.191_300_162E1)
    push!(dummy, SpecificCP(0.0035u"MPa", 300.0u"K") ≈ 0.191_300_162E1u"kJ/kg/K")

    push!(dummy, SpeedOfSound(0.0035, 300.0) ≈ 0.427_920_172E3)
    push!(dummy, SpeedOfSound(0.0035u"MPa", 300.0u"K") ≈ 0.427_920_172E3u"m/s")

    push!(dummy, SpecificG(0.001, 0.534_433_241E3) ≈ SpecificG_Ph(0.001, 3000))
    push!(dummy, SpecificG(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificG_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificF(0.001, 0.534_433_241E3) ≈ SpecificF_Ph(0.001, 3000))
    push!(dummy, SpecificF(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificF_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificV(0.001, 0.534_433_241E3) ≈ SpecificV_Ph(0.001, 3000))
    push!(dummy, SpecificV(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificV_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificU(0.001, 0.534_433_241E3) ≈ SpecificU_Ph(0.001, 3000))
    push!(dummy, SpecificU(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificU_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificS(0.001, 0.534_433_241E3) ≈ SpecificS_Ph(0.001, 3000))
    push!(dummy, SpecificS(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificS_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificCP(0.001, 0.534_433_241E3) ≈ SpecificCP_Ph(0.001, 3000))
    push!(dummy, SpecificCP(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificCP_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificCV(0.001, 0.534_433_241E3) ≈ SpecificCV_Ph(0.001, 3000))
    push!(dummy, SpecificCV(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificCV_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpeedOfSound(0.001, 0.534_433_241E3) ≈ SpeedOfSound_Ph(0.001, 3000))
    push!(dummy, SpeedOfSound(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpeedOfSound_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificG(5.0, 0.801_299_102E3) ≈ SpecificG_Ph(5.0, 3500))
    push!(dummy, SpecificG(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificG_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificF(5.0, 0.801_299_102E3) ≈ SpecificF_Ph(5.0, 3500))
    push!(dummy, SpecificF(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificF_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificV(5.0, 0.801_299_102E3) ≈ SpecificV_Ph(5.0, 3500))
    push!(dummy, SpecificV(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificV_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificU(5.0, 0.801_299_102E3) ≈ SpecificU_Ph(5.0, 3500))
    push!(dummy, SpecificU(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificU_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificS(5.0, 0.801_299_102E3) ≈ SpecificS_Ph(5.0, 3500))
    push!(dummy, SpecificS(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificS_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificCP(5.0, 0.801_299_102E3) ≈ SpecificCP_Ph(5.0, 3500))
    push!(dummy, SpecificCP(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificCP_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificCV(5.0, 0.801_299_102E3) ≈ SpecificCV_Ph(5.0, 3500))
    push!(dummy, SpecificCV(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificCV_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpeedOfSound(5.0, 0.801_299_102E3) ≈ SpeedOfSound_Ph(5.0, 3500))
    push!(dummy, SpeedOfSound(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpeedOfSound_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificG(8.0, 0.106_495_556E4) ≈ SpecificG_Ps(8.0, 7.5))
    push!(dummy, SpecificG(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificG_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificF(8.0, 0.600_484_040E3) ≈ SpecificF_Ps(8.0, 6.0))
    push!(dummy, SpecificF(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificF_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificV(8.0, 0.600_484_040E3) ≈ SpecificV_Ps(8.0, 6.0))
    push!(dummy, SpecificV(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificV_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificU(8.0, 0.600_484_040E3) ≈ SpecificU_Ps(8.0, 6.0))
    push!(dummy, SpecificU(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificU_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificH(8.0, 0.600_484_040E3) ≈ SpecificH_Ps(8.0, 6.0))
    push!(dummy, SpecificH(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificH_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificCP(8.0, 0.600_484_040E3) ≈ SpecificCP_Ps(8.0, 6.0))
    push!(dummy, SpecificCP(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificCP_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificCV(8.0, 0.600_484_040E3) ≈ SpecificCV_Ps(8.0, 6.0))
    push!(dummy, SpecificCV(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificCV_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(8.0, 0.600_484_040E3) ≈ SpeedOfSound_Ps(8.0, 6.0))
    push!(dummy, SpeedOfSound(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpeedOfSound_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificG(40.0, 0.743_056_411E3) ≈ SpecificG_Ph(40.0, 2700.0))
    push!(dummy, SpecificG(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificG_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificF(40.0, 0.743_056_411E3) ≈ SpecificF_Ph(40.0, 2700.0))
    push!(dummy, SpecificF(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificF_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificV(40.0, 0.743_056_411E3) ≈ SpecificV_Ph(40.0, 2700.0))
    push!(dummy, SpecificV(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificV_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificU(40.0, 0.743_056_411E3) ≈ SpecificU_Ph(40.0, 2700.0))
    push!(dummy, SpecificU(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificU_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificS(40.0, 0.743_056_411E3) ≈ SpecificS_Ph(40.0, 2700.0))
    push!(dummy, SpecificS(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificS_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificCP(40.0, 0.743_056_411E3) ≈ SpecificCP_Ph(40.0, 2700.0))
    push!(dummy, SpecificCP(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificCP_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificCV(40.0, 0.743_056_411E3) ≈ SpecificCV_Ph(40.0, 2700.0))
    push!(dummy, SpecificCV(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificCV_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(40.0, 0.743_056_411E3) ≈ SpeedOfSound_Ph(40.0, 2700.0))
    push!(dummy, SpeedOfSound(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpeedOfSound_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificV(0.255_837_018E2, 650.0) ≈ 0.002)
    push!(dummy, SpecificV(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.002u"m^3/kg")

    push!(dummy, SpecificH(0.255_837_018E2, 650.0) ≈ 0.186_343_019E4)
    push!(dummy, SpecificH(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.186_343_019E4u"kJ/kg")

    push!(dummy, SpecificU(0.255_837_018E2, 650.0) ≈ 0.181_226_279E4)
    push!(dummy, SpecificU(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.181_226_279E4u"kJ/kg")

    push!(dummy, SpecificS(0.255_837_018E2, 650.0) ≈ 0.405_427_273E1)
    push!(dummy, SpecificS(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.405_427_273E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.255_837_018E2, 650.0) ≈ 0.138_935_717E2)
    push!(dummy, SpecificCP(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.138_935_717E2u"kJ/kg/K")

    push!(dummy, SpeedOfSound(0.255_837_018E2, 650.0) ≈ 0.502_005_554E3)
    push!(dummy, SpeedOfSound(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.502_005_554E3u"m/s")

    push!(dummy, Tsat(0.10) ≈ 0.372_755_919E3)
    push!(dummy, Tsat(0.10u"MPa") ≈ 0.372_755_919E3u"K")

    push!(dummy, Psat(300.0) ≈ 0.353_658_941E-2)
    push!(dummy, Psat(300.0u"K") ≈ 0.353_658_941E-2u"MPa")

    push!(dummy, SpecificV(0.5, 1500.0) ≈ 0.138_455_090E1)
    push!(dummy, SpecificV(0.5u"MPa", 1500.0u"K") ≈ 0.138_455_090E1u"m^3/kg")

    push!(dummy, SpecificH(0.5, 1500.0) ≈ 0.521_976_855E4)
    push!(dummy, SpecificH(0.5u"MPa", 1500.0u"K") ≈ 0.521_976_855E4u"kJ/kg")

    push!(dummy, SpecificU(0.5, 1500.0) ≈ 0.452_749_310E4)
    push!(dummy, SpecificU(0.5u"MPa", 1500.0u"K") ≈ 0.452_749_310E4u"kJ/kg")

    push!(dummy, SpecificS(0.5, 1500.0) ≈ 0.965_408_875E1)
    push!(dummy, SpecificS(0.5u"MPa", 1500.0u"K") ≈ 0.965_408_875E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.5, 1500.0) ≈ 0.261_609_445E1)
    push!(dummy, SpecificCP(0.5u"MPa", 1500.0u"K") ≈ 0.261_609_445E1u"kJ/kg/K")

    push!(dummy, SpeedOfSound(0.5, 1500.0) ≈ 0.917_068_690E3)
    push!(dummy, SpeedOfSound(0.5u"MPa", 1500.0u"K") ≈ 0.917_068_690E3u"m/s")

    push!(dummy, mysignif(SatDensL(273.16), 6) ≈ 999.789)
    push!(dummy, mysignif(SatDensL(273.16u"K"), 6) ≈ 999.789u"kg/m^3")

    push!(dummy, mysignif(SatDensL(373.1243), 6) ≈ 958.365)
    push!(dummy, mysignif(SatDensL(373.1243u"K"), 6) ≈ 958.365u"kg/m^3")

    push!(dummy, mysignif(SatDensL(647.096), 3) ≈ 322)
    push!(dummy, mysignif(SatDensL(647.096u"K"), 3) ≈ 322u"kg/m^3")

    push!(dummy, mysignif(SatDensV(273.16), 6) ≈ 0.00485426)
    push!(dummy, mysignif(SatDensV(273.16u"K"), 6) ≈ 0.00485426u"kg/m^3")

    push!(dummy, mysignif(SatDensV(373.1243), 6) ≈ 0.597586)
    push!(dummy, mysignif(SatDensV(373.1243u"K"), 6) ≈ 0.597586u"kg/m^3")

    push!(dummy, mysignif(SatDensV(647.096), 3) ≈ 322)
    push!(dummy, mysignif(SatDensV(647.096u"K"), 3) ≈ 322u"kg/m^3")

    push!(dummy, mysignif(SatHL(273.16), 6) ≈ mysignif(0.000611786, 6))
    push!(dummy, mysignif(SatHL(273.16u"K"), 6) ≈ mysignif(0.000611786u"kJ/kg", 6))

    push!(dummy, mysignif(SatHL(373.1243), 5) ≈ mysignif(419.05, 5))
    push!(dummy, mysignif(SatHL(373.1243u"K"), 5) ≈ mysignif(419.05u"kJ/kg", 5))

    push!(dummy, mysignif(SatHL(647.096), 5) ≈ mysignif(2086.6, 5))
    push!(dummy, mysignif(SatHL(647.096u"K"), 5) ≈ mysignif(2086.6u"kJ/kg", 5))

    push!(dummy, mysignif(SatHV(273.16), 5) ≈ 2500.5)
    push!(dummy, mysignif(SatHV(273.16u"K"), 5) ≈ 2500.5u"kJ/kg")

    push!(dummy, mysignif(SatHV(373.1243), 5) ≈ 2675.7)
    push!(dummy, mysignif(SatHV(373.1243u"K"), 5) ≈ 2675.7u"kJ/kg")

    push!(dummy, mysignif(SatHV(647.096), 5) ≈ 2086.6)
    push!(dummy, mysignif(SatHV(647.096u"K"), 5) ≈ 2086.6u"kJ/kg")

    push!(dummy, abs(SatSL(273.16)) < 1e-8)
    push!(dummy, abs(SatSL(273.16u"K")) < 1e-8u"kJ/kg/K")

    push!(dummy, mysignif(SatSL(373.1243), 4) ≈ 1.307)
    push!(dummy, mysignif(SatSL(373.1243u"K"), 4) ≈ 1.307u"kJ/kg/K")

    push!(dummy, mysignif(SatSL(647.096), 4) ≈ 4.410)
    push!(dummy, mysignif(SatSL(647.096u"K"), 4) ≈ 4.410u"kJ/kg/K")

    push!(dummy, mysignif(SatSV(273.16), 4) ≈ 9.154)
    push!(dummy, mysignif(SatSV(273.16u"K"), 4) ≈ 9.154u"kJ/kg/K")

    push!(dummy, mysignif(SatSV(373.1243), 4) ≈ 7.355)
    push!(dummy, mysignif(SatSV(373.1243u"K"), 4) ≈ 7.355u"kJ/kg/K")

    push!(dummy, mysignif(SatSV(647.096), 4) ≈ 4.410)
    push!(dummy, mysignif(SatSV(647.096u"K"), 4) ≈ 4.410u"kJ/kg/K")

    return dummy
end