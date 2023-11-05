mysignif(x, n) = round(x; sigdigits = n)
function mysignif(x::Q, n::Int) where Q <: Quantity
    z = x.val
    y = Quantity(round(z, sigdigits = n), unit(x))
    return y
end
dummy = Bool[]

function runprecompworkload()
    push!(dummy, SpecificV(3.0, 300.0) ≈ 0.100_215_168E-2)
    push!(dummy, SpecificV(3.0u"MPa", 300.0u"K") ≈ 0.100_215_168E-2u"m^3/kg")

    push!(dummy, SpecificV(80.0, 300.0) ≈ 0.971_180_894E-3)
    push!(dummy, SpecificV(80.0u"MPa", 300.0u"K") ≈ 0.971_180_894E-3u"m^3/kg")

    push!(dummy, SpecificV(3.0, 500.0) ≈ 0.120_241_800E-2)
    push!(dummy, SpecificV(3.0u"MPa", 500.0u"K") ≈ 0.120_241_800E-2u"m^3/kg")

    push!(dummy, SpecificH(3.0, 300.0) ≈ 0.115_331_273E3)
    push!(dummy, SpecificH(3.0u"MPa", 300.0u"K") ≈ 0.115_331_273E3u"kJ/kg")

    push!(dummy, SpecificH(80.0, 300.0) ≈ 0.184_142_828E3)
    push!(dummy, SpecificH(80.0u"MPa", 300.0u"K") ≈ 0.184_142_828E3u"kJ/kg")

    push!(dummy, SpecificH(3.0, 500.0) ≈ 0.975_542_239E3)
    push!(dummy, SpecificH(3.0u"MPa", 500.0u"K") ≈ 0.975_542_239E3u"kJ/kg")

    push!(dummy, SpecificU(3.0, 300.0) ≈ 0.112_324_818E3)
    push!(dummy, SpecificU(3.0u"MPa", 300.0u"K") ≈ 0.112_324_818E3u"kJ/kg")

    push!(dummy, SpecificU(80.0, 300.0) ≈ 0.106_448_356E3)
    push!(dummy, SpecificU(80.0u"MPa", 300.0u"K") ≈ 0.106_448_356E3u"kJ/kg")

    push!(dummy, SpecificU(3.0, 500.0) ≈ 0.971_934_985E3)
    push!(dummy, SpecificU(3.0u"MPa", 500.0u"K") ≈ 0.971_934_985E3u"kJ/kg")

    push!(dummy, SpecificS(3.0, 300.0) ≈ 0.392_294_792)
    push!(dummy, SpecificS(3.0u"MPa", 300.0u"K") ≈ 0.392_294_792u"kJ/kg/K")

    push!(dummy, SpecificS(80.0, 300.0) ≈ 0.368_563_852)
    push!(dummy, SpecificS(80.0u"MPa", 300.0u"K") ≈ 0.368_563_852u"kJ/kg/K")

    push!(dummy, SpecificS(3.0, 500.0) ≈ 0.258_041_912E1)
    push!(dummy, SpecificS(3.0u"MPa", 500.0u"K") ≈ 0.258_041_912E1u"kJ/kg/K")

    push!(dummy, SpecificCP(3.0, 300.0) ≈ 0.417_301_218E1)
    push!(dummy, SpecificCP(3.0u"MPa", 300.0u"K") ≈ 0.417_301_218E1u"kJ/kg/K")

    push!(dummy, SpecificCP(80.0, 300.0) ≈ 0.401_008_987E1)
    push!(dummy, SpecificCP(80.0u"MPa", 300.0u"K") ≈ 0.401_008_987E1u"kJ/kg/K")

    push!(dummy, SpecificCP(3.0, 500.0) ≈ 0.465_580_682E1)
    push!(dummy, SpecificCP(3.0u"MPa", 500.0u"K") ≈ 0.465_580_682E1u"kJ/kg/K")

    push!(dummy, SpeedOfSound(3.0, 300.0) ≈ 0.150_773_921E4)
    push!(dummy, SpeedOfSound(3.0u"MPa", 300.0u"K") ≈ 0.150_773_921E4u"m/s")

    push!(dummy, SpeedOfSound(80.0, 300.0) ≈ 0.163_469_054E4)
    push!(dummy, SpeedOfSound(80.0u"MPa", 300.0u"K") ≈ 0.163_469_054E4u"m/s")

    push!(dummy, SpeedOfSound(3.0, 500.0) ≈ 0.124_071_337E4)
    push!(dummy, SpeedOfSound(3.0u"MPa", 500.0u"K") ≈ 0.124_071_337E4u"m/s")

    push!(dummy, mysignif(SpecificG(3.0, 0.391_798_509E3), 6) ≈ mysignif(SpecificG_Ph(3.0, 500.0), 6))
    push!(dummy, mysignif(SpecificG(3.0u"MPa", 0.391_798_509E3u"K"), 6) ≈ mysignif(SpecificG_Ph(3.0u"MPa", 500.0u"kJ/kg"), 6))

    push!(dummy, mysignif(SpecificG(80.0, 0.378_108_626E3), 6) ≈ mysignif(SpecificG_Ph(80.0, 500.0), 6))
    push!(dummy, mysignif(SpecificG(80.0u"MPa", 0.378_108_626E3u"K"), 6) ≈ mysignif(SpecificG_Ph(80.0u"MPa", 500.0u"kJ/kg"), 6))

    push!(dummy, SpecificG(80.0, 0.611_041_229E3) ≈ SpecificG_Ph(80.0, 1500.0))
    push!(dummy, SpecificG(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificG_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpecificF(3.0, 0.391_798_509E3) ≈ SpecificF_Ph(3.0, 500.0))
    push!(dummy, SpecificF(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificF_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, mysignif(SpecificF(80.0, 0.378_108_626E3), 6) ≈ mysignif(SpecificF_Ph(80.0, 500.0), 6))
    push!(dummy, mysignif(SpecificF(80.0u"MPa", 0.378_108_626E3u"K"), 6) ≈ mysignif(SpecificF_Ph(80.0u"MPa", 500.0u"kJ/kg"), 6))

    push!(dummy, SpecificF(80.0, 0.611_041_229E3) ≈ SpecificF_Ph(80.0, 1500.0))
    push!(dummy, SpecificF(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificF_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpecificV(3.0, 0.391_798_509E3) ≈ SpecificV_Ph(3.0, 500.0))
    push!(dummy, SpecificV(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificV_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificV(80.0, 0.378_108_626E3) ≈ SpecificV_Ph(80.0, 500.0))
    push!(dummy, SpecificV(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificV_Ph(80.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificV(80.0, 0.611_041_229E3) ≈ SpecificV_Ph(80.0, 1500.0))
    push!(dummy, SpecificV(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificV_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpecificU(3.0, 0.391_798_509E3) ≈ SpecificU_Ph(3.0, 500.0))
    push!(dummy, SpecificU(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificU_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificU(80.0, 0.378_108_626E3) ≈ SpecificU_Ph(80.0, 500.0))
    push!(dummy, SpecificU(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificU_Ph(80.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificU(80.0, 0.611_041_229E3) ≈ SpecificU_Ph(80.0, 1500.0))
    push!(dummy, SpecificU(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificU_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpecificS(3.0, 0.391_798_509E3) ≈ SpecificS_Ph(3.0, 500.0))
    push!(dummy, SpecificS(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificS_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificS(80.0, 0.378_108_626E3) ≈ SpecificS_Ph(80.0, 500.0))
    push!(dummy, SpecificS(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificS_Ph(80.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificS(80.0, 0.611_041_229E3) ≈ SpecificS_Ph(80.0, 1500.0))
    push!(dummy, SpecificS(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificS_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpecificCP(3.0, 0.391_798_509E3) ≈ SpecificCP_Ph(3.0, 500.0))
    push!(dummy, SpecificCP(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificCP_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificCP(80.0, 0.378_108_626E3) ≈ SpecificCP_Ph(80.0, 500.0))
    push!(dummy, SpecificCP(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificCP_Ph(80.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificCP(80.0, 0.611_041_229E3) ≈ SpecificCP_Ph(80.0, 1500.0))
    push!(dummy, SpecificCP(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificCP_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpecificCV(3.0, 0.391_798_509E3) ≈ SpecificCV_Ph(3.0, 500.0))
    push!(dummy, SpecificCV(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificCV_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificCV(80.0, 0.378_108_626E3) ≈ SpecificCV_Ph(80.0, 500.0))
    push!(dummy, SpecificCV(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificCV_Ph(80.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpecificCV(80.0, 0.611_041_229E3) ≈ SpecificCV_Ph(80.0, 1500.0))
    push!(dummy, SpecificCV(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificCV_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(3.0, 0.391_798_509E3) ≈ SpeedOfSound_Ph(3.0, 500.0))
    push!(dummy, SpeedOfSound(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(80.0, 0.378_108_626E3) ≈ SpeedOfSound_Ph(80.0, 500.0))
    push!(dummy, SpeedOfSound(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpeedOfSound_Ph(80.0u"MPa", 500.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(80.0, 0.611_041_229E3) ≈ SpeedOfSound_Ph(80.0, 1500.0))
    push!(dummy, SpeedOfSound(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpeedOfSound_Ph(80.0u"MPa", 1500.0u"kJ/kg"))

    push!(dummy, mysignif(SpecificG(3.0, 0.307_842_258E3), 6) ≈ mysignif(SpecificG_Ps(3.0, 0.5), 6))
    push!(dummy, mysignif(SpecificG(3.0u"MPa", 0.307_842_258E3u"K"), 6) ≈ mysignif(SpecificG_Ps(3.0u"MPa", 0.5u"kJ/kg/K"), 6))

    push!(dummy, mysignif(SpecificG(80.0, 0.309_979_785E3), 6) ≈ mysignif(SpecificG_Ps(80.0, 0.5), 6))
    push!(dummy, mysignif(SpecificG(80.0u"MPa", 0.309_979_785E3u"K"), 6) ≈ mysignif(SpecificG_Ps(80.0u"MPa", 0.5u"kJ/kg/K"), 6))

    push!(dummy, SpecificG(80.0, 0.565_899_909E3) ≈ SpecificG_Ps(80.0, 3.0))
    push!(dummy, SpecificG(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificG_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, mysignif(SpecificF(3.0, 0.307_842_258E3), 6) ≈ mysignif(SpecificF_Ps(3.0, 0.5), 6))
    push!(dummy, mysignif(SpecificF(3.0u"MPa", 0.307_842_258E3u"K"), 6) ≈ mysignif(SpecificF_Ps(3.0u"MPa", 0.5u"kJ/kg/K"), 6))

    push!(dummy, SpecificF(80.0, 0.309_979_785E3) ≈ SpecificF_Ps(80.0, 0.5))
    push!(dummy, SpecificF(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificF_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificF(80.0, 0.565_899_909E3) ≈ SpecificF_Ps(80.0, 3.0))
    push!(dummy, SpecificF(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificF_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpecificV(3.0, 0.307_842_258E3) ≈ SpecificV_Ps(3.0, 0.5))
    push!(dummy, SpecificV(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificV_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificV(80.0, 0.309_979_785E3) ≈ SpecificV_Ps(80.0, 0.5))
    push!(dummy, SpecificV(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificV_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificV(80.0, 0.565_899_909E3) ≈ SpecificV_Ps(80.0, 3.0))
    push!(dummy, SpecificV(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificV_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpecificU(3.0, 0.307_842_258E3) ≈ SpecificU_Ps(3.0, 0.5))
    push!(dummy, SpecificU(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificU_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificU(80.0, 0.309_979_785E3) ≈ SpecificU_Ps(80.0, 0.5))
    push!(dummy, SpecificU(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificU_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificU(80.0, 0.565_899_909E3) ≈ SpecificU_Ps(80.0, 3.0))
    push!(dummy, SpecificU(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificU_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpecificH(3.0, 0.307_842_258E3) ≈ SpecificH_Ps(3.0, 0.5))
    push!(dummy, SpecificH(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificH_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificH(80.0, 0.309_979_785E3) ≈ SpecificH_Ps(80.0, 0.5))
    push!(dummy, SpecificH(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificH_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificH(80.0, 0.565_899_909E3) ≈ SpecificH_Ps(80.0, 3.0))
    push!(dummy, SpecificH(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificH_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpecificCP(3.0, 0.307_842_258E3) ≈ SpecificCP_Ps(3.0, 0.5))
    push!(dummy, SpecificCP(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificCP_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificCP(80.0, 0.309_979_785E3) ≈ SpecificCP_Ps(80.0, 0.5))
    push!(dummy, SpecificCP(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificCP_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificCP(80.0, 0.565_899_909E3) ≈ SpecificCP_Ps(80.0, 3.0))
    push!(dummy, SpecificCP(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificCP_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpecificCV(3.0, 0.307_842_258E3) ≈ SpecificCV_Ps(3.0, 0.5))
    push!(dummy, SpecificCV(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificCV_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificCV(80.0, 0.309_979_785E3) ≈ SpecificCV_Ps(80.0, 0.5))
    push!(dummy, SpecificCV(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificCV_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpecificCV(80.0,0.565_899_909E3) ≈ SpecificCV_Ps(80.0, 3.0))
    push!(dummy, SpecificCV(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificCV_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(3.0, 0.307_842_258E3) ≈ SpeedOfSound_Ps(3.0, 0.5))
    push!(dummy, SpeedOfSound(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpeedOfSound_Ps(3.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(80.0, 0.309_979_785E3) ≈ SpeedOfSound_Ps(80.0, 0.5))
    push!(dummy, SpeedOfSound(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpeedOfSound_Ps(80.0u"MPa", 0.5u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(80.0, 0.565_899_909E3) ≈ SpeedOfSound_Ps(80.0, 3.0))
    push!(dummy, SpeedOfSound(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpeedOfSound_Ps(80.0u"MPa", 3.0u"kJ/kg/K"))

    push!(dummy, SpecificV(0.0035, 300.0) ≈ 0.394_913_866E2)
    push!(dummy, SpecificV(0.0035u"MPa", 300.0u"K") ≈ 0.394_913_866E2u"m^3/kg")

    push!(dummy, SpecificV(0.0035, 700.0) ≈ 0.923_015_898E2)
    push!(dummy, SpecificV(0.0035u"MPa", 700.0u"K") ≈ 0.923_015_898E2u"m^3/kg")

    push!(dummy, SpecificV(30.0, 700.0) ≈ 0.542_946_619E-2)
    push!(dummy, SpecificV(30.0u"MPa", 700.0u"K") ≈ 0.542_946_619E-2u"m^3/kg")

    push!(dummy, SpecificH(0.0035, 300.0) ≈ 0.254_991_145E4)
    push!(dummy, SpecificH(0.0035u"MPa", 300.0u"K") ≈ 0.254_991_145E4u"kJ/kg")

    push!(dummy, SpecificH(0.0035, 700.0) ≈ 0.333_568_375E4)
    push!(dummy, SpecificH(0.0035u"MPa", 700.0u"K") ≈ 0.333_568_375E4u"kJ/kg")

    push!(dummy, SpecificH(30.0, 700.0) ≈ 0.263_149_474E4)
    push!(dummy, SpecificH(30.0u"MPa", 700.0u"K") ≈ 0.263_149_474E4u"kJ/kg")

    push!(dummy, SpecificU(0.0035, 300.0) ≈ 0.241_169_160E4)
    push!(dummy, SpecificU(0.0035u"MPa", 300.0u"K") ≈ 0.241_169_160E4u"kJ/kg")

    push!(dummy, SpecificU(0.0035, 700.0) ≈ 0.301_262_819E4)
    push!(dummy, SpecificU(0.0035u"MPa", 700.0u"K") ≈ 0.301_262_819E4u"kJ/kg")

    push!(dummy, SpecificU(30.0, 700.0) ≈ 0.246_861_076E4)
    push!(dummy, SpecificU(30.0u"MPa", 700.0u"K") ≈ 0.246_861_076E4u"kJ/kg")

    push!(dummy, SpecificS(0.0035, 300.0) ≈ 0.852_238_967E1)
    push!(dummy, SpecificS(0.0035u"MPa", 300.0u"K") ≈ 0.852_238_967E1u"kJ/kg/K")

    push!(dummy, SpecificS(0.0035, 700.0) ≈ 0.101_749_996E2)
    push!(dummy, SpecificS(0.0035u"MPa", 700.0u"K") ≈ 0.101_749_996E2u"kJ/kg/K")

    push!(dummy, SpecificS(30.0, 700.0) ≈ 0.517_540_298E1)
    push!(dummy, SpecificS(30.0u"MPa", 700.0u"K") ≈ 0.517_540_298E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.0035, 300.0) ≈ 0.191_300_162E1)
    push!(dummy, SpecificCP(0.0035u"MPa", 300.0u"K") ≈ 0.191_300_162E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.0035, 700.0) ≈ 0.208_141_274E1)
    push!(dummy, SpecificCP(0.0035u"MPa", 700.0u"K") ≈ 0.208_141_274E1u"kJ/kg/K")

    push!(dummy, SpecificCP(30.0, 700.0) ≈ 0.103_505_092E2)
    push!(dummy, SpecificCP(30.0u"MPa", 700.0u"K") ≈ 0.103_505_092E2u"kJ/kg/K")

    push!(dummy, SpeedOfSound(0.0035, 300.0) ≈ 0.427_920_172E3)
    push!(dummy, SpeedOfSound(0.0035u"MPa", 300.0u"K") ≈ 0.427_920_172E3u"m/s")

    push!(dummy, SpeedOfSound(0.0035, 700.0) ≈ 0.644_289_068E3)
    push!(dummy, SpeedOfSound(0.0035u"MPa", 700.0u"K") ≈ 0.644_289_068E3u"m/s")

    push!(dummy, SpeedOfSound(30.0, 700.0) ≈ 0.480_386_523E3)
    push!(dummy, SpeedOfSound(30.0u"MPa", 700.0u"K") ≈ 0.480_386_523E3u"m/s")

    push!(dummy, SpecificG(0.001, 0.534_433_241E3) ≈ SpecificG_Ph(0.001, 3000))
    push!(dummy, SpecificG(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificG_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificG(3.0, 0.575_373_370E3) ≈ SpecificG_Ph(3.0, 3000))
    push!(dummy, SpecificG(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificG_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificG(3.0, 0.101_077_577E4) ≈ SpecificG_Ph(3.0, 4000))
    push!(dummy, SpecificG(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificG_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificF(0.001, 0.534_433_241E3) ≈ SpecificF_Ph(0.001, 3000))
    push!(dummy, SpecificF(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificF_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificF(3.0, 0.575_373_370E3) ≈ SpecificF_Ph(3.0, 3000))
    push!(dummy, SpecificF(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificF_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificF(3.0, 0.101_077_577E4) ≈ SpecificF_Ph(3.0, 4000))
    push!(dummy, SpecificF(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificF_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificV(0.001, 0.534_433_241E3) ≈ SpecificV_Ph(0.001, 3000))
    push!(dummy, SpecificV(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificV_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificV(3.0, 0.575_373_370E3) ≈ SpecificV_Ph(3.0, 3000))
    push!(dummy, SpecificV(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificV_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificV(3.0, 0.101_077_577E4) ≈ SpecificV_Ph(3.0, 4000))
    push!(dummy, SpecificV(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificV_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificU(0.001, 0.534_433_241E3) ≈ SpecificU_Ph(0.001, 3000))
    push!(dummy, SpecificU(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificU_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificU(3.0, 0.575_373_370E3) ≈ SpecificU_Ph(3.0, 3000))
    push!(dummy, SpecificU(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificU_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificU(3.0, 0.101_077_577E4) ≈ SpecificU_Ph(3.0, 4000))
    push!(dummy, SpecificU(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificU_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificS(0.001, 0.534_433_241E3) ≈ SpecificS_Ph(0.001, 3000))
    push!(dummy, SpecificS(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificS_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificS(3.0, 0.575_373_370E3) ≈ SpecificS_Ph(3.0, 3000))
    push!(dummy, SpecificS(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificS_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificS(3.0, 0.101_077_577E4) ≈ SpecificS_Ph(3.0, 4000))
    push!(dummy, SpecificS(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificS_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificCP(0.001, 0.534_433_241E3) ≈ SpecificCP_Ph(0.001, 3000))
    push!(dummy, SpecificCP(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificCP_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificCP(3.0, 0.575_373_370E3) ≈ SpecificCP_Ph(3.0, 3000))
    push!(dummy, SpecificCP(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificCP_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificCP(3.0, 0.101_077_577E4) ≈ SpecificCP_Ph(3.0, 4000))
    push!(dummy, SpecificCP(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificCP_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificCV(0.001, 0.534_433_241E3) ≈ SpecificCV_Ph(0.001, 3000))
    push!(dummy, SpecificCV(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificCV_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificCV(3.0, 0.575_373_370E3) ≈ SpecificCV_Ph(3.0, 3000))
    push!(dummy, SpecificCV(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificCV_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpecificCV(3.0, 0.101_077_577E4) ≈ SpecificCV_Ph(3.0, 4000))
    push!(dummy, SpecificCV(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificCV_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpeedOfSound(0.001, 0.534_433_241E3) ≈ SpeedOfSound_Ph(0.001, 3000))
    push!(dummy, SpeedOfSound(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpeedOfSound_Ph(0.001u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpeedOfSound(3.0, 0.575_373_370E3) ≈ SpeedOfSound_Ph(3.0, 3000))
    push!(dummy, SpeedOfSound(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 3000u"kJ/kg"))

    push!(dummy, SpeedOfSound(3.0, 0.101_077_577E4) ≈ SpeedOfSound_Ph(3.0, 4000))
    push!(dummy, SpeedOfSound(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificG(5.0, 0.801_299_102E3) ≈ SpecificG_Ph(5.0, 3500))
    push!(dummy, SpecificG(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificG_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificG(5.0, 0.101_531_583E4) ≈ SpecificG_Ph(5.0, 4000))
    push!(dummy, SpecificG(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificG_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificG(25.0, 0.875_279_054E3) ≈ SpecificG_Ph(25.0, 3500))
    push!(dummy, SpecificG(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificG_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificF(5.0, 0.801_299_102E3) ≈ SpecificF_Ph(5.0, 3500))
    push!(dummy, SpecificF(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificF_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificF(5.0, 0.101_531_583E4) ≈ SpecificF_Ph(5.0, 4000))
    push!(dummy, SpecificF(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificF_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificF(25.0, 0.875_279_054E3) ≈ SpecificF_Ph(25.0, 3500))
    push!(dummy, SpecificF(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificF_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificV(5.0, 0.801_299_102E3) ≈ SpecificV_Ph(5.0, 3500))
    push!(dummy, SpecificV(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificV_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificV(5.0, 0.101_531_583E4) ≈ SpecificV_Ph(5.0, 4000))
    push!(dummy, SpecificV(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificV_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificV(25.0, 0.875_279_054E3) ≈ SpecificV_Ph(25.0, 3500))
    push!(dummy, SpecificV(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificV_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificU(5.0, 0.801_299_102E3) ≈ SpecificU_Ph(5.0, 3500))
    push!(dummy, SpecificU(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificU_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificU(5.0, 0.101_531_583E4) ≈ SpecificU_Ph(5.0, 4000))
    push!(dummy, SpecificU(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificU_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificU(25.0, 0.875_279_054E3) ≈ SpecificU_Ph(25.0, 3500))
    push!(dummy, SpecificU(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificU_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificS(5.0, 0.801_299_102E3) ≈ SpecificS_Ph(5.0, 3500))
    push!(dummy, SpecificS(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificS_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificS(5.0, 0.101_531_583E4) ≈ SpecificS_Ph(5.0, 4000))
    push!(dummy, SpecificS(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificS_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificS(25.0, 0.875_279_054E3) ≈ SpecificS_Ph(25.0, 3500))
    push!(dummy, SpecificS(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificS_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificCP(5.0, 0.801_299_102E3) ≈ SpecificCP_Ph(5.0, 3500))
    push!(dummy, SpecificCP(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificCP_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificCP(5.0, 0.101_531_583E4) ≈ SpecificCP_Ph(5.0, 4000))
    push!(dummy, SpecificCP(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificCP_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificCP(25.0, 0.875_279_054E3) ≈ SpecificCP_Ph(25.0, 3500))
    push!(dummy, SpecificCP(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificCP_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificCV(5.0, 0.801_299_102E3) ≈ SpecificCV_Ph(5.0, 3500))
    push!(dummy, SpecificCV(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificCV_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificCV(5.0, 0.101_531_583E4) ≈ SpecificCV_Ph(5.0, 4000))
    push!(dummy, SpecificCV(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificCV_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpecificCV(25.0, 0.875_279_054E3) ≈ SpecificCV_Ph(25.0, 3500))
    push!(dummy, SpecificCV(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificCV_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpeedOfSound(5.0, 0.801_299_102E3) ≈ SpeedOfSound_Ph(5.0, 3500))
    push!(dummy, SpeedOfSound(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpeedOfSound_Ph(5.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpeedOfSound(5.0, 0.101_531_583E4) ≈ SpeedOfSound_Ph(5.0, 4000))
    push!(dummy, SpeedOfSound(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpeedOfSound_Ph(5.0u"MPa", 4000u"kJ/kg"))

    push!(dummy, SpeedOfSound(25.0, 0.875_279_054E3) ≈ SpeedOfSound_Ph(25.0, 3500))
    push!(dummy, SpeedOfSound(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpeedOfSound_Ph(25.0u"MPa", 3500u"kJ/kg"))

    push!(dummy, SpecificG(8.0, 0.106_495_556E4) ≈ SpecificG_Ps(8.0, 7.5))
    push!(dummy, SpecificG(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificG_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificG(90.0, 0.103_801_126E4) ≈ SpecificG_Ps(90.0, 6.0))
    push!(dummy, SpecificG(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificG_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificF(8.0, 0.600_484_040E3) ≈ SpecificF_Ps(8.0, 6.0))
    push!(dummy, SpecificF(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificF_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificF(8.0, 0.106_495_556E4) ≈ SpecificF_Ps(8.0, 7.5))
    push!(dummy, SpecificF(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificF_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificF(90.0, 0.103_801_126E4) ≈ SpecificF_Ps(90.0, 6.0))
    push!(dummy, SpecificF(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificF_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificV(8.0, 0.600_484_040E3) ≈ SpecificV_Ps(8.0, 6.0))
    push!(dummy, SpecificV(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificV_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificV(8.0, 0.106_495_556E4) ≈ SpecificV_Ps(8.0, 7.5))
    push!(dummy, SpecificV(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificV_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificV(90.0, 0.103_801_126E4) ≈ SpecificV_Ps(90.0, 6.0))
    push!(dummy, SpecificV(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificV_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificU(8.0, 0.600_484_040E3) ≈ SpecificU_Ps(8.0, 6.0))
    push!(dummy, SpecificU(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificU_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificU(8.0, 0.106_495_556E4) ≈ SpecificU_Ps(8.0, 7.5))
    push!(dummy, SpecificU(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificU_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificU(90.0, 0.103_801_126E4) ≈ SpecificU_Ps(90.0, 6.0))
    push!(dummy, SpecificU(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificU_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificH(8.0, 0.600_484_040E3) ≈ SpecificH_Ps(8.0, 6.0))
    push!(dummy, SpecificH(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificH_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificH(8.0, 0.106_495_556E4) ≈ SpecificH_Ps(8.0, 7.5))
    push!(dummy, SpecificH(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificH_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificH(90.0, 0.103_801_126E4) ≈ SpecificH_Ps(90.0, 6.0))
    push!(dummy, SpecificH(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificH_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificCP(8.0, 0.600_484_040E3) ≈ SpecificCP_Ps(8.0, 6.0))
    push!(dummy, SpecificCP(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificCP_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificCP(8.0, 0.106_495_556E4) ≈ SpecificCP_Ps(8.0, 7.5))
    push!(dummy, SpecificCP(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificCP_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificCP(90.0, 0.103_801_126E4) ≈ SpecificCP_Ps(90.0, 6.0))
    push!(dummy, SpecificCP(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificCP_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificCV(8.0, 0.600_484_040E3) ≈ SpecificCV_Ps(8.0, 6.0))
    push!(dummy, SpecificCV(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpecificCV_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificCV(8.0, 0.106_495_556E4) ≈ SpecificCV_Ps(8.0, 7.5))
    push!(dummy, SpecificCV(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpecificCV_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpecificCV(90.0, 0.103_801_126E4) ≈ SpecificCV_Ps(90.0, 6.0))
    push!(dummy, SpecificCV(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpecificCV_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(8.0, 0.600_484_040E3) ≈ SpeedOfSound_Ps(8.0, 6.0))
    push!(dummy, SpeedOfSound(8.0u"MPa", 0.600_484_040E3u"K") ≈ SpeedOfSound_Ps(8.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(8.0, 0.106_495_556E4) ≈ SpeedOfSound_Ps(8.0, 7.5))
    push!(dummy, SpeedOfSound(8.0u"MPa", 0.106_495_556E4u"K") ≈ SpeedOfSound_Ps(8.0u"MPa", 7.5u"kJ/kg/K"))

    push!(dummy, SpeedOfSound(90.0, 0.103_801_126E4) ≈ SpeedOfSound_Ps(90.0, 6.0))
    push!(dummy, SpeedOfSound(90.0u"MPa", 0.103_801_126E4u"K") ≈ SpeedOfSound_Ps(90.0u"MPa", 6.0u"kJ/kg/K"))

    push!(dummy, SpecificG(40.0, 0.743_056_411E3) ≈ SpecificG_Ph(40.0, 2700.0))
    push!(dummy, SpecificG(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificG_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificG(60.0, 0.791_137_067E3) ≈ SpecificG_Ph(60.0, 2700.0))
    push!(dummy, SpecificG(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificG_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificG(60.0, 0.882_756_860E3) ≈ SpecificG_Ph(60.0, 3200.0))
    push!(dummy, SpecificG(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificG_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificF(40.0, 0.743_056_411E3) ≈ SpecificF_Ph(40.0, 2700.0))
    push!(dummy, SpecificF(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificF_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificF(60.0, 0.791_137_067E3) ≈ SpecificF_Ph(60.0, 2700.0))
    push!(dummy, SpecificF(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificF_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificF(60.0, 0.882_756_860E3) ≈ SpecificF_Ph(60.0, 3200.0))
    push!(dummy, SpecificF(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificF_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificV(40.0, 0.743_056_411E3) ≈ SpecificV_Ph(40.0, 2700.0))
    push!(dummy, SpecificV(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificV_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificV(60.0, 0.791_137_067E3) ≈ SpecificV_Ph(60.0, 2700.0))
    push!(dummy, SpecificV(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificV_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificV(60.0, 0.882_756_860E3) ≈ SpecificV_Ph(60.0, 3200.0))
    push!(dummy, SpecificV(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificV_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificU(40.0, 0.743_056_411E3) ≈ SpecificU_Ph(40.0, 2700.0))
    push!(dummy, SpecificU(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificU_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificU(60.0, 0.791_137_067E3) ≈ SpecificU_Ph(60.0, 2700.0))
    push!(dummy, SpecificU(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificU_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificU(60.0, 0.882_756_860E3) ≈ SpecificU_Ph(60.0, 3200.0))
    push!(dummy, SpecificU(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificU_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificS(40.0, 0.743_056_411E3) ≈ SpecificS_Ph(40.0, 2700.0))
    push!(dummy, SpecificS(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificS_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificS(60.0, 0.791_137_067E3) ≈ SpecificS_Ph(60.0, 2700.0))
    push!(dummy, SpecificS(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificS_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificS(60.0, 0.882_756_860E3) ≈ SpecificS_Ph(60.0, 3200.0))
    push!(dummy, SpecificS(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificS_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificCP(40.0, 0.743_056_411E3) ≈ SpecificCP_Ph(40.0, 2700.0))
    push!(dummy, SpecificCP(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificCP_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificCP(60.0, 0.791_137_067E3) ≈ SpecificCP_Ph(60.0, 2700.0))
    push!(dummy, SpecificCP(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificCP_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificCP(60.0, 0.882_756_860E3) ≈ SpecificCP_Ph(60.0, 3200.0))
    push!(dummy, SpecificCP(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificCP_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificCV(40.0, 0.743_056_411E3) ≈ SpecificCV_Ph(40.0, 2700.0))
    push!(dummy, SpecificCV(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificCV_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificCV(60.0, 0.791_137_067E3) ≈ SpecificCV_Ph(60.0, 2700.0))
    push!(dummy, SpecificCV(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificCV_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpecificCV(60.0, 0.882_756_860E3) ≈ SpecificCV_Ph(60.0, 3200.0))
    push!(dummy, SpecificCV(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificCV_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(40.0, 0.743_056_411E3) ≈ SpeedOfSound_Ph(40.0, 2700.0))
    push!(dummy, SpeedOfSound(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpeedOfSound_Ph(40.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(60.0, 0.791_137_067E3) ≈ SpeedOfSound_Ph(60.0, 2700.0))
    push!(dummy, SpeedOfSound(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpeedOfSound_Ph(60.0u"MPa", 2700.0u"kJ/kg"))

    push!(dummy, SpeedOfSound(60.0, 0.882_756_860E3) ≈ SpeedOfSound_Ph(60.0, 3200.0))
    push!(dummy, SpeedOfSound(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpeedOfSound_Ph(60.0u"MPa", 3200.0u"kJ/kg"))

    push!(dummy, SpecificV(0.255_837_018E2, 650.0) ≈ 0.002)
    push!(dummy, SpecificV(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.002u"m^3/kg")

    push!(dummy, mysignif(SpecificV(0.222_930_643E2, 650.0), 3) ≈ 0.005)
    push!(dummy, mysignif(SpecificV(0.222_930_643E2u"MPa", 650.0u"K"), 3) ≈ 0.005u"m^3/kg")

    push!(dummy, SpecificV(0.783_095_639E2, 750.0) ≈ 0.002)
    push!(dummy, SpecificV(0.783_095_639E2u"MPa", 750.0u"K") ≈ 0.002u"m^3/kg")

    push!(dummy, SpecificH(0.255_837_018E2, 650.0) ≈ 0.186_343_019E4)
    push!(dummy, SpecificH(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.186_343_019E4u"kJ/kg")

    push!(dummy, SpecificH(0.222_930_643E2, 650.0) ≈ 0.237_512_401E4)
    push!(dummy, SpecificH(0.222_930_643E2u"MPa", 650.0u"K") ≈ 0.237_512_401E4u"kJ/kg")

    push!(dummy, SpecificH(0.783_095_639E2, 750.0) ≈ 0.225_868_845E4)
    push!(dummy, SpecificH(0.783_095_639E2u"MPa", 750.0u"K") ≈ 0.225_868_845E4u"kJ/kg")

    push!(dummy, SpecificU(0.255_837_018E2, 650.0) ≈ 0.181_226_279E4)
    push!(dummy, SpecificU(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.181_226_279E4u"kJ/kg")

    push!(dummy, SpecificU(0.222_930_643E2, 650.0) ≈ 0.226_365_868E4)
    push!(dummy, SpecificU(0.222_930_643E2u"MPa", 650.0u"K") ≈ 0.226_365_868E4u"kJ/kg")

    push!(dummy, SpecificU(0.783_095_639E2, 750.0) ≈ 0.210_206_932E4)
    push!(dummy, SpecificU(0.783_095_639E2u"MPa", 750.0u"K") ≈ 0.210_206_932E4u"kJ/kg")

    push!(dummy, SpecificS(0.255_837_018E2, 650.0) ≈ 0.405_427_273E1)
    push!(dummy, SpecificS(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.405_427_273E1u"kJ/kg/K")

    push!(dummy, SpecificS(0.222_930_643E2, 650.0) ≈ 0.485_438_792E1)
    push!(dummy, SpecificS(0.222_930_643E2u"MPa", 650.0u"K") ≈ 0.485_438_792E1u"kJ/kg/K")

    push!(dummy, SpecificS(0.783_095_639E2, 750.0) ≈ 0.446_971_906E1)
    push!(dummy, SpecificS(0.783_095_639E2u"MPa", 750.0u"K") ≈ 0.446_971_906E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.255_837_018E2, 650.0) ≈ 0.138_935_717E2)
    push!(dummy, SpecificCP(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.138_935_717E2u"kJ/kg/K")

    push!(dummy, SpecificCP(0.222_930_643E2, 650.0) ≈ 0.446_579_373E2)
    push!(dummy, SpecificCP(0.222_930_643E2u"MPa", 650.0u"K") ≈ 0.446_579_373E2u"kJ/kg/K")

    push!(dummy, SpecificCP(0.783_095_639E2, 750.0) ≈ 0.634_165_359E1)
    push!(dummy, SpecificCP(0.783_095_639E2u"MPa", 750.0u"K") ≈ 0.634_165_359E1u"kJ/kg/K")

    push!(dummy, SpeedOfSound(0.255_837_018E2, 650.0) ≈ 0.502_005_554E3)
    push!(dummy, SpeedOfSound(0.255_837_018E2u"MPa", 650.0u"K") ≈ 0.502_005_554E3u"m/s")

    push!(dummy, SpeedOfSound(0.222_930_643E2, 650.0) ≈ 0.383_444_594E3)
    push!(dummy, SpeedOfSound(0.222_930_643E2u"MPa", 650.0u"K") ≈ 0.383_444_594E3u"m/s")

    push!(dummy, SpeedOfSound(0.783_095_639E2, 750.0) ≈ 0.760_696_041E3)
    push!(dummy, SpeedOfSound(0.783_095_639E2u"MPa", 750.0u"K") ≈ 0.760_696_041E3u"m/s")

    push!(dummy, Tsat(0.10) ≈ 0.372_755_919E3)
    push!(dummy, Tsat(0.10u"MPa") ≈ 0.372_755_919E3u"K")

    push!(dummy, Tsat(1.00) ≈ 0.453_035_632E3)
    push!(dummy, Tsat(1.00u"MPa") ≈ 0.453_035_632E3u"K")

    push!(dummy, Tsat(10.0) ≈ 0.584_149_488E3)
    push!(dummy, Tsat(10.0u"MPa") ≈ 0.584_149_488E3u"K")

    push!(dummy, Psat(300.0) ≈ 0.353_658_941E-2)
    push!(dummy, Psat(300.0u"K") ≈ 0.353_658_941E-2u"MPa")

    push!(dummy, Psat(500.0) ≈ 0.263_889_776E1)
    push!(dummy, Psat(500.0u"K") ≈ 0.263_889_776E1u"MPa")

    push!(dummy, Psat(600.0) ≈ 0.123_443_146E2)
    push!(dummy, Psat(600.0u"K") ≈ 0.123_443_146E2u"MPa")

    push!(dummy, SpecificV(0.5, 1500.0) ≈ 0.138_455_090E1)
    push!(dummy, SpecificV(0.5u"MPa", 1500.0u"K") ≈ 0.138_455_090E1u"m^3/kg")

    push!(dummy, SpecificV(30.0, 1500.0) ≈ 0.230_761_299E-1)
    push!(dummy, SpecificV(30.0u"MPa", 1500.0u"K") ≈ 0.230_761_299E-1u"m^3/kg")

    push!(dummy, SpecificV(30.0, 2000.0) ≈ 0.311_385_219E-1)
    push!(dummy, SpecificV(30.0u"MPa", 2000.0u"K") ≈ 0.311_385_219E-1u"m^3/kg")

    push!(dummy, SpecificH(0.5, 1500.0) ≈ 0.521_976_855E4)
    push!(dummy, SpecificH(0.5u"MPa", 1500.0u"K") ≈ 0.521_976_855E4u"kJ/kg")

    push!(dummy, SpecificH(30.0, 1500.0) ≈ 0.516_723_514E4)
    push!(dummy, SpecificH(30.0u"MPa", 1500.0u"K") ≈ 0.516_723_514E4u"kJ/kg")

    push!(dummy, SpecificH(30.0, 2000.0) ≈ 0.657_122_604E4)
    push!(dummy, SpecificH(30.0u"MPa", 2000.0u"K") ≈ 0.657_122_604E4u"kJ/kg")

    push!(dummy, SpecificU(0.5, 1500.0) ≈ 0.452_749_310E4)
    push!(dummy, SpecificU(0.5u"MPa", 1500.0u"K") ≈ 0.452_749_310E4u"kJ/kg")

    push!(dummy, SpecificU(30.0, 1500.0) ≈ 0.447_495_124E4)
    push!(dummy, SpecificU(30.0u"MPa", 1500.0u"K") ≈ 0.447_495_124E4u"kJ/kg")

    push!(dummy, SpecificU(30.0, 2000.0) ≈ 0.563_707_038E4)
    push!(dummy, SpecificU(30.0u"MPa", 2000.0u"K") ≈ 0.563_707_038E4u"kJ/kg")

    push!(dummy, SpecificS(0.5, 1500.0) ≈ 0.965_408_875E1)
    push!(dummy, SpecificS(0.5u"MPa", 1500.0u"K") ≈ 0.965_408_875E1u"kJ/kg/K")

    push!(dummy, SpecificS(30.0, 1500.0) ≈ 0.772_970_133E1)
    push!(dummy, SpecificS(30.0u"MPa", 1500.0u"K") ≈ 0.772_970_133E1u"kJ/kg/K")

    push!(dummy, SpecificS(30.0, 2000.0) ≈ 0.853_640_523E1)
    push!(dummy, SpecificS(30.0u"MPa", 2000.0u"K") ≈ 0.853_640_523E1u"kJ/kg/K")

    push!(dummy, SpecificCP(0.5, 1500.0) ≈ 0.261_609_445E1)
    push!(dummy, SpecificCP(0.5u"MPa", 1500.0u"K") ≈ 0.261_609_445E1u"kJ/kg/K")

    push!(dummy, SpecificCP(30.0, 1500.0) ≈ 0.272_724_317E1)
    push!(dummy, SpecificCP(30.0u"MPa", 1500.0u"K") ≈ 0.272_724_317E1u"kJ/kg/K")

    push!(dummy, SpecificCP(30.0, 2000.0) ≈ 0.288_569_882E1)
    push!(dummy, SpecificCP(30.0u"MPa", 2000.0u"K") ≈ 0.288_569_882E1u"kJ/kg/K")

    push!(dummy, SpeedOfSound(0.5, 1500.0) ≈ 0.917_068_690E3)
    push!(dummy, SpeedOfSound(0.5u"MPa", 1500.0u"K") ≈ 0.917_068_690E3u"m/s")

    push!(dummy, SpeedOfSound(30.0, 1500.0) ≈ 0.928_548_002E3)
    push!(dummy, SpeedOfSound(30.0u"MPa", 1500.0u"K") ≈ 0.928_548_002E3u"m/s")

    push!(dummy, SpeedOfSound(30.0, 2000.0) ≈ 0.106_736_948E4)
    push!(dummy, SpeedOfSound(30.0u"MPa", 2000.0u"K") ≈ 0.106_736_948E4u"m/s")

    return dummy
end