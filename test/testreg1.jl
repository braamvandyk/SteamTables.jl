@test SpecificV(3.0, 300.0) ≈ 0.100_215_168E-2
@test SpecificV(3.0u"MPa", 300.0u"K") ≈ 0.100_215_168E-2u"m^3/kg"

@test SpecificV(80.0, 300.0) ≈ 0.971_180_894E-3
@test SpecificV(80.0u"MPa", 300.0u"K") ≈ 0.971_180_894E-3u"m^3/kg"

@test SpecificV(3.0, 500.0) ≈ 0.120_241_800E-2
@test SpecificV(3.0u"MPa", 500.0u"K") ≈ 0.120_241_800E-2u"m^3/kg"

@test SpecificH(3.0, 300.0) ≈ 0.115_331_273E3
@test SpecificH(3.0u"MPa", 300.0u"K") ≈ 0.115_331_273E3u"kJ/kg"

@test SpecificH(80.0, 300.0) ≈ 0.184_142_828E3
@test SpecificH(80.0u"MPa", 300.0u"K") ≈ 0.184_142_828E3u"kJ/kg"

@test SpecificH(3.0, 500.0) ≈ 0.975_542_239E3
@test SpecificH(3.0u"MPa", 500.0u"K") ≈ 0.975_542_239E3u"kJ/kg"

@test SpecificU(3.0, 300.0) ≈ 0.112_324_818E3
@test SpecificU(3.0u"MPa", 300.0u"K") ≈ 0.112_324_818E3u"kJ/kg"

@test SpecificU(80.0, 300.0) ≈ 0.106_448_356E3
@test SpecificU(80.0u"MPa", 300.0u"K") ≈ 0.106_448_356E3u"kJ/kg"

@test SpecificU(3.0, 500.0) ≈ 0.971_934_985E3
@test SpecificU(3.0u"MPa", 500.0u"K") ≈ 0.971_934_985E3u"kJ/kg"

@test SpecificS(3.0, 300.0) ≈ 0.392_294_792
@test SpecificS(3.0u"MPa", 300.0u"K") ≈ 0.392_294_792u"kJ/kg/K"

@test SpecificS(80.0, 300.0) ≈ 0.368_563_852
@test SpecificS(80.0u"MPa", 300.0u"K") ≈ 0.368_563_852u"kJ/kg/K"

@test SpecificS(3.0, 500.0) ≈ 0.258_041_912E1
@test SpecificS(3.0u"MPa", 500.0u"K") ≈ 0.258_041_912E1u"kJ/kg/K"

@test SpecificCP(3.0, 300.0) ≈ 0.417_301_218E1
@test SpecificCP(3.0u"MPa", 300.0u"K") ≈ 0.417_301_218E1u"kJ/kg/K"

@test SpecificCP(80.0, 300.0) ≈ 0.401_008_987E1
@test SpecificCP(80.0u"MPa", 300.0u"K") ≈ 0.401_008_987E1u"kJ/kg/K"

@test SpecificCP(3.0, 500.0) ≈ 0.465_580_682E1
@test SpecificCP(3.0u"MPa", 500.0u"K") ≈ 0.465_580_682E1u"kJ/kg/K"

@test SpeedOfSound(3.0, 300.0) ≈ 0.150_773_921E4
@test SpeedOfSound(3.0u"MPa", 300.0u"K") ≈ 0.150_773_921E4u"m/s"

@test SpeedOfSound(80.0, 300.0) ≈ 0.163_469_054E4
@test SpeedOfSound(80.0u"MPa", 300.0u"K") ≈ 0.163_469_054E4u"m/s"

@test SpeedOfSound(3.0, 500.0) ≈ 0.124_071_337E4
@test SpeedOfSound(3.0u"MPa", 500.0u"K") ≈ 0.124_071_337E4u"m/s"
