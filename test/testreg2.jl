@test SpecificV(0.0035, 300.0) ≈ 0.394_913_866E2
@test SpecificV(0.0035u"MPa", 300.0u"K") ≈ 0.394_913_866E2u"m^3/kg"

@test SpecificV(0.0035, 700.0) ≈ 0.923_015_898E2
@test SpecificV(0.0035u"MPa", 700.0u"K") ≈ 0.923_015_898E2u"m^3/kg"

@test SpecificV(30.0, 700.0) ≈ 0.542_946_619E-2
@test SpecificV(30.0u"MPa", 700.0u"K") ≈ 0.542_946_619E-2u"m^3/kg"

@test SpecificH(0.0035, 300.0) ≈ 0.254_991_145E4
@test SpecificH(0.0035u"MPa", 300.0u"K") ≈ 0.254_991_145E4u"kJ/kg"

@test SpecificH(0.0035, 700.0) ≈ 0.333_568_375E4
@test SpecificH(0.0035u"MPa", 700.0u"K") ≈ 0.333_568_375E4u"kJ/kg"

@test SpecificH(30.0, 700.0) ≈ 0.263_149_474E4
@test SpecificH(30.0u"MPa", 700.0u"K") ≈ 0.263_149_474E4u"kJ/kg"

@test SpecificU(0.0035, 300.0) ≈ 0.241_169_160E4
@test SpecificU(0.0035u"MPa", 300.0u"K") ≈ 0.241_169_160E4u"kJ/kg"

@test SpecificU(0.0035, 700.0) ≈ 0.301_262_819E4
@test SpecificU(0.0035u"MPa", 700.0u"K") ≈ 0.301_262_819E4u"kJ/kg"

@test SpecificU(30.0, 700.0) ≈ 0.246_861_076E4
@test SpecificU(30.0u"MPa", 700.0u"K") ≈ 0.246_861_076E4u"kJ/kg"

@test SpecificS(0.0035, 300.0) ≈ 0.852_238_967E1
@test SpecificS(0.0035u"MPa", 300.0u"K") ≈ 0.852_238_967E1u"kJ/kg/K"

@test SpecificS(0.0035, 700.0) ≈ 0.101_749_996E2
@test SpecificS(0.0035u"MPa", 700.0u"K") ≈ 0.101_749_996E2u"kJ/kg/K"

@test SpecificS(30.0, 700.0) ≈ 0.517_540_298E1
@test SpecificS(30.0u"MPa", 700.0u"K") ≈ 0.517_540_298E1u"kJ/kg/K"

@test SpecificCP(0.0035, 300.0) ≈ 0.191_300_162E1
@test SpecificCP(0.0035u"MPa", 300.0u"K") ≈ 0.191_300_162E1u"kJ/kg/K"

@test SpecificCP(0.0035, 700.0) ≈ 0.208_141_274E1
@test SpecificCP(0.0035u"MPa", 700.0u"K") ≈ 0.208_141_274E1u"kJ/kg/K"

@test SpecificCP(30.0, 700.0) ≈ 0.103_505_092E2
@test SpecificCP(30.0u"MPa", 700.0u"K") ≈ 0.103_505_092E2u"kJ/kg/K"

@test SpeedOfSound(0.0035, 300.0) ≈ 0.427_920_172E3
@test SpeedOfSound(0.0035u"MPa", 300.0u"K") ≈ 0.427_920_172E3u"m/s"

@test SpeedOfSound(0.0035, 700.0) ≈ 0.644_289_068E3
@test SpeedOfSound(0.0035u"MPa", 700.0u"K") ≈ 0.644_289_068E3u"m/s"

@test SpeedOfSound(30.0, 700.0) ≈ 0.480_386_523E3
@test SpeedOfSound(30.0u"MPa", 700.0u"K") ≈ 0.480_386_523E3u"m/s"
