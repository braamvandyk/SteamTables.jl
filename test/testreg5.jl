@test SpecificV(0.5, 1500.0) ≈ 0.138_455_090E1
@test SpecificV(0.5u"MPa", 1500.0u"K") ≈ 0.138_455_090E1u"m^3/kg"

@test SpecificV(30.0, 1500.0) ≈ 0.230_761_299E-1
@test SpecificV(30.0u"MPa", 1500.0u"K") ≈ 0.230_761_299E-1u"m^3/kg"

@test SpecificV(30.0, 2000.0) ≈ 0.311_385_219E-1
@test SpecificV(30.0u"MPa", 2000.0u"K") ≈ 0.311_385_219E-1u"m^3/kg"

@test SpecificH(0.5, 1500.0) ≈ 0.521_976_855E4
@test SpecificH(0.5u"MPa", 1500.0u"K") ≈ 0.521_976_855E4u"kJ/kg"

@test SpecificH(30.0, 1500.0) ≈ 0.516_723_514E4
@test SpecificH(30.0u"MPa", 1500.0u"K") ≈ 0.516_723_514E4u"kJ/kg"

@test SpecificH(30.0, 2000.0) ≈ 0.657_122_604E4
@test SpecificH(30.0u"MPa", 2000.0u"K") ≈ 0.657_122_604E4u"kJ/kg"

@test SpecificU(0.5, 1500.0) ≈ 0.452_749_310E4
@test SpecificU(0.5u"MPa", 1500.0u"K") ≈ 0.452_749_310E4u"kJ/kg"

@test SpecificU(30.0, 1500.0) ≈ 0.447_495_124E4
@test SpecificU(30.0u"MPa", 1500.0u"K") ≈ 0.447_495_124E4u"kJ/kg"

@test SpecificU(30.0, 2000.0) ≈ 0.563_707_038E4
@test SpecificU(30.0u"MPa", 2000.0u"K") ≈ 0.563_707_038E4u"kJ/kg"

@test SpecificS(0.5, 1500.0) ≈ 0.965_408_875E1
@test SpecificS(0.5u"MPa", 1500.0u"K") ≈ 0.965_408_875E1u"kJ/kg/K"

@test SpecificS(30.0, 1500.0) ≈ 0.772_970_133E1
@test SpecificS(30.0u"MPa", 1500.0u"K") ≈ 0.772_970_133E1u"kJ/kg/K"

@test SpecificS(30.0, 2000.0) ≈ 0.853_640_523E1
@test SpecificS(30.0u"MPa", 2000.0u"K") ≈ 0.853_640_523E1u"kJ/kg/K"

@test SpecificCP(0.5, 1500.0) ≈ 0.261_609_445E1
@test SpecificCP(0.5u"MPa", 1500.0u"K") ≈ 0.261_609_445E1u"kJ/kg/K"

@test SpecificCP(30.0, 1500.0) ≈ 0.272_724_317E1
@test SpecificCP(30.0u"MPa", 1500.0u"K") ≈ 0.272_724_317E1u"kJ/kg/K"

@test SpecificCP(30.0, 2000.0) ≈ 0.288_569_882E1
@test SpecificCP(30.0u"MPa", 2000.0u"K") ≈ 0.288_569_882E1u"kJ/kg/K"

@test SpeedOfSound(0.5, 1500.0) ≈ 0.917_068_690E3
@test SpeedOfSound(0.5u"MPa", 1500.0u"K") ≈ 0.917_068_690E3u"m/s"

@test SpeedOfSound(30.0, 1500.0) ≈ 0.928_548_002E3
@test SpeedOfSound(30.0u"MPa", 1500.0u"K") ≈ 0.928_548_002E3u"m/s"

@test SpeedOfSound(30.0, 2000.0) ≈ 0.106_736_948E4
@test SpeedOfSound(30.0u"MPa", 2000.0u"K") ≈ 0.106_736_948E4u"m/s"

#issue #23
let P = 30,T = 2000
    H = SpecificH(P,T)
    S = SpecificS(P,T)
    @test S ≈ SpecificS_Ph(P,H)
    @test H ≈ SpecificH_Ps(P,S)
end