@test mysignif(SpecificG(20.0, 0.697_992_849E3), TestDigits) ≈ mysignif(SpecificG_Ps(20.0, 5.75), TestDigits)
@test mysignif(SpecificG(20.0u"MPa", 0.697_992_849E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(20.0u"MPa", 5.75u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificG(80.0, 0.854_011_484E3), TestDigits) ≈ mysignif(SpecificG_Ps(80.0, 5.25), TestDigits)
@test mysignif(SpecificG(80.0u"MPa", 0.854_011_484E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(80.0u"MPa", 5.25u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificG(80.0, 0.949_017_998E3), TestDigits) ≈ mysignif(SpecificG_Ps(80.0, 5.75), TestDigits)
@test mysignif(SpecificG(80.0u"MPa", 0.949_017_998E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(80.0u"MPa", 5.75u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(20.0, 0.697_992_849E3), TestDigits) ≈ mysignif(SpecificF_Ps(20.0, 5.75), TestDigits)
@test mysignif(SpecificF(20.0u"MPa", 0.697_992_849E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(20.0u"MPa", 5.75u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(80.0, 0.854_011_484E3), TestDigits) ≈ mysignif(SpecificF_Ps(80.0, 5.25), TestDigits)
@test mysignif(SpecificF(80.0u"MPa", 0.854_011_484E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(80.0u"MPa", 5.25u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(80.0, 0.949_017_998E3), TestDigits) ≈ mysignif(SpecificF_Ps(80.0, 5.75), TestDigits)
@test mysignif(SpecificF(80.0u"MPa", 0.949_017_998E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(80.0u"MPa", 5.75u"kJ/kg/K"), TestDigits)

@test SpecificV(20.0, 0.697_992_849E3) ≈ SpecificV_Ps(20.0, 5.75)
@test SpecificV(20.0u"MPa", 0.697_992_849E3u"K") ≈ SpecificV_Ps(20.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificV(80.0, 0.854_011_484E3) ≈ SpecificV_Ps(80.0, 5.25)
@test SpecificV(80.0u"MPa", 0.854_011_484E3u"K") ≈ SpecificV_Ps(80.0u"MPa", 5.25u"kJ/kg/K")

@test SpecificV(80.0, 0.949_017_998E3) ≈ SpecificV_Ps(80.0, 5.75)
@test SpecificV(80.0u"MPa", 0.949_017_998E3u"K") ≈ SpecificV_Ps(80.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificU(20.0, 0.697_992_849E3) ≈ SpecificU_Ps(20.0, 5.75)
@test SpecificU(20.0u"MPa", 0.697_992_849E3u"K") ≈ SpecificU_Ps(20.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificU(80.0, 0.854_011_484E3) ≈ SpecificU_Ps(80.0, 5.25)
@test SpecificU(80.0u"MPa", 0.854_011_484E3u"K") ≈ SpecificU_Ps(80.0u"MPa", 5.25u"kJ/kg/K")

@test SpecificU(80.0, 0.949_017_998E3) ≈ SpecificU_Ps(80.0, 5.75)
@test SpecificU(80.0u"MPa", 0.949_017_998E3u"K") ≈ SpecificU_Ps(80.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificH(20.0, 0.697_992_849E3) ≈ SpecificH_Ps(20.0, 5.75)
@test SpecificH(20.0u"MPa", 0.697_992_849E3u"K") ≈ SpecificH_Ps(20.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificH(80.0, 0.854_011_484E3) ≈ SpecificH_Ps(80.0, 5.25)
@test SpecificH(80.0u"MPa", 0.854_011_484E3u"K") ≈ SpecificH_Ps(80.0u"MPa", 5.25u"kJ/kg/K")

@test SpecificH(80.0, 0.949_017_998E3) ≈ SpecificH_Ps(80.0, 5.75)
@test SpecificH(80.0u"MPa", 0.949_017_998E3u"K") ≈ SpecificH_Ps(80.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificCP(20.0, 0.697_992_849E3) ≈ SpecificCP_Ps(20.0, 5.75)
@test SpecificCP(20.0u"MPa", 0.697_992_849E3u"K") ≈ SpecificCP_Ps(20.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificCP(80.0, 0.854_011_484E3) ≈ SpecificCP_Ps(80.0, 5.25)
@test SpecificCP(80.0u"MPa", 0.854_011_484E3u"K") ≈ SpecificCP_Ps(80.0u"MPa", 5.25u"kJ/kg/K")

@test SpecificCP(80.0, 0.949_017_998E3) ≈ SpecificCP_Ps(80.0, 5.75)
@test SpecificCP(80.0u"MPa", 0.949_017_998E3u"K") ≈ SpecificCP_Ps(80.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificCV(20.0, 0.697_992_849E3) ≈ SpecificCV_Ps(20.0, 5.75)
@test SpecificCV(20.0u"MPa", 0.697_992_849E3u"K") ≈ SpecificCV_Ps(20.0u"MPa", 5.75u"kJ/kg/K")

@test SpecificCV(80.0, 0.854_011_484E3) ≈ SpecificCV_Ps(80.0, 5.25)
@test SpecificCV(80.0u"MPa", 0.854_011_484E3u"K") ≈ SpecificCV_Ps(80.0u"MPa", 5.25u"kJ/kg/K")

@test SpecificCV(80.0, 0.949_017_998E3) ≈ SpecificCV_Ps(80.0, 5.75)
@test SpecificCV(80.0u"MPa", 0.949_017_998E3u"K") ≈ SpecificCV_Ps(80.0u"MPa", 5.75u"kJ/kg/K")

@test SpeedOfSound(20.0, 0.697_992_849E3) ≈ SpeedOfSound_Ps(20.0, 5.75)
@test SpeedOfSound(20.0u"MPa", 0.697_992_849E3u"K") ≈ SpeedOfSound_Ps(20.0u"MPa", 5.75u"kJ/kg/K")

@test SpeedOfSound(80.0, 0.854_011_484E3) ≈ SpeedOfSound_Ps(80.0, 5.25)
@test SpeedOfSound(80.0u"MPa", 0.854_011_484E3u"K") ≈ SpeedOfSound_Ps(80.0u"MPa", 5.25u"kJ/kg/K")

@test SpeedOfSound(80.0, 0.949_017_998E3) ≈ SpeedOfSound_Ps(80.0, 5.75)
@test SpeedOfSound(80.0u"MPa", 0.949_017_998E3u"K") ≈ SpeedOfSound_Ps(80.0u"MPa", 5.75u"kJ/kg/K")
