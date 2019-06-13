@test mysignif(SpecificG(3.0, 0.391_798_509E3), TestDigits) ≈ mysignif(SpecificG_Ph(3.0, 500.0), TestDigits)
@test mysignif(SpecificG(3.0u"MPa", 0.391_798_509E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(3.0u"MPa", 500.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(80.0, 0.378_108_626E3), TestDigits) ≈ mysignif(SpecificG_Ph(80.0, 500.0), TestDigits)
@test mysignif(SpecificG(80.0u"MPa", 0.378_108_626E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(80.0u"MPa", 500.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(80.0, 0.611_041_229E3), TestDigits) ≈ mysignif(SpecificG_Ph(80.0, 1500.0), TestDigits)
@test mysignif(SpecificG(80.0u"MPa", 0.611_041_229E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(80.0u"MPa", 1500.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(3.0, 0.391_798_509E3), TestDigits) ≈ mysignif(SpecificF_Ph(3.0, 500.0), TestDigits)
@test mysignif(SpecificF(3.0u"MPa", 0.391_798_509E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(3.0u"MPa", 500.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(80.0, 0.378_108_626E3), TestDigits) ≈ mysignif(SpecificF_Ph(80.0, 500.0), TestDigits)
@test mysignif(SpecificF(80.0u"MPa", 0.378_108_626E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(80.0u"MPa", 500.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(80.0, 0.611_041_229E3), TestDigits) ≈ mysignif(SpecificF_Ph(80.0, 1500.0), TestDigits)
@test mysignif(SpecificF(80.0u"MPa", 0.611_041_229E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(80.0u"MPa", 1500.0u"kJ/kg"), TestDigits)

@test SpecificV(3.0, 0.391_798_509E3) ≈ SpecificV_Ph(3.0, 500.0)
@test SpecificV(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificV_Ph(3.0u"MPa", 500.0u"kJ/kg")

@test SpecificV(80.0, 0.378_108_626E3) ≈ SpecificV_Ph(80.0, 500.0)
@test SpecificV(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificV_Ph(80.0u"MPa", 500.0u"kJ/kg")

@test SpecificV(80.0, 0.611_041_229E3) ≈ SpecificV_Ph(80.0, 1500.0)
@test SpecificV(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificV_Ph(80.0u"MPa", 1500.0u"kJ/kg")

@test SpecificU(3.0, 0.391_798_509E3) ≈ SpecificU_Ph(3.0, 500.0)
@test SpecificU(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificU_Ph(3.0u"MPa", 500.0u"kJ/kg")

@test SpecificU(80.0, 0.378_108_626E3) ≈ SpecificU_Ph(80.0, 500.0)
@test SpecificU(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificU_Ph(80.0u"MPa", 500.0u"kJ/kg")

@test SpecificU(80.0, 0.611_041_229E3) ≈ SpecificU_Ph(80.0, 1500.0)
@test SpecificU(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificU_Ph(80.0u"MPa", 1500.0u"kJ/kg")

@test SpecificS(3.0, 0.391_798_509E3) ≈ SpecificS_Ph(3.0, 500.0)
@test SpecificS(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificS_Ph(3.0u"MPa", 500.0u"kJ/kg")

@test SpecificS(80.0, 0.378_108_626E3) ≈ SpecificS_Ph(80.0, 500.0)
@test SpecificS(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificS_Ph(80.0u"MPa", 500.0u"kJ/kg")

@test SpecificS(80.0, 0.611_041_229E3) ≈ SpecificS_Ph(80.0, 1500.0)
@test SpecificS(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificS_Ph(80.0u"MPa", 1500.0u"kJ/kg")

@test SpecificCP(3.0, 0.391_798_509E3) ≈ SpecificCP_Ph(3.0, 500.0)
@test SpecificCP(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificCP_Ph(3.0u"MPa", 500.0u"kJ/kg")

@test SpecificCP(80.0, 0.378_108_626E3) ≈ SpecificCP_Ph(80.0, 500.0)
@test SpecificCP(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificCP_Ph(80.0u"MPa", 500.0u"kJ/kg")

@test SpecificCP(80.0, 0.611_041_229E3) ≈ SpecificCP_Ph(80.0, 1500.0)
@test SpecificCP(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificCP_Ph(80.0u"MPa", 1500.0u"kJ/kg")

@test SpecificCV(3.0, 0.391_798_509E3) ≈ SpecificCV_Ph(3.0, 500.0)
@test SpecificCV(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpecificCV_Ph(3.0u"MPa", 500.0u"kJ/kg")

@test SpecificCV(80.0, 0.378_108_626E3) ≈ SpecificCV_Ph(80.0, 500.0)
@test SpecificCV(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpecificCV_Ph(80.0u"MPa", 500.0u"kJ/kg")

@test SpecificCV(80.0, 0.611_041_229E3) ≈ SpecificCV_Ph(80.0, 1500.0)
@test SpecificCV(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpecificCV_Ph(80.0u"MPa", 1500.0u"kJ/kg")

@test SpeedOfSound(3.0, 0.391_798_509E3) ≈ SpeedOfSound_Ph(3.0, 500.0)
@test SpeedOfSound(3.0u"MPa", 0.391_798_509E3u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 500.0u"kJ/kg")

@test SpeedOfSound(80.0, 0.378_108_626E3) ≈ SpeedOfSound_Ph(80.0, 500.0)
@test SpeedOfSound(80.0u"MPa", 0.378_108_626E3u"K") ≈ SpeedOfSound_Ph(80.0u"MPa", 500.0u"kJ/kg")

@test SpeedOfSound(80.0, 0.611_041_229E3) ≈ SpeedOfSound_Ph(80.0, 1500.0)
@test SpeedOfSound(80.0u"MPa", 0.611_041_229E3u"K") ≈ SpeedOfSound_Ph(80.0u"MPa", 1500.0u"kJ/kg")
