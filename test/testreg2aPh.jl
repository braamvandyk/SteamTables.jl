@test mysignif(SpecificG(0.001, 0.534_433_241E3), TestDigits) ≈ mysignif(SpecificG_Ph(0.001, 3000), TestDigits)
@test mysignif(SpecificG(0.001u"MPa", 0.534_433_241E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(0.001u"MPa", 3000u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(3.0, 0.575_373_370E3), TestDigits) ≈ mysignif(SpecificG_Ph(3.0, 3000), TestDigits)
@test mysignif(SpecificG(3.0u"MPa", 0.575_373_370E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(3.0u"MPa", 3000u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(3.0, 0.101_077_577E4), TestDigits) ≈ mysignif(SpecificG_Ph(3.0, 4000), TestDigits)
@test mysignif(SpecificG(3.0u"MPa", 0.101_077_577E4u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(3.0u"MPa", 4000u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(0.001, 0.534_433_241E3), TestDigits) ≈ mysignif(SpecificF_Ph(0.001, 3000), TestDigits)
@test mysignif(SpecificF(0.001u"MPa", 0.534_433_241E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(0.001u"MPa", 3000u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(3.0, 0.575_373_370E3), TestDigits) ≈ mysignif(SpecificF_Ph(3.0, 3000), TestDigits)
@test mysignif(SpecificF(3.0u"MPa", 0.575_373_370E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(3.0u"MPa", 3000u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(3.0, 0.101_077_577E4), TestDigits) ≈ mysignif(SpecificF_Ph(3.0, 4000), TestDigits)
@test mysignif(SpecificF(3.0u"MPa", 0.101_077_577E4u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(3.0u"MPa", 4000u"kJ/kg"), TestDigits)

@test SpecificV(0.001, 0.534_433_241E3) ≈ SpecificV_Ph(0.001, 3000)
@test SpecificV(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificV_Ph(0.001u"MPa", 3000u"kJ/kg")

@test SpecificV(3.0, 0.575_373_370E3) ≈ SpecificV_Ph(3.0, 3000)
@test SpecificV(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificV_Ph(3.0u"MPa", 3000u"kJ/kg")

@test SpecificV(3.0, 0.101_077_577E4) ≈ SpecificV_Ph(3.0, 4000)
@test SpecificV(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificV_Ph(3.0u"MPa", 4000u"kJ/kg")

@test SpecificU(0.001, 0.534_433_241E3) ≈ SpecificU_Ph(0.001, 3000)
@test SpecificU(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificU_Ph(0.001u"MPa", 3000u"kJ/kg")

@test SpecificU(3.0, 0.575_373_370E3) ≈ SpecificU_Ph(3.0, 3000)
@test SpecificU(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificU_Ph(3.0u"MPa", 3000u"kJ/kg")

@test SpecificU(3.0, 0.101_077_577E4) ≈ SpecificU_Ph(3.0, 4000)
@test SpecificU(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificU_Ph(3.0u"MPa", 4000u"kJ/kg")

@test SpecificS(0.001, 0.534_433_241E3) ≈ SpecificS_Ph(0.001, 3000)
@test SpecificS(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificS_Ph(0.001u"MPa", 3000u"kJ/kg")

@test SpecificS(3.0, 0.575_373_370E3) ≈ SpecificS_Ph(3.0, 3000)
@test SpecificS(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificS_Ph(3.0u"MPa", 3000u"kJ/kg")

@test SpecificS(3.0, 0.101_077_577E4) ≈ SpecificS_Ph(3.0, 4000)
@test SpecificS(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificS_Ph(3.0u"MPa", 4000u"kJ/kg")

@test SpecificCP(0.001, 0.534_433_241E3) ≈ SpecificCP_Ph(0.001, 3000)
@test SpecificCP(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificCP_Ph(0.001u"MPa", 3000u"kJ/kg")

@test SpecificCP(3.0, 0.575_373_370E3) ≈ SpecificCP_Ph(3.0, 3000)
@test SpecificCP(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificCP_Ph(3.0u"MPa", 3000u"kJ/kg")

@test SpecificCP(3.0, 0.101_077_577E4) ≈ SpecificCP_Ph(3.0, 4000)
@test SpecificCP(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificCP_Ph(3.0u"MPa", 4000u"kJ/kg")

@test SpecificCV(0.001, 0.534_433_241E3) ≈ SpecificCV_Ph(0.001, 3000)
@test SpecificCV(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpecificCV_Ph(0.001u"MPa", 3000u"kJ/kg")

@test SpecificCV(3.0, 0.575_373_370E3) ≈ SpecificCV_Ph(3.0, 3000)
@test SpecificCV(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpecificCV_Ph(3.0u"MPa", 3000u"kJ/kg")

@test SpecificCV(3.0, 0.101_077_577E4) ≈ SpecificCV_Ph(3.0, 4000)
@test SpecificCV(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpecificCV_Ph(3.0u"MPa", 4000u"kJ/kg")

@test SpeedOfSound(0.001, 0.534_433_241E3) ≈ SpeedOfSound_Ph(0.001, 3000)
@test SpeedOfSound(0.001u"MPa", 0.534_433_241E3u"K") ≈ SpeedOfSound_Ph(0.001u"MPa", 3000u"kJ/kg")

@test SpeedOfSound(3.0, 0.575_373_370E3) ≈ SpeedOfSound_Ph(3.0, 3000)
@test SpeedOfSound(3.0u"MPa", 0.575_373_370E3u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 3000u"kJ/kg")

@test SpeedOfSound(3.0, 0.101_077_577E4) ≈ SpeedOfSound_Ph(3.0, 4000)
@test SpeedOfSound(3.0u"MPa", 0.101_077_577E4u"K") ≈ SpeedOfSound_Ph(3.0u"MPa", 4000u"kJ/kg")
