@test mysignif(SpecificG(5.0, 0.801_299_102E3), TestDigits) ≈ mysignif(SpecificG_Ph(5.0, 3500), TestDigits)
@test mysignif(SpecificG(5.0u"MPa", 0.801_299_102E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(5.0u"MPa", 3500u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(5.0, 0.101_531_583E4), TestDigits) ≈ mysignif(SpecificG_Ph(5.0, 4000), TestDigits)
@test mysignif(SpecificG(5.0u"MPa", 0.101_531_583E4u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(5.0u"MPa", 4000u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(25.0, 0.875_279_054E3), TestDigits) ≈ mysignif(SpecificG_Ph(25.0, 3500), TestDigits)
@test mysignif(SpecificG(25.0u"MPa", 0.875_279_054E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(25.0u"MPa", 3500u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(5.0, 0.801_299_102E3), TestDigits) ≈ mysignif(SpecificF_Ph(5.0, 3500), TestDigits)
@test mysignif(SpecificF(5.0u"MPa", 0.801_299_102E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(5.0u"MPa", 3500u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(5.0, 0.101_531_583E4), TestDigits) ≈ mysignif(SpecificF_Ph(5.0, 4000), TestDigits)
@test mysignif(SpecificF(5.0u"MPa", 0.101_531_583E4u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(5.0u"MPa", 4000u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(25.0, 0.875_279_054E3), TestDigits) ≈ mysignif(SpecificF_Ph(25.0, 3500), TestDigits)
@test mysignif(SpecificF(25.0u"MPa", 0.875_279_054E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(25.0u"MPa", 3500u"kJ/kg"), TestDigits)

@test SpecificV(5.0, 0.801_299_102E3) ≈ SpecificV_Ph(5.0, 3500)
@test SpecificV(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificV_Ph(5.0u"MPa", 3500u"kJ/kg")

@test SpecificV(5.0, 0.101_531_583E4) ≈ SpecificV_Ph(5.0, 4000)
@test SpecificV(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificV_Ph(5.0u"MPa", 4000u"kJ/kg")

@test SpecificV(25.0, 0.875_279_054E3) ≈ SpecificV_Ph(25.0, 3500)
@test SpecificV(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificV_Ph(25.0u"MPa", 3500u"kJ/kg")

@test SpecificU(5.0, 0.801_299_102E3) ≈ SpecificU_Ph(5.0, 3500)
@test SpecificU(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificU_Ph(5.0u"MPa", 3500u"kJ/kg")

@test SpecificU(5.0, 0.101_531_583E4) ≈ SpecificU_Ph(5.0, 4000)
@test SpecificU(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificU_Ph(5.0u"MPa", 4000u"kJ/kg")

@test SpecificU(25.0, 0.875_279_054E3) ≈ SpecificU_Ph(25.0, 3500)
@test SpecificU(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificU_Ph(25.0u"MPa", 3500u"kJ/kg")

@test SpecificS(5.0, 0.801_299_102E3) ≈ SpecificS_Ph(5.0, 3500)
@test SpecificS(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificS_Ph(5.0u"MPa", 3500u"kJ/kg")

@test SpecificS(5.0, 0.101_531_583E4) ≈ SpecificS_Ph(5.0, 4000)
@test SpecificS(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificS_Ph(5.0u"MPa", 4000u"kJ/kg")

@test SpecificS(25.0, 0.875_279_054E3) ≈ SpecificS_Ph(25.0, 3500)
@test SpecificS(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificS_Ph(25.0u"MPa", 3500u"kJ/kg")

@test SpecificCP(5.0, 0.801_299_102E3) ≈ SpecificCP_Ph(5.0, 3500)
@test SpecificCP(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificCP_Ph(5.0u"MPa", 3500u"kJ/kg")

@test SpecificCP(5.0, 0.101_531_583E4) ≈ SpecificCP_Ph(5.0, 4000)
@test SpecificCP(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificCP_Ph(5.0u"MPa", 4000u"kJ/kg")

@test SpecificCP(25.0, 0.875_279_054E3) ≈ SpecificCP_Ph(25.0, 3500)
@test SpecificCP(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificCP_Ph(25.0u"MPa", 3500u"kJ/kg")

@test SpecificCV(5.0, 0.801_299_102E3) ≈ SpecificCV_Ph(5.0, 3500)
@test SpecificCV(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpecificCV_Ph(5.0u"MPa", 3500u"kJ/kg")

@test SpecificCV(5.0, 0.101_531_583E4) ≈ SpecificCV_Ph(5.0, 4000)
@test SpecificCV(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpecificCV_Ph(5.0u"MPa", 4000u"kJ/kg")

@test SpecificCV(25.0, 0.875_279_054E3) ≈ SpecificCV_Ph(25.0, 3500)
@test SpecificCV(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpecificCV_Ph(25.0u"MPa", 3500u"kJ/kg")

@test SpeedOfSound(5.0, 0.801_299_102E3) ≈ SpeedOfSound_Ph(5.0, 3500)
@test SpeedOfSound(5.0u"MPa", 0.801_299_102E3u"K") ≈ SpeedOfSound_Ph(5.0u"MPa", 3500u"kJ/kg")

@test SpeedOfSound(5.0, 0.101_531_583E4) ≈ SpeedOfSound_Ph(5.0, 4000)
@test SpeedOfSound(5.0u"MPa", 0.101_531_583E4u"K") ≈ SpeedOfSound_Ph(5.0u"MPa", 4000u"kJ/kg")

@test SpeedOfSound(25.0, 0.875_279_054E3) ≈ SpeedOfSound_Ph(25.0, 3500)
@test SpeedOfSound(25.0u"MPa", 0.875_279_054E3u"K") ≈ SpeedOfSound_Ph(25.0u"MPa", 3500u"kJ/kg")
