@test mysignif(SpecificV(0.255_837_018E2, 650.0), TestDigits)    ≈ 0.002
@test mysignif(SpecificV(0.222_930_643E2, 650.0), TestDigits)    ≈ 0.005
@test mysignif(SpecificV(0.783_095_639E2, 750.0), TestDigits)    ≈ 0.002
@test SpecificH(0.255_837_018E2, 650.0)    ≈ 0.186_343_019E4
@test SpecificH(0.222_930_643E2, 650.0)    ≈ 0.237_512_401E4
@test SpecificH(0.783_095_639E2, 750.0)    ≈ 0.225_868_845E4
@test SpecificU(0.255_837_018E2, 650.0)    ≈ 0.181_226_279E4
@test SpecificU(0.222_930_643E2, 650.0)    ≈ 0.226_365_868E4
@test SpecificU(0.783_095_639E2, 750.0)    ≈ 0.210_206_932E4
@test SpecificS(0.255_837_018E2, 650.0)    ≈ 0.405_427_273E1
@test SpecificS(0.222_930_643E2, 650.0)    ≈ 0.485_438_792E1
@test SpecificS(0.783_095_639E2, 750.0)    ≈ 0.446_971_906E1
@test SpecificCP(0.255_837_018E2, 650.0)   ≈ 0.138_935_717E2
@test SpecificCP(0.222_930_643E2, 650.0)   ≈ 0.446_579_373E2
@test SpecificCP(0.783_095_639E2, 750.0)   ≈ 0.634_165_359E1
@test SpeedOfSound(0.255_837_018E2, 650.0) ≈ 0.502_005_554E3
@test SpeedOfSound(0.222_930_643E2, 650.0) ≈ 0.383_444_594E3
@test SpeedOfSound(0.783_095_639E2, 750.0) ≈ 0.760_696_041E3