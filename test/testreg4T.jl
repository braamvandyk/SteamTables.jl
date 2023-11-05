@test Psat(300.0) ≈ 0.353_658_941E-2
@test Psat(300.0u"K") ≈ 0.353_658_941E-2u"MPa"

@test Psat(500.0) ≈ 0.263_889_776E1
@test Psat(500.0u"K") ≈ 0.263_889_776E1u"MPa"

@test Psat(600.0) ≈ 0.123_443_146E2
@test Psat(600.0u"K") ≈ 0.123_443_146E2u"MPa"

@test mysignif(SatDensL(273.16), 6) ≈ 999.789
@test mysignif(SatDensL(273.16u"K"), 6) ≈ 999.789u"kg/m^3"

@test mysignif(SatDensL(373.1243), 6) ≈ 958.365
@test mysignif(SatDensL(373.1243u"K"), 6) ≈ 958.365u"kg/m^3"

@test mysignif(SatDensL(647.096), 3) ≈ 322
@test mysignif(SatDensL(647.096u"K"), 3) ≈ 322u"kg/m^3"

@test mysignif(SatDensV(273.16), 6) ≈ 0.00485426
@test mysignif(SatDensV(273.16u"K"), 6) ≈ 0.00485426u"kg/m^3"

@test mysignif(SatDensV(373.1243), 6) ≈ 0.597586
@test mysignif(SatDensV(373.1243u"K"), 6) ≈ 0.597586u"kg/m^3"

@test mysignif(SatDensV(647.096), 3) ≈ 322
@test mysignif(SatDensV(647.096u"K"), 3) ≈ 322u"kg/m^3"

@test mysignif(SatHL(273.16), 6) ≈ mysignif(0.000611786, 6)
@test mysignif(SatHL(273.16u"K"), 6) ≈ mysignif(0.000611786u"kJ/kg", 6)

@test mysignif(SatHL(373.1243), 5) ≈ mysignif(419.05, 5)
@test mysignif(SatHL(373.1243u"K"), 5) ≈ mysignif(419.05u"kJ/kg", 5)

@test mysignif(SatHL(647.096), 5) ≈ mysignif(2086.6, 5)
@test mysignif(SatHL(647.096u"K"), 5) ≈ mysignif(2086.6u"kJ/kg", 5)


@test mysignif(SatHV(273.16), 5) ≈ 2500.5
@test mysignif(SatHV(273.16u"K"), 5) ≈ 2500.5u"kJ/kg"

@test mysignif(SatHV(373.1243), 5) ≈ 2675.7
@test mysignif(SatHV(373.1243u"K"), 5) ≈ 2675.7u"kJ/kg"

@test mysignif(SatHV(647.096), 5) ≈ 2086.6
@test mysignif(SatHV(647.096u"K"), 5) ≈ 2086.6u"kJ/kg"

@test abs(SatSL(273.16)) < 1e-8
@test abs(SatSL(273.16u"K")) < 1e-8u"kJ/kg/K"

@test mysignif(SatSL(373.1243), 4) ≈ 1.307
@test mysignif(SatSL(373.1243u"K"), 4) ≈ 1.307u"kJ/kg/K"

@test mysignif(SatSL(647.096), 4) ≈ 4.410
@test mysignif(SatSL(647.096u"K"), 4) ≈ 4.410u"kJ/kg/K"

@test mysignif(SatSV(273.16), 4) ≈ 9.154
@test mysignif(SatSV(273.16u"K"), 4) ≈ 9.154u"kJ/kg/K"

@test mysignif(SatSV(373.1243), 4) ≈ 7.355
@test mysignif(SatSV(373.1243u"K"), 4) ≈ 7.355u"kJ/kg/K"

@test mysignif(SatSV(647.096), 4) ≈ 4.410
@test mysignif(SatSV(647.096u"K"), 4) ≈ 4.410u"kJ/kg/K"