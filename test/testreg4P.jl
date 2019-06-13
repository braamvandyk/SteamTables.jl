@test Tsat(0.10) ≈ 0.372_755_919E3
@test Tsat(0.10u"MPa") ≈ 0.372_755_919E3u"K"

@test Tsat(1.00) ≈ 0.453_035_632E3
@test Tsat(1.00u"MPa") ≈ 0.453_035_632E3u"K"

@test Tsat(10.0) ≈ 0.584_149_488E3
@test Tsat(10.0u"MPa") ≈ 0.584_149_488E3u"K"
