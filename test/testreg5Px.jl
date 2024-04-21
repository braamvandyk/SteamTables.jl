let P = 30,T = 2000
    H = SpecificH(P,T)
    S = SpecificS(P,T)
    @test S ≈ SpecificS_Ph(P,H)
    @test H ≈ SpecificH_Ps(P,S)
end
