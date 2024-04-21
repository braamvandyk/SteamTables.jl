P = 30
T = 650
H = SpecificH(P,T)
S = SpecificS(P,T)
@test S ≈ SpecificS_Ph(P,H)
@test H ≈ SpecificH_Ps(P,S)
