using SteamTables
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

const TestDigits = 6 # number of significant digits to test selected functions to

println("Region 1, forwards equations for P and T")
include("testreg1.jl")

println("Region 1, backwards equations for P and h")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg1Ph.jl")

println("Region 1, backwards equations for P and s")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg1Ps.jl")

println("Region 2, forwards equations for P and T")
include("testreg2.jl")

println("Region 2a, backwards equations for P and h")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg2aPh.jl")

println("Region 2b, backwards equations for P and h")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg2bPh.jl")

println("Region 2c, backwards equations for P and h")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg2cPh.jl")

println("Region 2a, backwards equations for P and s")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg2aPs.jl")

println("Region 2b, backwards equations for P and s")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
#@test signif(SpecificG(8.0,  0.600_484_040E3), TestDigits)  ≈ signif(SpecificG_Ps(8.0,  6.0), TestDigits)
include("testreg2bPs.jl")

println("Region 2c, backwards equations for P and s")
# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
include("testreg2cPs.jl")

println("Region 3, forwards equations for P and T")
# Region 3 is specified in terms of density rather than pressure. The utility functions fix this by solving
# for the pressure. This does introduce some small error, which means we should check to fewer significant digits.
# For consistency the checks here are done to the same number of digits as for the free energies, although the
# accuracy here will be better.
include("testreg3.jl")

println("Region 4 - saturation pressure from temperature")
include("testreg4T.jl")

println("Region 4 - saturation temperature from pressure")
include("testreg4P.jl")

println("Region 5 - forwards equations for P and T")
include("testreg5.jl")
