using SteamTables
using Unitful
using Test

mysignif(x, n) = round(x; sigdigits = n)
function mysignif(x::Q, n::Int) where Q <: Quantity
    y = Quantity(round(x.val, sigdigits = n), unit(x))
    return y
end


const TestDigits = 6 # number of significant digits to test selected functions to

@testset "Region 1, forwards equations for P and T              " begin
    include("testreg1.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 1, backwards equations for P and h             " begin
    include("testreg1Ph.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 1, backwards equations for P and s             " begin
    include("testreg1Ps.jl")
end
@testset "Region 2, forwards equations for P and T              " begin
    include("testreg2.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 2a, backwards equations for P and h            " begin
    include("testreg2aPh.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 2b, backwards equations for P and h            " begin
    include("testreg2bPh.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 2c, backwards equations for P and h            " begin
    include("testreg2cPh.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 2a, backwards equations for P and s            " begin
    include("testreg2aPs.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
#@test mysignif(SpecificG(8.0,  0.600_484_040E3), TestDigits)  ≈ mysignif(SpecificG_Ps(8.0,  6.0), TestDigits)
@testset "Region 2b, backwards equations for P and s            " begin
    include("testreg2bPs.jl")
end

# The free energies have lower consistency between the forward and backwards equations. Check to fewer significant digits.
@testset "Region 2c, backwards equations for P and s            " begin
    include("testreg2cPs.jl")
end

# Region 3 is specified in terms of density rather than pressure. The utility functions fix this by solving
# for the pressure. This does introduce some small error, which means we should check to fewer significant digits.
# For consistency the checks here are done to the same number of digits as for the free energies, although the
# accuracy here will be better.
@testset "Region 3, forwards equations for P and T              " begin
    include("testreg3.jl")
end

@testset "Region 4 - saturation pressure from temperature       " begin
    include("testreg4T.jl")
end

@testset "Region 4 - saturation temperature from pressure       " begin
    include("testreg4P.jl")
end

@testset "Region 5 - forwards equations for P and T             " begin
    include("testreg5.jl")
end
