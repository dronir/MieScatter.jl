using MieScatter
using Test


@testset "Cross sections" begin

# These tests are made compared to the FORTRAN77 Mie code
# written by Karri Muinonen to the precision of the output
# of that code.

# Test computation of scattering, extinction and absorbsion
# cross sections.

# Size parameter 2, refractive index 1.33
S, Qsca, Qext, Qback = compute_mie(2.0, 1.33, 181)
Qabs = Qext-Qsca
@test isapprox(Qsca, 0.71295 ; atol = 6)
@test isapprox(Qext, 0.71295 ; atol = 6)
@test isapprox(Qabs, 0.0 ; atol = 10)

# Size parameter 2, refractive index 1.6
S, Qsca, Qext, Qback = compute_mie(2.0, 1.6, 181)
Qabs = Qext-Qsca
@test isapprox(Qsca, 2.4834 ; atol = 5)
@test isapprox(Qext, 2.4834 ; atol = 5)
@test isapprox(Qabs, 0.0 ; atol = 10)

# Size parameter 20, refractive index 1.33
S, Qsca, Qext, Qback = compute_mie(20.0, 1.33, 181)
Qabs = Qext-Qsca
@test isapprox(Qsca, 0.21401E+01 ; atol = 6)
@test isapprox(Qext, 0.21401E+01 ; atol = 6)
@test isapprox(Qabs, 0.66613E-14 ; atol = 6)

# Size parameter 20, refractive index 1.6
S, Qsca, Qext, Qback = compute_mie(20.0, 1.6, 181)
Qabs = Qext-Qsca
@test isapprox(Qsca, 0.26115E+01 ; atol = 6)
@test isapprox(Qext, 0.26115E+01 ; atol = 6)
@test isapprox(Qabs, 0.39968E-14 ; atol = 6)

end # cross section testset
