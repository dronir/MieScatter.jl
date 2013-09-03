#!/usr/bin/env julia

using Base.Test
using MieScatter

# These tests are made compared to the FORTRAN77 Mie code
# written by Karri Muinonen to the precision of the output
# of that code.

# Test computation of scattering, extinction and absorbsion
# cross sections.

# Size parameter 2, refractive index 1.33
S, Qsca, Qext, Qback = compute_mie(2.0, 1.33, 181)
Qabs = Qext-Qsca
@test_approx_eq_eps Qsca 0.71295 1e-5
@test_approx_eq_eps Qext 0.71295 1e-5
@test_approx_eq_eps Qabs 0.0 1e-5

# Size parameter 2, refractive index 1.6
S, Qsca, Qext, Qback = compute_mie(2.0, 1.6, 181)
Qabs = Qext-Qsca
@test_approx_eq_eps Qsca 2.4834 1e-4
@test_approx_eq_eps Qext 2.4834 1e-4
@test_approx_eq_eps Qabs 0.0 1e-4

# Size parameter 20, refractive index 1.33
S, Qsca, Qext, Qback = compute_mie(20.0, 1.33, 181)
Qabs = Qext-Qsca
@test_approx_eq_eps Qsca 0.21401E+01 1e-4
@test_approx_eq_eps Qext 0.21401E+01 1e-4
@test_approx_eq_eps Qabs 0.66613E-14 1e-4

# Size parameter 20, refractive index 1.6
S, Qsca, Qext, Qback = compute_mie(20.0, 1.6, 181)
Qabs = Qext-Qsca
@test_approx_eq_eps Qsca 0.26115E+01 1e-4
@test_approx_eq_eps Qext 0.26115E+01 1e-4
@test_approx_eq_eps Qabs 0.39968E-14 1e-4
