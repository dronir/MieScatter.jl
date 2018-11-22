
# MieScatter

Compute Mie scattering in Julia. Mie scattering is the scattering of an electromagnetic
plane wave by a homogeneous sphere.

Based on a Fortran code by Karri Muinonen.

Now updated for Julia 1.0.

## Usage

The work is done by the `compute_mie` function. There are two ways to call it:

```
using MieScatter
S, Qsca, Qext, Qback = compute_mie(x, m, N)
S, Qsca, Qext, Qback = compute_mie(x, m, list_of_angles)
```

where `x` is the size parameter (two pi times the sphere radius divided by the wavelength)
and `m` is the (possibly complex-valued) relative refractive index of the sphere (the
refractive index of the sphere divided by the refractive index of the infinite medium
containing it).

In the first form, `N` is the number of different scattering angles to compute. The range
from 0 to 180 degrees in scattering angle will be split in `N` steps. N needs to be at
least 2.

In the second form, the third parameter can be a list of numbers, which should be between 0
and pi. The Mie computation is performed using these numbers as scattering angles. The list
can contain only one number. For example, to compute only backscattering, you can call
`compute_mie(x, m, [pi])`.

The returned matrix `S` is an `N`x4 matrix containing the scattering matrix elements S11,
S12, S33 and S34 (in this order), for all `N` scattering angles. `Qsca` is the scattering
coefficient, `Qext` the extinction coefficient and `Qback` backscattering coefficient.
Multiplying these coefficients by pi times the radius of the particle squared, you get the
cross sections. The absorbtion coefficient `Qabs = Qext-Qsca`.

Note that if the list of scattering angles is provided by the user, Qext and Qback can only
be computed if the list includes zero (forward scattering) and pi (backscattering). These
variables will return `NaN` if this is not true.

For convenience, there is also a function `size_parameter(r, lambda)`, which just returns
the size parameter `x = 2*pi*r / lambda`.
