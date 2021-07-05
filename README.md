# Pollard-Rho and General Factorization of Primes
This implements a complete prime factorization of any composite N using trial division, pollard p-1, and different versions of pollard rho.
It uses the [GMP](https://gmplib.org/) library for efficient large scale calculations, in particular [gmpxx](https://gmplib.org/manual/C_002b_002b-Interface-General) for integers and [mpreal](https://github.com/advanpix/mpreal) for real number calculations.
Additionally, premade tests are implemented, some of which use the [gnuplot-iostream](https://github.com/dstahlke/gnuplot-iostream) library for graphic output.

A short rundown of the consecutive algorithms used when calling ```Factorize::findFactors()```:
1. Remove small prime factors using consecutive trial division
2. Remove (and factor in 3.) cryptographically weak factors using Pollard (p-1)
3. Recursively find all remaining prime factors using a given variant of Pollard Rho (Floyd's improved algorithm by default)
   - Current Pollard Rho implementations: ```Floyd()```
