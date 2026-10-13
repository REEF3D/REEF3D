/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef WAVE_LIB_FASTEXP_H_
#define WAVE_LIB_FASTEXP_H_

#include<cstdint>
#include<cstring>

/*--------------------------------------------------------------------
exp(x) for the vertical profiles of the irregular-wave components
(wave_lib_irregular_1st, wave_lib_irregular_2nd_cache.h).

Branch-free, so that a loop over the components vectorises with any
compiler and instruction set (SSE2 / AVX / NEON), unlike a call of the
libm exp. x = n ln2 + r, |r| <= ln2/2, exp(r) from its Taylor series to
r^13 (truncation below 1e-17), 2^n from the exponent bits. Relative error
within 2.3e-16 (1 ulp) of exp() for x in [-708, 709]; x must be in that
range (no check: a clamp keeps GCC from vectorising the loops without
AVX-512). The callers check the range per column and use exp() outside.
With AVX2 0.8 ns against 4 ns for exp(), with SSE2 about 3 ns.
--------------------------------------------------------------------*/

inline double wave_exp(double x)
{
    // n = round(x / ln2): the shift puts n into the low mantissa bits of kd
    const double shift = 6755399441055744.0;           // 1.5 * 2^52
    const double kd = x*1.4426950408889634074 + shift;
    const double n = kd - shift;

    // r = x - n ln2 with ln2 in two parts (n ln2_hi exact)
    const double r = (x - n*6.93147180369123816490e-01) - n*1.90821492927058770002e-10;

    // exp(r), Horner
    double s = 1.0/6227020800.0;                       // 1/13!
    s = s*r + 1.0/479001600.0;
    s = s*r + 1.0/39916800.0;
    s = s*r + 1.0/3628800.0;
    s = s*r + 1.0/362880.0;
    s = s*r + 1.0/40320.0;
    s = s*r + 1.0/5040.0;
    s = s*r + 1.0/720.0;
    s = s*r + 1.0/120.0;
    s = s*r + 1.0/24.0;
    s = s*r + 1.0/6.0;
    s = s*r + 0.5;
    s = s*r + 1.0;
    s = s*r + 1.0;

    // 2^n: n + 1023 in the exponent field (the bits of the shift drop out of the shift by 52)
    std::uint64_t ki;
    std::memcpy(&ki,&kd,sizeof(ki));
    const std::uint64_t eb = (ki + 1023u) << 52;
    double sc;
    std::memcpy(&sc,&eb,sizeof(sc));

    return s*sc;
}

#endif
