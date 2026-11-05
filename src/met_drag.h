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

#ifndef MET_DRAG_H_
#define MET_DRAG_H_

#include<algorithm>
#include<cmath>

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200): drag coefficient of the 10 m wind,
tau = rho_air Cd(|U10|) |U10| U10.  No lexer, unit test met_test.cpp.

W 207  1  Garratt (1977), capped at 3.5e-3 above 41 m/s (as A 570 1)
       2  Wu (1982)                                        (as A 570 2)
       3  Smith and Banke (1975)                           (as A 570 3)
       4  Zijlema et al. (2012), Uref = 31.5 m/s, >= 0     (as A 570 4)
       5  Garratt capped at W 213 (default 2.5e-3): the saturation of
          the drag in hurricane winds (Powell et al. 2003, Donelan et
          al. 2004), reached at 26.1 m/s; the default of the vortex W 200 3
       6  Zijlema et al. (2012), beyond its peak (2.03e-3 at 31.4 m/s)
          not below W 214 (default 1.5e-3): the decline of Powell et
          al. (2003) without the unphysical zero at 68 m/s
       7  wave-dependent, with the sea state of the SEASTATE coupling
          (A 750 1 or 2; met_drag_wave): the roughness length of Taylor
          and Yelland (2001) from Hs and the peak wave length,
            z0 = 1200 Hs (Hs/Lp)^4.5,
          Lp from the linear dispersion at the local water depth (Guo
          2002), with W 229 1 (default) not below the Charnock (1955)
          roughness 0.011 u*^2/g (Smith 1988; W 229 0: no floor),
          Cd = (kappa/ln(10/z0))^2 capped at W 228 (default 3.5e-3);
          where there are no waves (Hs < 0.1 m) law 5
W 207 0 (default): 1 for W 200 1/2 (as before), 5 for the vortex and the
          blend (W 200 3, 4).

Laws 1-4 are those of wind_f (A 570), evaluated per cell with the local
wind speed.
--------------------------------------------------------------------*/

// the law of W 207 (0: auto)
inline int met_drag_law(int W207, int W200) {return W207>0 ? W207 : (W200>=3 ? 5 : 1);}

inline double met_drag(int law, double U, double cap=2.5e-3, double floor=1.5e-3)
{
    double Cd = 0.0;

    if(law==1)
    Cd = (U<41.0) ? 0.001*(0.75 + 0.067*U) : 0.001*3.5;

    if(law==2)
    Cd = (U<7.5) ? 1.2875e-3 : 0.001*(0.8 + 0.065*U);

    if(law==3)
    Cd = 0.001*(0.63 + 0.066*U);

    if(law==4)
    {
    const double x = U/31.5;
    Cd = 0.001*(0.55 + 2.97*x - 1.49*x*x);
    Cd = Cd>0.0 ? Cd : 0.0;
    }

    // law 7 without a sea state: law 5
    if(law==5 || law==7)
    {
    Cd = 0.001*(0.75 + 0.067*U);
    Cd = Cd<cap ? Cd : cap;
    }

    if(law==6)
    {
    const double x = U/31.5;
    Cd = 0.001*(0.55 + 2.97*x - 1.49*x*x);

        if(x>2.97/2.98)
        Cd = Cd>floor ? Cd : floor;
    }

    return Cd;
}

// wave number [1/m] of the angular frequency w [1/s] at the depth d [m]: Guo (2002), explicit, within 0.75 %
inline double met_wavenumber(double w, double d, double g=9.81)
{
    if(!(d>0.0))
    return w*w/g;

    const double x = w*std::sqrt(d/g);
    const double kh = x*x*std::pow(1.0 - std::exp(-std::pow(x,2.5)),-0.4);

    return kh/d;
}

// law 7, the parts: the roughness length [m] of Taylor and Yelland (2001) of the sea state Hs [m], Tp [s] at the
// depth d [m] (-1: no waves, Hs < 0.1 m or no Tp), the Charnock roughness 0.011 u*^2/g of the wind U [m/s] (fixed point
// of u*^2 = Cd U^2), and the drag of a roughness length capped at cap
inline double met_z0_waves(double Hs, double Tp, double d)
{
    const double pi = 3.14159265358979323846;

    if(!(Hs>=0.1) || !(Tp>0.5))
    return -1.0;

    const double Lp = 2.0*pi/met_wavenumber(2.0*pi/Tp,d);

    return 1200.0*Hs*std::pow(Hs/Lp,4.5);
}

inline double met_z0_charnock(double U)
{
    const double g = 9.81, kappa = 0.4;
    double Cd = 1.2e-3;

    if(!(U>0.0))
    return 0.0;

    for(int n=0; n<30; ++n)
    {
    const double c = std::pow(kappa/std::log(10.0*g/(0.011*Cd*U*U)),2.0);
    const bool conv = std::fabs(c-Cd)<1.0e-12;
    Cd = c;

        if(conv)
        break;
    }

    return 0.011*Cd*U*U/g;
}

inline double met_drag_z0(double z0, double cap)
{
    z0 = std::min(std::max(z0,1.0e-7),1.0);
    const double Cd = std::pow(0.4/std::log(10.0/z0),2.0);

    return Cd<cap ? Cd : cap;
}

// law 7: drag of the 10 m wind U [m/s] over the sea state Hs [m], Tp [s] at the water depth d [m]; floor: Charnock
inline double met_drag_wave(double U, double Hs, double Tp, double d, double cap=3.5e-3, double cap5=2.5e-3, bool floor=true)
{
    const double z0 = met_z0_waves(Hs,Tp,d);

    if(z0<0.0 || !(U>0.0))
    return met_drag(5,U,cap5);

    return met_drag_z0(floor ? std::max(z0,met_z0_charnock(U)) : z0,cap);
}

#endif
