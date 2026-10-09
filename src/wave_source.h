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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef WAVE_SOURCE_H_
#define WAVE_SOURCE_H_

class lexer;
class wave_lib;

/*--------------------------------------------------------------------
wave_lexer_context: the lexer fields a wave_lib reads per wave or writes
back during construction (wT, wN, wHs, ...). The wave libraries keep
reading these from the lexer, so every source stores its own copy and
swaps it in while it is evaluated. The legacy source (B 92) owns the
copy that stays in the lexer between evaluations, so all other code
sees exactly the single-source state it saw before.
--------------------------------------------------------------------*/

struct wave_lexer_context
{
    // written by the wave libraries
    double wT,wC,wA,wk,wL,wwp,ww,wd,wTp,wHs,wH,wAs,phiin,wts,wte,wLp,ww_s,ww_e;
    int wN;

    // wave input that may differ per source
    int B84,B85,B86,B87,B91,B92,B93,B94,B130,B133,B136,B138,B138_1,B138_2,B139;
    double B87_1,B87_2,B88,B91_1,B91_2,B93_1,B93_2,B94_wdt,B105_1,B131,B132_s,B132_e,B134;

    void save(const lexer*);
    void load(lexer*) const;
};

/*--------------------------------------------------------------------
wave_source: one wave input. Wraps an unchanged wave_lib plus what the
global B 9x parameters cannot express per wave.

Frame: positions arrive in the generation frame of the legacy source
(iowave xgen/ygen; the global frame for B 105 0 0 0). A source has an
origin (x0,y0) and a direction rot (degrees) relative to that frame;
its wave_lib turns velocities into the global frame with its own B 105.
A source given in the global frame (B 505 id 1) is turned into this
frame when it is read (wave_field::read).
--------------------------------------------------------------------*/

class wave_source
{
public:
    wave_source(int, int);
    ~wave_source();

    bool active(const lexer*) const;     // inside [ts,te] and the library's own window
    double ramp(const lexer*) const;     // cosine ramp up after ts, down before te
    void local(double, double, double&, double&) const;

    int id;
    int type;                // B 92 code
    double H,T;              // height (Hs) and period (Tp)
    double rot;              // direction relative to the legacy frame [deg]
    double phase;            // phase lag [deg], applied as a time shift phase/360*T
    double ts,te,t_ramp;     // time window [s] and ramp duration [s]
    double x0,y0;            // origin in the legacy generation frame
    int seed;                // random phases of irregular sources (B 139, B 138)
    bool global = false;     // B 505 id 1: direction and origin given in the global frame
    double dir_in=0.0, x0_in=0.0, y0_in=0.0;   // as given (global frame), for the log

    double tshift;           // phase/360*T
    double cr,sr;            // cos/sin of rot

    wave_lib *lib;
    wave_lexer_context ctx;
};

#endif
