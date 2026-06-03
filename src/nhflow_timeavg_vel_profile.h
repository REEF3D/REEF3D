/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>
--------------------------------------------------------------------
Author: Thomas Becker
--------------------------------------------------------------------*/

#ifndef NHFLOW_TIMEAVG_VEL_PROFILE_H_
#define NHFLOW_TIMEAVG_VEL_PROFILE_H_

#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

class nhflow_timeavg_vel_profile : public increment
{
public:
    nhflow_timeavg_vel_profile(lexer*,fdm_nhf*);
    virtual ~nhflow_timeavg_vel_profile();

    void start(lexer*, fdm_nhf*, ghostcell*);

private:
    void ini_location(lexer*, fdm_nhf*);

    char name[100];
    int *iloc,*jloc,*flag;
    int n,q;
    const int probenum;
    ofstream *pout;

    double *time_accum,*start_time;
    bool *started,*printed;
    int max_points;

    double **z_timeint;
    double **u_timeint;
    double **v_timeint;
    double **vel_timeint;
};

#endif
